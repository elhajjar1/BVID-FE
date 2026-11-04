"""Solid-shell Hex8: incompatible modes plus assumed transverse shear (issue #132).

A buckling-capable 8-node brick for thin laminated plates, one element per
ply. It is the internal building block of the planned fe3d sublaminate
buckling path (issue #132); nothing in the fe3d tier uses it yet, and the
fe3d first-ply-failure solve keeps the plain :class:`~bvidfe.elements.hex8.Hex8Element`.

Formulation
-----------
* **Incompatible modes** (Wilson-Taylor bubbles, as in ``hex8i.py``):
  ``M_1 = 1 - xi^2``, ``M_2 = 1 - eta^2``, ``M_3 = 1 - zeta^2`` for each of
  the three displacement components (9 internal parameters ``alpha``),
  gradients mapped with the centre Jacobian ``J0`` and scaled by
  ``det J0 / det J`` (Taylor, Beresford & Wilson 1976; Simo & Rifai 1990), so
  the enhanced strains integrate to zero over any element shape. On a cuboid
  this is exactly ``Hex8iElement``. The bubbles carry in-plane bending and
  the Poisson thickness strain; ``alpha`` is condensed out per element.
* **Assumed natural transverse shear** (Dvorkin & Bathe 1984; Bathe &
  Dvorkin 1985 tying, applied at each ``zeta`` level; ``zeta`` must be the
  thickness direction): the covariant shear ``E_xi,zeta`` is sampled at
  ``(xi, eta) = (0, -1)`` and ``(0, +1)`` and interpolated linearly in
  ``eta``; ``E_eta,zeta`` is sampled at ``(-1, 0)`` and ``(+1, 0)`` and
  interpolated linearly in ``xi``. Covariant components are in the local
  basis at each point (``E = J eps J^T``), so the element is invariant to
  rigid rotations and keeps exactly six zero-energy (rigid-body) modes. On
  the axis-aligned cuboids of the structured fe3d mesh they are the
  Cartesian ``gamma_xz`` / ``gamma_yz`` up to a scale. This removes the
  transverse-shear locking that Hex8i alone shows when the ply's transverse
  shear is stiff (the Kirchhoff limit).
* **Patch test.** Constant strain is reproduced exactly, and a multi-element
  patch passes, for elements whose thickness edges are parallel to one
  direction with a uniform thickness (any in-plane distortion: every element
  of a layered plate mesh). Like MITC shells, the tied shear does not
  reproduce constant strain on a warped element (thickness varying, or
  through-thickness edges not parallel).
* Transverse-normal (thickness) ANS is not used: thickness locking needs
  curvature, and the laminate mesh is flat.
* The **geometric stiffness** is the standard compatible form
  ``K_g = sum_gp kron(gradN^T S_gp gradN, I_3) detJ w`` (Bathe sec. 6.8) with
  a full 3x3 Cauchy stress ``S_gp`` given at each of the 8 Gauss points, so a
  prestress from an actual static solution can be used.
* **Stress recovery** includes the enhanced strain
  (``alpha = -Kaa^-1 Kau u``), without which the bending stress of an
  unsymmetric sublaminate's prestress would be wrong.

Verification (``tests/elements/test_hex8s.py``): constant-strain patch tests
on in-plane distorted elements and 2x2x2 patches, exactly 6 rigid-body zero modes, a 100:1
cantilever, and SSSS plate buckling against the exact cross-ply Navier value
and a converged D16/D26 polynomial Ritz, including transverse shear stiffened
50x.

The batched functions (``solid_shell_stiffness_batch`` and friends) work on a
:class:`SolidShellTable` shared by every element of one shape (on the
structured mesh: one per ply thickness) with a per-element constitutive
matrix, which is how a mesh with per-element damage would be built.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np

from bvidfe.core.material import OrthotropicMaterial
from bvidfe.elements.gauss import gauss_points_hex
from bvidfe.elements.hex8 import _b_at_point, _shape_derivatives, _validate_jacobian
from bvidfe.elements.hex8i import Hex8iElement

# Number of 2x2x2 Gauss points.
_N_GP = 8

# Elements per chunk in the batched builds: bounds the (chunk, 8, 6, 24)
# temporaries to ~20 MB.
_BATCH_CHUNK = 2048

# Voigt index pairs, engineering shear: [xx, yy, zz, 2yz, 2xz, 2xy].
_VOIGT_PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))


@dataclass(frozen=True)
class SolidShellTable:
    """Per-Gauss-point operators for one element shape (2x2x2 rule).

    Attributes
    ----------
    B : np.ndarray
        (8, 6, 24) assumed-strain strain-displacement matrices (on a
        cuboid: the compatible membrane and normal rows and the tied
        transverse-shear rows).
    G : np.ndarray
        (8, 6, 9) enhanced (incompatible-mode) strain operators, already
        scaled by ``det J0 / det J``.
    dV : np.ndarray
        (8,) Gauss weight times ``det J``.
    gradN : np.ndarray
        (8, 3, 8) compatible shape-function spatial gradients, for ``K_g``.
    """

    B: np.ndarray
    G: np.ndarray
    dV: np.ndarray
    gradN: np.ndarray


def _strain_to_natural(J: np.ndarray) -> np.ndarray:
    """6x6 Voigt map ``T`` with ``E_voigt = T eps_voigt``, ``E = J eps J^T``.

    ``J[i, j] = dx_j / dxi_i``, so ``E_ij = g_i . eps . g_j`` with ``g_i``
    the rows of ``J`` (covariant components in the local basis). Both
    sides use engineering shear.
    """
    T = np.empty((6, 6))
    for k, (a, b) in enumerate(_VOIGT_PAIRS):
        eps = np.zeros((3, 3))
        val = 1.0 if a == b else 0.5
        eps[a, b] = val
        eps[b, a] = val
        E = J @ eps @ J.T
        for r, (i, j) in enumerate(_VOIGT_PAIRS):
            T[r, k] = E[i, j] if i == j else 2.0 * E[i, j]
    return T


def _bubble_operator(xi: float, eta: float, zeta: float, J0_inv: np.ndarray) -> np.ndarray:
    """Unscaled enhanced strain operator (6, 9), bubble gradients via ``J0``.

    Same layout as ``Hex8iElement._Benh``: columns ``3k + d`` are bubble k
    (xi, eta, zeta) in displacement direction d (x, y, z).
    """
    dMn = np.zeros((3, 3))
    dMn[0, 0] = -2.0 * xi
    dMn[1, 1] = -2.0 * eta
    dMn[2, 2] = -2.0 * zeta
    dM = J0_inv @ dMn  # (3, 3): rows x/y/z, columns bubble index
    G = np.zeros((6, 9))
    for k in range(3):
        Mx, My, Mz = dM[0, k], dM[1, k], dM[2, k]
        c = 3 * k
        G[0, c + 0] = Mx
        G[1, c + 1] = My
        G[2, c + 2] = Mz
        G[3, c + 1] = Mz
        G[3, c + 2] = My
        G[4, c + 0] = Mz
        G[4, c + 2] = Mx
        G[5, c + 0] = My
        G[5, c + 1] = Mx
    return G


def _natural_b(node_coords: np.ndarray, xi: float, eta: float, zeta: float) -> np.ndarray:
    """Covariant strain-displacement rows ``T(J) B`` (6, 24) at a point, in
    the local covariant basis there (``E = J eps J^T``)."""
    dN = _shape_derivatives(xi, eta, zeta)
    B, _, _ = _b_at_point(node_coords, dN, xi, eta, zeta)
    return _strain_to_natural(dN @ node_coords) @ B


def build_solid_shell_table(node_coords: np.ndarray) -> SolidShellTable:
    """Build the :class:`SolidShellTable` for one element's node coordinates.

    The thickness direction must be ``zeta`` (nodes 0-3 one face, 4-7 the
    other, as in the fe3d mesh), since the tied components are the
    ``xi-zeta`` and ``eta-zeta`` shears.

    Raises :class:`~bvidfe.elements.hex8.DegenerateElementError` if ``det J``
    is non-positive at the centre, a Gauss point or a tying point.
    """
    X = np.asarray(node_coords, dtype=float)
    if X.shape != (8, 3):
        raise ValueError(f"node_coords must be (8,3), got {X.shape}")
    gp, wt = gauss_points_hex(order=2)
    J0 = _shape_derivatives(0.0, 0.0, 0.0) @ X
    detJ0 = float(np.linalg.det(J0))
    _validate_jacobian(detJ0, 0.0, 0.0, 0.0, X)
    J0_inv = np.linalg.inv(J0)

    B_all = np.empty((_N_GP, 6, 24))
    G_all = np.empty((_N_GP, 6, 9))
    dV = np.empty(_N_GP)
    gradN_all = np.empty((_N_GP, 3, 8))
    tied: dict[tuple[float, float, float], np.ndarray] = {}

    def tied_b(xi: float, eta: float, zeta: float) -> np.ndarray:
        key = (xi, eta, zeta)
        if key not in tied:
            tied[key] = _natural_b(X, xi, eta, zeta)
        return tied[key]

    for ig in range(_N_GP):
        xi, eta, zeta = (float(v) for v in gp[ig])
        dN = _shape_derivatives(xi, eta, zeta)
        B, detJ, gradN = _b_at_point(X, dN, xi, eta, zeta)
        T = _strain_to_natural(dN @ X)
        Bn = T @ B
        # E_eta,zeta (row 3): tied at (xi, eta) = (+-1, 0), linear in xi.
        Bn[3] = 0.5 * (1.0 - xi) * tied_b(-1.0, 0.0, zeta)[3]
        Bn[3] += 0.5 * (1.0 + xi) * tied_b(1.0, 0.0, zeta)[3]
        # E_xi,zeta (row 4): tied at (xi, eta) = (0, +-1), linear in eta.
        Bn[4] = 0.5 * (1.0 - eta) * tied_b(0.0, -1.0, zeta)[4]
        Bn[4] += 0.5 * (1.0 + eta) * tied_b(0.0, 1.0, zeta)[4]
        B_all[ig] = np.linalg.solve(T, Bn)
        G_all[ig] = (detJ0 / detJ) * _bubble_operator(xi, eta, zeta, J0_inv)
        dV[ig] = detJ * wt[ig]
        gradN_all[ig] = gradN
    for arr in (B_all, G_all, dV, gradN_all):
        # Tables are shared across elements (and cached): keep them immutable.
        arr.setflags(write=False)
    return SolidShellTable(B=B_all, G=G_all, dV=dV, gradN=gradN_all)


def _as_stack(C: np.ndarray) -> np.ndarray:
    C = np.asarray(C, dtype=float)
    if C.ndim != 3 or C.shape[1:] != (6, 6):
        raise ValueError(f"C must be (n, 6, 6), got {C.shape}")
    return C


def _condensation_blocks(
    table: SolidShellTable, C: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(Kuu (n,24,24), Kua (n,24,9), Kaa (n,9,9))`` for a stack of ``C``."""
    CB = np.einsum("nij,gjk->ngik", C, table.B)  # (n, 8, 6, 24)
    CG = np.einsum("nij,gjk->ngik", C, table.G)  # (n, 8, 6, 9)
    Kuu = np.einsum("gji,ngjk,g->nik", table.B, CB, table.dV)
    Kua = np.einsum("gji,ngjk,g->nik", table.B, CG, table.dV)
    Kaa = np.einsum("gji,ngjk,g->nik", table.G, CG, table.dV)
    return Kuu, Kua, Kaa


def solid_shell_stiffness_batch(table: SolidShellTable, C: np.ndarray) -> np.ndarray:
    """Condensed 24x24 stiffness of ``n`` elements sharing ``table``.

    Parameters
    ----------
    table : the elements' common :class:`SolidShellTable`.
    C : (n, 6, 6) global-frame constitutive matrices, one per element.

    Returns
    -------
    (n, 24, 24) symmetric element stiffness matrices.
    """
    C = _as_stack(C)
    out = np.empty((C.shape[0], 24, 24))
    for s in range(0, C.shape[0], _BATCH_CHUNK):
        Kuu, Kua, Kaa = _condensation_blocks(table, C[s : s + _BATCH_CHUNK])
        K = Kuu - Kua @ np.linalg.solve(Kaa, np.swapaxes(Kua, 1, 2))
        out[s : s + _BATCH_CHUNK] = 0.5 * (K + np.swapaxes(K, 1, 2))
    return out


def solid_shell_strain_batch(table: SolidShellTable, C: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Total (assumed + enhanced) Voigt strain (n, 8, 6) at the Gauss points.

    ``u`` is (n, 24) element displacements; ``C`` (n, 6, 6) is needed to
    recover the condensed parameters ``alpha = -Kaa^-1 Kau u``.
    """
    C = _as_stack(C)
    u = np.asarray(u, dtype=float)
    if u.shape != (C.shape[0], 24):
        raise ValueError(f"u must be ({C.shape[0]}, 24), got {u.shape}")
    out = np.empty((C.shape[0], _N_GP, 6))
    for s in range(0, C.shape[0], _BATCH_CHUNK):
        sl = slice(s, s + _BATCH_CHUNK)
        _, Kua, Kaa = _condensation_blocks(table, C[sl])
        rhs = np.einsum("nia,ni->na", Kua, u[sl])[..., None]
        alpha = -np.linalg.solve(Kaa, rhs)[..., 0]
        out[sl] = np.einsum("gij,nj->ngi", table.B, u[sl]) + np.einsum(
            "gij,nj->ngi", table.G, alpha
        )
    return out


def solid_shell_stress_batch(table: SolidShellTable, C: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Global-frame Voigt stress (n, 8, 6) at the Gauss points, enhanced
    strain included (see :func:`solid_shell_strain_batch`)."""
    C = _as_stack(C)
    eps = solid_shell_strain_batch(table, C, u)
    return np.einsum("nij,ngj->ngi", C, eps)


def stress_voigt_to_tensor(sigma: np.ndarray) -> np.ndarray:
    """Voigt stress ``[..., 6]`` ``[sxx, syy, szz, tyz, txz, txy]`` to
    symmetric tensors ``[..., 3, 3]``."""
    s = np.asarray(sigma, dtype=float)
    if s.shape[-1] != 6:
        raise ValueError(f"Voigt stress must have last axis 6, got {s.shape}")
    S = np.empty(s.shape[:-1] + (3, 3))
    for r, (i, j) in enumerate(_VOIGT_PAIRS):
        S[..., i, j] = s[..., r]
        S[..., j, i] = s[..., r]
    return S


def geometric_stiffness_batch(table: SolidShellTable, sigma_gp: np.ndarray) -> np.ndarray:
    """Geometric stiffness (n, 24, 24) from a Gauss-point stress field.

    ``K_g = sum_gp kron(gradN^T S_gp gradN, I_3) detJ w`` with ``sigma_gp``
    (n, 8, 3, 3) the global-frame Cauchy stress at each element's 8 Gauss
    points (MPa). Symmetric, linear in the stress, and possibly indefinite.
    With a uniform stress it equals ``Hex8Element.geometric_stiffness_matrix``.
    The buckling multiplier ``lambda`` solves ``K phi = -lambda K_g phi``
    for a compressive prestress, i.e. ``(K + lambda K_g) phi = 0``.
    """
    S = np.asarray(sigma_gp, dtype=float)
    if S.ndim != 4 or S.shape[1:] != (_N_GP, 3, 3):
        raise ValueError(f"sigma_gp must be (n, 8, 3, 3), got {S.shape}")
    H = np.einsum("gai,ngab,gbj,g->nij", table.gradN, S, table.gradN, table.dV)  # (n, 8, 8)
    n = H.shape[0]
    return np.einsum("nij,de->nidje", H, np.eye(3)).reshape(n, 24, 24)


@functools.lru_cache(maxsize=256)
def _cached_table(key: bytes) -> SolidShellTable:
    return build_solid_shell_table(np.frombuffer(key, dtype=float).reshape(8, 3))


class Hex8SolidShellElement(Hex8iElement):
    """Hex8 with incompatible modes and ANS transverse shear (solid shell).

    Drop-in 24-DOF element with the node ordering and material handling of
    :class:`~bvidfe.elements.hex8.Hex8Element`. Its stiffness, strain and
    stress recovery use the assumed-strain operators of
    :class:`SolidShellTable`; :meth:`geometric_stiffness_matrix` accepts a
    uniform (3, 3) stress or a per-Gauss-point (8, 3, 3) field.

    Parameters
    ----------
    node_coords, material, ply_angle_deg : as ``Hex8Element``.
    shell_table : optional precomputed :class:`SolidShellTable` shared by
        elements of the same shape; built (and cached) from ``node_coords``
        when omitted.
    """

    def __init__(
        self,
        node_coords: np.ndarray,
        material: OrthotropicMaterial,
        ply_angle_deg: float = 0.0,
        shell_table: SolidShellTable | None = None,
    ) -> None:
        super().__init__(node_coords, material, ply_angle_deg)
        if shell_table is None:
            shell_table = _cached_table(np.ascontiguousarray(self.node_coords).tobytes())
        self.shell_table = shell_table

    def stiffness_matrix(self) -> np.ndarray:
        """Condensed 24x24 stiffness (incompatible modes eliminated)."""
        return solid_shell_stiffness_batch(self.shell_table, self._C_global[None])[0]

    def strain_at_gauss_points(self, u_elem: np.ndarray) -> np.ndarray:
        """Total Voigt strain (8, 6) at the Gauss points, enhanced part included."""
        u = np.asarray(u_elem, dtype=float).reshape(1, 24)
        return solid_shell_strain_batch(self.shell_table, self._C_global[None], u)[0]

    def stress_at_gauss_points(self, u_elem: np.ndarray, C: np.ndarray | None = None) -> np.ndarray:
        """Global-frame Voigt stress (8, 6) at the Gauss points, enhanced
        strain included. ``C`` replaces the element's constitutive matrix in
        the stress (``alpha`` is always condensed with the element's own)."""
        eps = self.strain_at_gauss_points(u_elem)
        C_use = self._C_global if C is None else np.asarray(C, dtype=float)
        return eps @ C_use.T

    def geometric_stiffness_matrix(self, sigma_bar_3x3: np.ndarray) -> np.ndarray:
        """Geometric stiffness 24x24 for a uniform (3, 3) stress or a
        per-Gauss-point (8, 3, 3) stress field (global frame, MPa).

        See :func:`geometric_stiffness_batch`. Uses the compatible gradients,
        so a uniform stress gives exactly ``Hex8Element``'s ``K_g``.
        """
        S = np.asarray(sigma_bar_3x3, dtype=float)
        if S.shape == (3, 3):
            S = np.broadcast_to(S, (_N_GP, 3, 3))
        elif S.shape != (_N_GP, 3, 3):
            raise ValueError(f"sigma must be (3,3) or ({_N_GP},3,3), got {S.shape}")
        return geometric_stiffness_batch(self.shell_table, S[None])[0]
