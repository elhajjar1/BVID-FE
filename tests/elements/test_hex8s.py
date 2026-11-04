"""Verification of the solid-shell element ``Hex8SolidShellElement``
(Hex8i incompatible modes + Dvorkin-Bathe ANS transverse shear), issue #132
PR1. Numbered as in the #132 design plan, sec. 5:

1. constant-strain patch tests (distorted element, 2x2x2 patch) and the
   rigid-body modes;
2. a 100:1 cantilever against beam theory;
3./4. pristine SSSS plate buckling against the exact cross-ply Navier value
   and the converged D16/D26 Ritz value for QI;
5. the Kirchhoff limit (G13, G23 x 50), where Hex8i alone locks;
6. h-refinement.

The 5 mm plate solves (34k DOFs, ~15 s each) are ``slow``.
"""

from __future__ import annotations

import functools
from typing import Literal

import numpy as np
import pytest

from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.elements.gauss import gauss_points_hex
from bvidfe.elements.hex8 import DegenerateElementError, Hex8Element
from bvidfe.elements.hex8i import Hex8iElement
from bvidfe.elements.hex8s import (
    Hex8SolidShellElement,
    build_solid_shell_table,
    geometric_stiffness_batch,
    solid_shell_stiffness_batch,
    solid_shell_stress_batch,
    stress_voigt_to_tensor,
)
from tests.elements._plate_buckling import (
    CP8,
    NAVIER_CP_MPA,
    NAVIER_QI_MPA,
    QI8,
    RITZ_QI_MPA,
    fe_ssss_buckling_stress,
    ritz_ssss_buckling_stress,
)

MAT = MATERIAL_LIBRARY["IM7/8552"]

_UNIT = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ],
    dtype=float,
)


def _box(dx: float, dy: float, dz: float, origin=(0.0, 0.0, 0.0)) -> np.ndarray:
    return _UNIT * np.array([dx, dy, dz]) + np.asarray(origin, dtype=float)


def _distorted() -> np.ndarray:
    """In-plane distorted element: a general quadrilateral extruded through a
    uniform thickness (the shape class of a layered plate mesh)."""
    quad = np.array([[0.0, 0.0], [4.3, 0.4], [3.6, 3.5], [-0.5, 2.7]])
    return np.vstack(
        [np.column_stack([quad, np.full(4, 0.2)]), np.column_stack([quad, np.full(4, 1.0)])]
    )


def _warped() -> np.ndarray:
    """Every node moved: thickness and through-thickness edges vary."""
    rng = np.random.default_rng(7)
    return _box(4.0, 3.0, 1.0) + rng.uniform(-0.3, 0.3, (8, 3))


def _linear_field(A: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nodal displacements of ``u = A x`` and its Voigt strain."""
    eps = np.array(
        [A[0, 0], A[1, 1], A[2, 2], A[1, 2] + A[2, 1], A[0, 2] + A[2, 0], A[0, 1] + A[1, 0]]
    )
    return (X @ A.T).ravel(), eps


# ---------------------------------------------------------------- 1. patch tests


@pytest.mark.parametrize("angle", [0.0, 30.0, -45.0, 90.0])
def test_constant_strain_reproduced_on_distorted_element(angle):
    X = _distorted()
    elem = Hex8SolidShellElement(X, MAT, angle)
    A = np.random.default_rng(3).standard_normal((3, 3))
    u, eps = _linear_field(A, X)
    # Strain at every Gauss point (enhanced part included) is the imposed one,
    # so the incompatible modes stay at zero and the ANS rows are exact.
    assert np.allclose(elem.strain_at_gauss_points(u), eps, rtol=0, atol=1e-12)
    sigma = elem.stress_at_gauss_points(u)
    assert np.allclose(sigma, elem._C_global @ eps, rtol=1e-12, atol=1e-9)
    # Nodal forces K u are the consistent forces of the constant stress.
    f_consistent = np.einsum("gij,i,g->j", elem.shell_table.B, sigma[0], elem.shell_table.dV)
    assert np.allclose(elem.stiffness_matrix() @ u, f_consistent, rtol=1e-10, atol=1e-8)


def _patch_mesh(distort: bool) -> tuple[np.ndarray, np.ndarray]:
    """2 x 2 x 2 element patch with uneven spacing; optionally the interior
    through-thickness node column moved in-plane (every element then an
    in-plane distorted prism)."""
    xs, ys, zs = np.array([0, 1.3, 3.0]), np.array([0, 0.8, 2.0]), np.array([0, 0.4, 1.0])
    nodes = np.array([[x, y, z] for z in zs for y in ys for x in xs])

    def nid(i, j, k):
        return i + 3 * j + 9 * k

    conn = []
    for k in range(2):
        for j in range(2):
            for i in range(2):
                conn.append(
                    [
                        nid(i, j, k),
                        nid(i + 1, j, k),
                        nid(i + 1, j + 1, k),
                        nid(i, j + 1, k),
                        nid(i, j, k + 1),
                        nid(i + 1, j, k + 1),
                        nid(i + 1, j + 1, k + 1),
                        nid(i, j + 1, k + 1),
                    ]
                )
    if distort:
        for k in range(3):
            nodes[nid(1, 1, k)] += [0.21, -0.13, 0.0]
    return nodes, np.array(conn)


@pytest.mark.parametrize("angle", [0.0, 30.0])
@pytest.mark.parametrize("distort", [False, True], ids=["cuboid", "distorted"])
def test_patch_test_2x2x2(distort, angle):
    """Boundary nodes on a linear field: the free interior node lands on it
    and every Gauss point of every element has the imposed strain (one
    material, so the constant stress is in equilibrium)."""
    nodes, conn = _patch_mesh(distort)
    elems = [Hex8SolidShellElement(nodes[c], MAT, angle) for c in conn]
    n = 3 * len(nodes)
    K = np.zeros((n, n))
    for c, e in zip(conn, elems):
        dofs = (3 * c[:, None] + np.arange(3)).ravel()
        K[np.ix_(dofs, dofs)] += e.stiffness_matrix()
    A = np.array([[1e-3, 2e-4, -3e-4], [5e-4, -7e-4, 1e-4], [-2e-4, 4e-4, 6e-4]])
    u_exact, eps = _linear_field(A, nodes)
    inner = 13  # node (1, 1, 1)
    free = np.arange(3 * inner, 3 * inner + 3)
    fixed = np.setdiff1d(np.arange(n), free)
    u = u_exact.copy()
    u[free] = np.linalg.solve(K[np.ix_(free, free)], -K[np.ix_(free, fixed)] @ u_exact[fixed])
    assert np.allclose(u[free], u_exact[free], rtol=0, atol=1e-13)
    for c, e in zip(conn, elems):
        dofs = (3 * c[:, None] + np.arange(3)).ravel()
        assert np.allclose(e.strain_at_gauss_points(u[dofs]), eps, rtol=0, atol=1e-12)


@pytest.mark.parametrize(
    "nodes",
    [_box(1.0, 1.0, 1.0), _box(10.0, 8.0, 1.0), _distorted(), _warped()],
    ids=["cube", "plate-like", "distorted", "warped"],
)
@pytest.mark.parametrize("angle", [0.0, 45.0])
def test_free_element_has_exactly_six_rigid_body_modes(nodes, angle):
    K = Hex8SolidShellElement(nodes, MAT, angle).stiffness_matrix()
    assert np.allclose(K, K.T, rtol=0, atol=1e-12 * np.abs(K).max())
    w = np.linalg.eigvalsh(K) / np.abs(K).max()
    assert np.all(np.abs(w[:6]) < 1e-12), w[:7]
    # no spurious zero-energy mode: the 7th eigenvalue is a real deformation
    assert w[6] > 1e-6, w[:7]
    # the rigid motions are in its null space
    c = nodes.mean(axis=0)
    for axis in range(3):
        t = np.zeros(3)
        t[axis] = 1.0
        assert np.allclose(K @ np.tile(t, 8), 0.0, atol=1e-10 * np.abs(K).max())
        rot = np.cross(t, nodes - c).ravel()
        assert np.allclose(K @ rot, 0.0, atol=1e-10 * np.abs(K).max())


def test_cuboid_reduces_to_hex8i_when_transverse_shear_is_constant():
    """On a cuboid the bubbles are Hex8i's: the two elements give the same
    strain energy for a field with no transverse shear (in-plane bending
    and membrane), and differ only through the tied transverse shear."""
    X = _box(10.0, 8.0, 0.152)
    K_s = Hex8SolidShellElement(X, MAT, 30.0).stiffness_matrix()
    K_i = Hex8iElement(X, MAT, 30.0).stiffness_matrix()
    A = np.array([[1e-3, 2e-4, 0.0], [5e-4, -7e-4, 0.0], [0.0, 0.0, 6e-4]])
    u, _ = _linear_field(A, X)
    assert u @ K_s @ u == pytest.approx(u @ K_i @ u, rel=1e-12)
    assert not np.allclose(K_s, K_i)


def test_pure_bending_exact_and_stress_includes_enhanced_strain():
    """Exact uniaxial pure bending of a 0-deg ply, u = -k x z, v = nu12 k y z,
    w = k/2 (x^2 - nu12 y^2 + nu13 z^2): every Gauss-point strain is exact,
    the strain energy is the beam's, and the recovered stress is
    sigma_xx = -E11 k z only. The thickness strain nu13 k z exists only
    through the enhanced (zeta-bubble) strain, so a recovery without it is
    wrong."""
    a, b, t, kap = 10.0, 8.0, 0.5, 1e-3
    X = _box(a, b, t, origin=(0.0, 0.0, -0.5 * t))
    x, y, z = X.T
    nu = MAT.nu12
    u = np.column_stack(
        [-kap * x * z, nu * kap * y * z, 0.5 * kap * (x**2 - nu * y**2 + nu * z**2)]
    )
    elem = Hex8SolidShellElement(X, MAT, 0.0)
    tbl = elem.shell_table
    gp, _ = gauss_points_hex(2)
    gp_z = 0.5 * t * gp[:, 2]
    eps_exact = np.zeros((8, 6))
    eps_exact[:, 0] = -kap * gp_z
    eps_exact[:, 1] = nu * kap * gp_z
    eps_exact[:, 2] = nu * kap * gp_z
    assert np.allclose(elem.strain_at_gauss_points(u.ravel()), eps_exact, rtol=0, atol=1e-15)
    sigma = elem.stress_at_gauss_points(u.ravel())
    assert np.allclose(sigma[:, 0], -MAT.E11 * kap * gp_z, rtol=1e-10)
    assert np.allclose(sigma[:, 1:], 0.0, atol=1e-9 * MAT.E11 * kap * t)
    energy = 0.5 * u.ravel() @ elem.stiffness_matrix() @ u.ravel()
    assert energy == pytest.approx(0.5 * MAT.E11 * kap**2 * a * b * t**3 / 12.0, rel=1e-10)
    # without the enhanced part the thickness strain is missing
    compatible = np.einsum("gij,j->gi", tbl.B, u.ravel())
    assert np.allclose(compatible[:, 2], 0.0, atol=1e-15)
    # plain Hex8 locks in this bending field
    hex8 = Hex8Element(X, MAT, 0.0).stiffness_matrix()
    assert 0.5 * u.ravel() @ hex8 @ u.ravel() > 2.0 * energy


# ------------------------------------------------- geometric stiffness and batches


def test_geometric_stiffness_uniform_stress_matches_hex8():
    X = _warped()
    S = np.array([[-2.0, 0.3, 0.1], [0.3, 0.5, -0.2], [0.1, -0.2, 0.05]])
    Kg = Hex8SolidShellElement(X, MAT, 45.0).geometric_stiffness_matrix(S)
    assert np.allclose(Kg, Hex8Element(X, MAT, 45.0).geometric_stiffness_matrix(S), atol=1e-12)
    assert np.allclose(
        Hex8SolidShellElement(X, MAT).geometric_stiffness_matrix(np.broadcast_to(S, (8, 3, 3))), Kg
    )


def test_geometric_stiffness_from_gauss_point_field():
    """Linear in the Gauss-point stress, symmetric, and each Gauss point
    contributes kron(gradN^T S gradN, I3) detJ w."""
    X = _warped()
    elem = Hex8SolidShellElement(X, MAT)
    rng = np.random.default_rng(11)
    S1 = stress_voigt_to_tensor(rng.standard_normal((8, 6)))
    S2 = stress_voigt_to_tensor(rng.standard_normal((8, 6)))
    K1, K2 = elem.geometric_stiffness_matrix(S1), elem.geometric_stiffness_matrix(S2)
    assert np.allclose(K1, K1.T, atol=1e-12)
    assert np.allclose(elem.geometric_stiffness_matrix(2.0 * S1 - S2), 2.0 * K1 - K2, atol=1e-11)
    one = np.zeros((8, 3, 3))
    one[5] = S1[5]
    tbl = elem.shell_table
    H = tbl.gradN[5].T @ S1[5] @ tbl.gradN[5] * tbl.dV[5]
    assert np.allclose(elem.geometric_stiffness_matrix(one), np.kron(H, np.eye(3)), atol=1e-12)
    with pytest.raises(ValueError):
        elem.geometric_stiffness_matrix(np.zeros((4, 3, 3)))


def test_voigt_to_tensor_ordering():
    S = stress_voigt_to_tensor(np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
    assert np.array_equal(S, [[1.0, 6.0, 5.0], [6.0, 2.0, 4.0], [5.0, 4.0, 3.0]])


def test_batched_builds_match_per_element():
    """Elements of one shape share a table and differ in C (a damaged mesh)."""
    X = _box(2.5, 2.5, 0.152)
    table = build_solid_shell_table(X)
    elems = [Hex8SolidShellElement(X, MAT, a, shell_table=table) for a in (0.0, 45.0, -45.0, 90.0)]
    C = np.stack([e._C_global for e in elems])
    C[2] *= 0.15  # a softened element
    K = solid_shell_stiffness_batch(table, C)
    U = np.random.default_rng(5).standard_normal((4, 24)) * 1e-3
    sig = solid_shell_stress_batch(table, C, U)
    for i, e in enumerate(elems):
        if i == 2:
            e._C_global = C[2]
        assert np.allclose(K[i], e.stiffness_matrix(), rtol=1e-13, atol=1e-9)
        assert np.allclose(sig[i], e.stress_at_gauss_points(U[i]), rtol=1e-12, atol=1e-12)
    Kg = geometric_stiffness_batch(table, stress_voigt_to_tensor(sig))
    for i, e in enumerate(elems):
        assert np.allclose(Kg[i], e.geometric_stiffness_matrix(stress_voigt_to_tensor(sig[i])))


def test_table_is_read_only_and_degenerate_element_rejected():
    table = build_solid_shell_table(_box(1.0, 1.0, 1.0))
    with pytest.raises(ValueError):
        table.B[0, 0, 0] = 1.0
    flat = _box(1.0, 1.0, 1.0)
    flat[4:, 2] = 0.0
    with pytest.raises(DegenerateElementError):
        Hex8SolidShellElement(flat, MAT)


# ---------------------------------------------------------- 2. thin cantilever


def _cantilever_tip_deflection(element_cls, L=100.0, b=10.0, n_el=10) -> tuple[float, float]:
    """[0/90]s strip, 1 element per ply (0.25 mm), clamped root, unit tip
    shear. Returns (FE mean tip w, beam theory with shear)."""
    layup, t_ply = [0.0, 90.0, 90.0, 0.0], 0.25
    xs = np.linspace(0.0, L, n_el + 1)
    ys = np.array([0.0, b])
    zs = t_ply * np.arange(len(layup) + 1)
    nx, ny = len(xs), len(ys)
    nodes = np.array([[x, y, z] for z in zs for y in ys for x in xs])

    def nid(i, j, k):
        return i + nx * (j + ny * k)

    n = 3 * len(nodes)
    K = np.zeros((n, n))
    for k, angle in enumerate(layup):
        for i in range(n_el):
            c = np.array(
                [
                    nid(i, 0, k),
                    nid(i + 1, 0, k),
                    nid(i + 1, 1, k),
                    nid(i, 1, k),
                    nid(i, 0, k + 1),
                    nid(i + 1, 0, k + 1),
                    nid(i + 1, 1, k + 1),
                    nid(i, 1, k + 1),
                ]
            )
            dofs = (3 * c[:, None] + np.arange(3)).ravel()
            K[np.ix_(dofs, dofs)] += element_cls(nodes[c], MAT, angle).stiffness_matrix()
    root = np.flatnonzero(np.abs(nodes[:, 0]) < 1e-9)
    tip = np.flatnonzero(np.abs(nodes[:, 0] - L) < 1e-9)
    f = np.zeros(n)
    # consistent nodal forces of a uniform shear traction on the tip face
    face_weight = np.where(np.isin(np.round(nodes[tip, 2] / t_ply), [0, len(layup)]), 0.5, 1.0) / (
        2 * len(layup)
    )
    f[3 * tip + 2] = -face_weight
    fixed = (3 * root[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(n), fixed)
    u = np.zeros(n)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], f[free])
    w_tip = -float(u[3 * tip + 2].mean())

    lam = Laminate(MAT, layup, t_ply)
    _, _, D = lam.abd_matrices()
    EI = b / np.linalg.inv(D)[0, 0]  # narrow beam: anticlastic curvature free
    h = lam.thickness_mm
    Gxz = np.mean([MAT.G13 if a == 0.0 else MAT.G23 for a in layup])
    beam = L**3 / (3.0 * EI) + L / (5.0 / 6.0 * Gxz * b * h)
    return w_tip, beam


def test_thin_cantilever_matches_beam_theory():
    """100:1 strip, 10 elements along it (element aspect 40:1)."""
    w, beam = _cantilever_tip_deflection(Hex8SolidShellElement)
    assert w == pytest.approx(beam, rel=0.01)
    # the locking regression: plain Hex8 is ~4x too stiff on the same mesh
    w8, _ = _cantilever_tip_deflection(Hex8Element)
    assert w8 < 0.3 * beam


# -------------------------------------------------- 3.-6. SSSS plate buckling


@functools.lru_cache(maxsize=None)
def _fe(
    layup: tuple, size_mm: float, element: Literal["hex8s", "hex8i"] = "hex8s", shear: float = 1.0
) -> float:
    """Cached: the slow h-refinement test reuses the 5 mm solve."""
    lam, _, _ = fe_ssss_buckling_stress(list(layup), size_mm, element, shear)
    return lam


def test_ritz_reference_values():
    """The references the plate tests use: the polynomial Ritz reproduces
    the exact cross-ply Navier value and is converged for QI with D16/D26."""
    assert ritz_ssss_buckling_stress(CP8, 8) == pytest.approx(NAVIER_CP_MPA, abs=5e-4)
    qi = [ritz_ssss_buckling_stress(QI8, m) for m in (8, 10, 12)]
    assert qi[0] >= qi[1] >= qi[2]
    assert qi[2] == pytest.approx(RITZ_QI_MPA, abs=1e-3)
    assert qi[2] < 0.975 * NAVIER_QI_MPA  # D16/D26 lower it by ~3%


def test_cross_ply_plate_10mm():
    assert _fe(tuple(CP8), 10.0) == pytest.approx(NAVIER_CP_MPA, rel=0.01)  # +0.5%


@pytest.mark.slow
def test_cross_ply_plate_5mm():
    assert _fe(tuple(CP8), 5.0) == pytest.approx(NAVIER_CP_MPA, rel=0.002)  # -0.05%


def test_quasi_isotropic_plate_10mm():
    lam = _fe(tuple(QI8), 10.0)
    assert lam == pytest.approx(RITZ_QI_MPA, rel=0.02)  # -0.3%
    assert lam < NAVIER_QI_MPA


@pytest.mark.slow
def test_quasi_isotropic_plate_5mm():
    lam = _fe(tuple(QI8), 5.0)
    assert lam == pytest.approx(RITZ_QI_MPA, rel=0.02)  # -1.0% (transverse-shear flexibility)
    assert lam < NAVIER_QI_MPA


def test_kirchhoff_limit_ans_removes_shear_locking_10mm():
    """G13, G23 x 50 (nu23 kept): the 3D plate approaches CLT. ANS stays
    within 1% of the Ritz value; Hex8i alone locks by ~+39% at 66:1."""
    assert _fe(tuple(QI8), 10.0, "hex8s", 50.0) == pytest.approx(RITZ_QI_MPA, rel=0.01)  # +0.4%
    assert _fe(tuple(QI8), 10.0, "hex8i", 50.0) > 1.3 * RITZ_QI_MPA


@pytest.mark.slow
def test_kirchhoff_limit_ans_5mm():
    assert _fe(tuple(QI8), 5.0, "hex8s", 50.0) == pytest.approx(RITZ_QI_MPA, rel=0.01)  # +0.06%


@pytest.mark.slow
def test_cross_ply_h_refinement():
    """20 / 10 / 5 mm: monotone from above, error shrinking (2.3% / 0.5% /
    0.05%), Richardson order >= 1.5."""
    f = [_fe(tuple(CP8), h) for h in (20.0, 10.0, 5.0)]
    err = [abs(v / NAVIER_CP_MPA - 1.0) for v in f]
    assert f[0] > f[1] > f[2]
    assert err[0] > err[1] > err[2]
    order = np.log2((f[0] - f[1]) / (f[1] - f[2]))
    assert order >= 1.5, order
