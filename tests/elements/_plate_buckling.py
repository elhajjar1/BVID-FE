"""Pristine SSSS plate buckling on the structured fe3d mesh, for verifying
the solid-shell element (issue #132, PR1).

Not a production path: nothing in ``bvidfe`` calls this. It assembles the
elastic stiffness of the chosen element and the geometric stiffness of a
uniform ``sigma_xx = -1 MPa`` prestress on ``build_fe_mesh_skeleton``'s mesh
(one element per ply), restrains ``u_z`` on the mid-surface node lines of the
four edges (simply supported), ``u_x`` on the ``x = 0`` mid-surface line and
``u_y`` at one node, and solves ``K phi = lambda (-K_g) phi`` with
``linear_buckling``. ``lambda`` is then the far-field buckling stress in MPa.

The reference is a Rayleigh-Ritz solution of the CLT plate with D16/D26,
``w = x (a - x) y (b - y) P_i(x) P_j(y)`` (Legendre polynomials; only the
geometric boundary conditions are imposed, so it converges from above to
the SSSS value). For a specially orthotropic plate it reproduces the Navier
value.
"""

from __future__ import annotations

import dataclasses
import time
from typing import Literal

import numpy as np
from numpy.polynomial import legendre as leg
from numpy.polynomial import polynomial as poly
from scipy.linalg import eigh

from bvidfe.analysis import AnalysisConfig, MeshParams
from bvidfe.analysis.fe_mesh import build_fe_mesh_skeleton
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY, OrthotropicMaterial
from bvidfe.damage.state import DamageState
from bvidfe.elements.hex8i import Hex8iElement
from bvidfe.elements.hex8s import (
    build_solid_shell_table,
    geometric_stiffness_batch,
    solid_shell_stiffness_batch,
)
from bvidfe.solver.assembler import assemble_coo
from bvidfe.solver.boundary import mid_surface_z
from bvidfe.solver.buckling import linear_buckling

QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
CP8 = [0, 90, 90, 0, 0, 90, 90, 0]
PLY_T = 0.152
LX, LY = 200.0, 150.0

# Exact Navier value for the canonical cross-ply panel (D16 = D26 = 0) and the
# converged D16/D26 Ritz value for the QI panel (#132 design plan, sec. 3.1).
NAVIER_CP_MPA = 11.583
RITZ_QI_MPA = 12.204
NAVIER_QI_MPA = 12.582  # D16 = D26 = 0: the closed form #129 compared against


def material(shear_factor: float = 1.0) -> OrthotropicMaterial:
    """IM7/8552, with G13 and G23 scaled by ``shear_factor`` (nu23 kept)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    if shear_factor == 1.0:
        return m
    return dataclasses.replace(
        m, G13=m.G13 * shear_factor, G23=m.G23 * shear_factor, nu23=m.nu23_resolved
    )


def _ritz_1d(a: float, M: int, n_gauss: int = 64):
    """Integrals ``I(r, s)[i, j] = int f_i^(r) f_j^(s) dx`` of the basis
    ``f_i = x (a - x) P_i(2x/a - 1)`` on ``[0, a]``."""
    xg, wg = leg.leggauss(n_gauss)
    x = 0.5 * (xg + 1.0) * a
    w = 0.5 * wg * a
    F = np.empty((M, 3, n_gauss))
    for i in range(M):
        c = leg.leg2poly(np.eye(M)[i])
        px: np.ndarray = np.zeros(1)
        for k, ck in enumerate(c):
            px = poly.polyadd(px, ck * poly.polypow([-1.0, 2.0 / a], k))
        px = poly.polymul(px, [0.0, a, -1.0])
        for r in range(3):
            F[i, r] = poly.polyval(x, poly.polyder(px, r) if r else px)

    def integral(r: int, s: int) -> np.ndarray:
        return np.einsum("ik,jk,k->ij", F[:, r], F[:, s], w)

    return integral


def ritz_ssss_buckling_stress(layup: list[float], M: int = 10) -> float:
    """Uniaxial (N_x) SSSS buckling stress (MPa) of the CLT plate over
    ``LX x LY``, D16/D26 included."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], layup, PLY_T)
    _, _, D = lam.abd_matrices()
    Ix, Iy = _ritz_1d(LX, M), _ritz_1d(LY, M)
    k = np.kron
    D11, D22, D12, D66, D16, D26 = D[0, 0], D[1, 1], D[0, 1], D[2, 2], D[0, 2], D[1, 2]
    K = (
        D11 * k(Ix(2, 2), Iy(0, 0))
        + D22 * k(Ix(0, 0), Iy(2, 2))
        + D12 * (k(Ix(2, 0), Iy(0, 2)) + k(Ix(0, 2), Iy(2, 0)))
        + 4.0 * D66 * k(Ix(1, 1), Iy(1, 1))
        + 2.0 * D16 * (k(Ix(2, 1), Iy(0, 1)) + k(Ix(1, 2), Iy(1, 0)))
        + 2.0 * D26 * (k(Ix(0, 1), Iy(2, 1)) + k(Ix(1, 0), Iy(1, 2)))
    )
    G = k(Ix(1, 1), Iy(0, 0))
    return float(eigh(K, G, eigvals_only=True)[0]) / lam.thickness_mm


def fe_ssss_buckling_stress(
    layup: list[float],
    size_mm: float,
    element: Literal["hex8s", "hex8i"] = "hex8s",
    shear_factor: float = 1.0,
) -> tuple[float, int, float]:
    """Lowest FE buckling stress (MPa) of the pristine SSSS panel.

    Returns ``(lambda_1, n_dof, seconds)``.
    """
    t0 = time.perf_counter()
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=layup,
        ply_thickness_mm=PLY_T,
        panel=PanelGeometry(LX, LY),
        loading="compression",
        tier="fe3d",
        damage=DamageState([], 0.0),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=size_mm),
    )
    mesh = build_fe_mesh_skeleton(cfg)
    mat = material(shear_factor)
    X = mesh.node_coords
    conn = mesh.element_connectivity

    # Elements of one shape and ply angle share a stiffness (pristine plate).
    ext = X[conn].max(axis=1) - X[conn].min(axis=1)
    keys = np.column_stack([np.round(ext, 9), mesh.ply_angles_deg])
    uniq, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = inverse.ravel()
    sigma = np.zeros((1, 8, 3, 3))
    sigma[..., 0, 0] = -1.0
    Ke = np.empty((len(uniq), 24, 24))
    Kge = np.empty((len(uniq), 24, 24))
    for g in range(len(uniq)):
        e = int(np.flatnonzero(inverse == g)[0])
        nodes = X[conn[e]]
        angle = float(mesh.ply_angles_deg[e])
        ref = Hex8iElement(nodes, mat, angle)
        table = build_solid_shell_table(nodes)
        if element == "hex8s":
            Ke[g] = solid_shell_stiffness_batch(table, ref._C_global[None])[0]
        else:
            Ke[g] = ref.stiffness_matrix()
        Kge[g] = geometric_stiffness_batch(table, sigma)[0]

    K = assemble_coo(mesh.element_dof_maps, 3 * len(X), lambda e: Ke[inverse[e]])
    Kg = assemble_coo(mesh.element_dof_maps, 3 * len(X), lambda e: Kge[inverse[e]])

    tol = 1e-9
    mid = np.abs(X[:, 2] - mid_surface_z(X)) < tol
    x0 = np.abs(X[:, 0]) < tol
    edge = (
        x0 | (np.abs(X[:, 0] - LX) < tol) | (np.abs(X[:, 1]) < tol) | (np.abs(X[:, 1] - LY) < tol)
    )
    fixed = set((3 * np.flatnonzero(mid & edge) + 2).tolist())
    fixed |= set((3 * np.flatnonzero(mid & x0)).tolist())
    line = np.flatnonzero(mid & x0)
    fixed.add(3 * int(line[np.argmin(np.abs(X[line, 1] - 0.5 * LY))]) + 1)
    free = np.setdiff1d(np.arange(3 * len(X)), np.fromiter(fixed, dtype=int))
    K_ff = K[free][:, free]
    G_ff = (-Kg)[free][:, free]
    eigs, _ = linear_buckling(K_ff, G_ff, n_modes=2)
    lam1 = float(np.min(eigs[eigs > 0]))
    return lam1, 3 * len(X), time.perf_counter() - t0
