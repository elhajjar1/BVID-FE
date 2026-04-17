"""3D finite-element tier for BVID-FE.

v0.2.0: Primary CAI path uses geometric-stiffness-based linear buckling.
First-ply-failure on damaged mesh is retained as a fallback / comparison.
"""

from __future__ import annotations

from typing import List

import numpy as np
import scipy.sparse as sp

from bvidfe.analysis.config import AnalysisConfig
from bvidfe.analysis.fe_mesh import FeMesh, build_fe_mesh
from bvidfe.core.laminate import Laminate
from bvidfe.damage.state import DamageState
from bvidfe.elements.hex8 import Hex8Element
from bvidfe.failure.larc05 import larc05_index
from bvidfe.failure.tsai_wu import tsai_wu_index
from bvidfe.solver.assembler import assemble_global_stiffness
from bvidfe.solver.boundary import (
    BoundaryCondition,
    apply_dirichlet_penalty,
    compression_bcs,
    tension_bcs,
)
from bvidfe.solver.buckling import linear_buckling
from bvidfe.solver.static import solve_linear_static


def _build_elements(mesh: FeMesh, lam: Laminate) -> List[Hex8Element]:
    """Build a Hex8Element for each mesh element, with damage-aware stiffness scaling."""
    elements: List[Hex8Element] = []
    for eidx in range(mesh.n_elements):
        node_ids = mesh.element_connectivity[eidx]
        node_coords = mesh.node_coords[node_ids]
        mat = lam.material
        ply_angle = float(mesh.ply_angles_deg[eidx])
        elem = Hex8Element(node_coords, mat, ply_angle_deg=ply_angle)
        # Scale the pre-computed global stiffness by damage factor
        factor = float(mesh.damage_factors[eidx])
        elem._C_global = elem._C_global * factor
        elements.append(elem)
    return elements


def _max_failure_index_at_strain(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    applied_strain: float,
    criterion: str,
) -> float:
    """Solve static FE at given applied strain; return max failure index across all Gauss points."""
    if criterion == "larc05":
        bcs = compression_bcs(mesh.node_coords, applied_strain)
    else:
        bcs = tension_bcs(mesh.node_coords, applied_strain)
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs)
    max_idx = 0.0
    material = cfg.material if not isinstance(cfg.material, str) else None
    # Need the actual OrthotropicMaterial for the failure criterion
    if material is None:
        from bvidfe.core.material import MATERIAL_LIBRARY

        material = MATERIAL_LIBRARY[cfg.material]

    for eidx, elem in enumerate(elements):
        dof_map = mesh.element_dof_maps[eidx]
        u_elem = u[dof_map]
        stress_field = elem.stress_at_gauss_points(u_elem)  # (n_gp, 6)
        for gp in range(stress_field.shape[0]):
            stress = stress_field[gp]
            if criterion == "larc05":
                idx = larc05_index(material, stress)
            else:
                idx = tsai_wu_index(material, stress)
            if idx > max_idx:
                max_idx = idx
    return max_idx


def _bisect_failure_strain(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    strain_sign: int,
    criterion: str,
    strain_lo: float = 1e-5,
    strain_hi: float = 0.05,
    max_iter: int = 20,
    tol: float = 0.02,
) -> float:
    """Bisect on |applied_strain| until max failure index ≈ 1."""
    idx_hi = _max_failure_index_at_strain(cfg, mesh, elements, strain_sign * strain_hi, criterion)
    if idx_hi < 1.0:
        return strain_hi  # even at strain_hi we don't fail; return upper bracket
    lo, hi = strain_lo, strain_hi
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        idx_mid = _max_failure_index_at_strain(cfg, mesh, elements, strain_sign * mid, criterion)
        if abs(idx_mid - 1.0) < tol:
            return mid
        if idx_mid > 1.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _effective_modulus(lam: Laminate) -> float:
    """Effective in-plane Young's modulus along x (Ex from CLT)."""
    Ex, _, _, _ = lam.effective_engineering_constants()
    return Ex


def _fe3d_cai_first_ply_failure(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
) -> float:
    """3D FE compression-after-impact residual strength via first-ply-failure (MPa).

    Original v0.1.0 implementation — bisects on applied strain until LaRC05 failure
    index reaches 1 on the damaged mesh. Retained as fallback / comparison path.
    """
    mesh = build_fe_mesh(cfg, damage)
    elements = _build_elements(mesh, lam)
    strain_at_failure = _bisect_failure_strain(
        cfg,
        mesh,
        elements,
        strain_sign=-1,
        criterion="larc05",
    )
    E = _effective_modulus(lam)
    sigma = strain_at_failure * E
    return min(sigma, sigma_pristine_MPa)


def fe3d_cai_buckling(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
    sigma_ref_MPa: float = 1.0,
) -> tuple[float, float]:
    """3D FE compression-after-impact via true linear buckling eigensolve.

    Assembles K and K_g under a constant uniaxial pre-stress sigma_ref along x
    (scaled by per-element damage factor), then solves K phi = lambda K_g phi for
    the smallest positive eigenvalue. Critical buckling stress = lambda * sigma_ref.

    Uses the "constant pre-stress" approximation (Cook §17.7 / Bathe §6.8):
    - sigma_0 = sigma_ref_MPa along x everywhere (unit reference compression).
    - Damaged elements carry damage_factor * sigma_0 (reduced stress-carrying capacity).
    - K_g is assembled element-by-element from Hex8Element.geometric_stiffness_matrix().
    - A minimal penalty-BC set suppresses rigid-body modes before the eigensolve.

    Returns
    -------
    (sigma_critical_MPa, lambda_crit)
        sigma_critical_MPa : min(lambda_crit * sigma_ref, sigma_pristine_MPa)
        lambda_crit        : smallest positive buckling load factor (0 if solve failed)
    """
    mesh = build_fe_mesh(cfg, damage)
    elements = _build_elements(mesh, lam)

    # Assemble elastic stiffness K
    K = assemble_global_stiffness(elements, mesh.element_dof_maps, mesh.n_dof)

    # Assemble geometric stiffness K_g under uniform uniaxial pre-stress sigma_ref along x.
    # Damaged elements carry damage_factor * sigma_ref (reduced load fraction).
    sigma_bar_ref = np.zeros((3, 3))
    sigma_bar_ref[0, 0] = sigma_ref_MPa

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    for eidx, elem in enumerate(elements):
        damage_f = float(mesh.damage_factors[eidx])
        sigma_elem = sigma_bar_ref * damage_f
        Kg_e = elem.geometric_stiffness_matrix(sigma_elem)
        dof_map = mesh.element_dof_maps[eidx]
        for i in range(24):
            gi = int(dof_map[i])
            for j in range(24):
                rows.append(gi)
                cols.append(int(dof_map[j]))
                data.append(Kg_e[i, j])

    Kg = sp.coo_matrix((data, (rows, cols)), shape=(K.shape[0], K.shape[0])).tocsc()

    # Apply penalty BCs to K (not Kg) to suppress rigid-body modes.
    # Fix node 0 fully (DOFs 0,1,2), plus translational DOFs on two
    # neighbour nodes to remove all six rigid-body modes.
    n_dof = K.shape[0]
    bcs = [
        BoundaryCondition(dof=0, value=0.0),
        BoundaryCondition(dof=1, value=0.0),
        BoundaryCondition(dof=2, value=0.0),
        BoundaryCondition(dof=4, value=0.0),
        BoundaryCondition(dof=5, value=0.0),
        BoundaryCondition(dof=7, value=0.0),
    ]
    F_dummy = np.zeros(n_dof)
    K_mod, _ = apply_dirichlet_penalty(K, F_dummy, bcs, penalty=1.0e10)

    # Solve generalised eigenproblem K phi = lambda K_g phi for smallest positive eig.
    try:
        n_req = min(6, n_dof - 1)
        eigs, _ = linear_buckling(K_mod, Kg, n_modes=n_req)
        positive_eigs = [float(e) for e in eigs if e > 1e-6]
        if not positive_eigs:
            return sigma_pristine_MPa, 0.0
        lambda_crit = min(positive_eigs)
    except Exception:
        return sigma_pristine_MPa, 0.0

    sigma_critical = lambda_crit * sigma_ref_MPa
    return min(sigma_critical, sigma_pristine_MPa), lambda_crit


def fe3d_cai(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
) -> float:
    """3D FE compression-after-impact residual strength (MPa).

    Primary path: true linear buckling eigensolve (fe3d_cai_buckling).
    Fallback: first-ply-failure on damaged mesh (_fe3d_cai_first_ply_failure).
    Returns the smaller of the two — whichever failure mode governs.
    """
    sigma_buckling, lambda_crit = fe3d_cai_buckling(cfg, damage, lam, sigma_pristine_MPa)
    sigma_fpf = _fe3d_cai_first_ply_failure(cfg, damage, lam, sigma_pristine_MPa)
    return min(sigma_buckling, sigma_fpf)


def fe3d_tai(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
) -> float:
    """3D FE tension-after-impact residual strength (MPa)."""
    mesh = build_fe_mesh(cfg, damage)
    elements = _build_elements(mesh, lam)
    strain_at_failure = _bisect_failure_strain(
        cfg,
        mesh,
        elements,
        strain_sign=+1,
        criterion="tsai_wu",
    )
    E = _effective_modulus(lam)
    sigma = strain_at_failure * E
    return min(sigma, sigma_pristine_MPa)
