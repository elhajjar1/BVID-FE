"""3D finite-element tier for BVID-FE.

v0.1.0 pragmatic simplification: first-ply-failure-under-damaged-mesh.
True buckling-based CAI is deferred to v0.2.0.
"""

from __future__ import annotations

from typing import List

from bvidfe.analysis.config import AnalysisConfig
from bvidfe.analysis.fe_mesh import FeMesh, build_fe_mesh
from bvidfe.core.laminate import Laminate
from bvidfe.damage.state import DamageState
from bvidfe.elements.hex8 import Hex8Element
from bvidfe.failure.larc05 import larc05_index
from bvidfe.failure.tsai_wu import tsai_wu_index
from bvidfe.solver.boundary import compression_bcs, tension_bcs
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


def fe3d_cai(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
) -> float:
    """3D FE compression-after-impact residual strength (MPa)."""
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
