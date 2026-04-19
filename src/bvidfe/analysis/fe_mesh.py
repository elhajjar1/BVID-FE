"""Structured hex mesh builder for the 3D FE tier.

Generates a regular brick mesh of the full panel with per-element ply
assignment and damage-aware stiffness reduction. The stiffness-reduction
approach (a scalar factor inside delamination footprints) is a simplification
of zero-thickness cohesive surfaces and is adequate for linear buckling
in v0.1.0. True cohesive surfaces are deferred.

DAMAGE_STIFFNESS_FACTOR = 0.3 represents typical residual in-plane
stiffness after delamination — the plies themselves are intact (so in-plane
load-carrying is mostly preserved), only the through-thickness coupling is
lost. Literature values for this "effective residual layer modulus"
fraction range from 0.1 to 0.5 for CFRP (see e.g. Bolotin 2001 review).
The prior value of 1e-4 was physically unrealistic — it treated damaged
elements as essentially null in the stress field, which meant the
failure-index criterion could never flag damaged regions, and the fe3d
residual-strength prediction actually *increased* with impact energy past
the point where damage exceeded ~15% of the panel area (because the peak
stress in undamaged elements drops as the damage footprint grows and
spreads the load over a larger periphery).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

import numpy as np

from bvidfe.analysis.config import AnalysisConfig, MeshParams
from bvidfe.damage.state import DamageState, DelaminationEllipse

DAMAGE_STIFFNESS_FACTOR = 0.30


def estimate_fe_mesh_size(config: AnalysisConfig) -> dict:
    """Return a dict with n_elements, n_nodes, n_dof for a would-be fe_mesh build.

    Use this to warn the user BEFORE running the analysis."""
    mesh = config.mesh if config.mesh is not None else MeshParams()
    n_plies = len(config.layup_deg)
    nx = max(1, math.ceil(config.panel.Lx_mm / mesh.in_plane_size_mm))
    ny = max(1, math.ceil(config.panel.Ly_mm / mesh.in_plane_size_mm))
    nz = n_plies * mesh.elements_per_ply
    n_elements = nx * ny * nz
    n_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    n_dof = 3 * n_nodes
    return {
        "n_elements": n_elements,
        "n_nodes": n_nodes,
        "n_dof": n_dof,
        "nx": nx,
        "ny": ny,
        "nz": nz,
    }


@dataclass
class FeMesh:
    """Structured hex mesh with per-element ply/damage metadata."""

    node_coords: np.ndarray  # (n_nodes, 3)
    element_connectivity: np.ndarray  # (n_elements, 8) node indices
    element_dof_maps: List[np.ndarray]  # length n_elements, each (24,)
    ply_indices: np.ndarray  # (n_elements,) int
    ply_angles_deg: np.ndarray  # (n_elements,) float
    damage_factors: np.ndarray  # (n_elements,) float in (0, 1]
    n_nodes: int = field(init=False)
    n_elements: int = field(init=False)
    n_dof: int = field(init=False)

    def __post_init__(self) -> None:
        self.n_nodes = self.node_coords.shape[0]
        self.n_elements = self.element_connectivity.shape[0]
        self.n_dof = 3 * self.n_nodes


def _point_in_ellipse(x: float, y: float, ellipse: DelaminationEllipse) -> bool:
    """Return True if (x, y) is inside the ellipse footprint."""
    c = math.cos(math.radians(-ellipse.orientation_deg))
    s = math.sin(math.radians(-ellipse.orientation_deg))
    dx = x - ellipse.centroid_mm[0]
    dy = y - ellipse.centroid_mm[1]
    xr = c * dx - s * dy
    yr = s * dx + c * dy
    return (xr / ellipse.major_mm) ** 2 + (yr / ellipse.minor_mm) ** 2 <= 1.0


def build_fe_mesh(config: AnalysisConfig, damage: DamageState) -> FeMesh:
    """Build a structured brick mesh for the panel with per-element metadata."""
    panel = config.panel
    layup = config.layup_deg
    h = config.ply_thickness_mm
    n_plies = len(layup)
    mesh = config.mesh if config.mesh is not None else MeshParams()

    Lx, Ly, Lz = panel.Lx_mm, panel.Ly_mm, n_plies * h
    nx = max(1, math.ceil(Lx / mesh.in_plane_size_mm))
    ny = max(1, math.ceil(Ly / mesh.in_plane_size_mm))
    nz = n_plies * mesh.elements_per_ply

    # Nodes on a regular grid
    x_nodes = np.linspace(0.0, Lx, nx + 1)
    y_nodes = np.linspace(0.0, Ly, ny + 1)
    z_nodes = np.linspace(0.0, Lz, nz + 1)
    node_coords = np.array(
        [[x, y, z] for z in z_nodes for y in y_nodes for x in x_nodes],
        dtype=float,
    )

    def _node_id(i: int, j: int, k: int) -> int:
        """Flat index into node_coords. i along x, j along y, k along z."""
        return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i

    # Build element connectivity (Abaqus hex convention)
    connectivity = np.zeros((nx * ny * nz, 8), dtype=int)
    ply_indices = np.zeros(nx * ny * nz, dtype=int)
    ply_angles = np.zeros(nx * ny * nz, dtype=float)
    damage_factors = np.ones(nx * ny * nz, dtype=float)
    element_dof_maps: List[np.ndarray] = []

    elem_idx = 0
    for k in range(nz):
        ply_i = k // mesh.elements_per_ply
        for j in range(ny):
            for i in range(nx):
                nodes_this_element = [
                    _node_id(i, j, k),
                    _node_id(i + 1, j, k),
                    _node_id(i + 1, j + 1, k),
                    _node_id(i, j + 1, k),
                    _node_id(i, j, k + 1),
                    _node_id(i + 1, j, k + 1),
                    _node_id(i + 1, j + 1, k + 1),
                    _node_id(i, j + 1, k + 1),
                ]
                connectivity[elem_idx] = nodes_this_element
                ply_indices[elem_idx] = ply_i
                ply_angles[elem_idx] = layup[ply_i]
                dof_map = np.array(
                    [3 * n + d for n in nodes_this_element for d in range(3)],
                    dtype=int,
                )
                element_dof_maps.append(dof_map)

                # Compute element centroid for damage check
                cx = 0.5 * (x_nodes[i] + x_nodes[i + 1])
                cy = 0.5 * (y_nodes[j] + y_nodes[j + 1])
                cz_top = z_nodes[k + 1]
                cz_bot = z_nodes[k]

                # Check delamination overlap: element straddles an interface at
                # z = (iface + 1) * h if cz_bot < z_iface < cz_top AND
                # (cx, cy) is inside the ellipse.
                for ell in damage.delaminations:
                    z_iface = (ell.interface_index + 1) * h
                    if cz_bot <= z_iface <= cz_top:
                        if _point_in_ellipse(cx, cy, ell):
                            damage_factors[elem_idx] = DAMAGE_STIFFNESS_FACTOR
                            break

                # Fiber-break core: any element within fiber_break_radius of
                # any delamination centroid (all through-thickness layers)
                if damage.fiber_break_radius_mm > 0:
                    for ell in damage.delaminations:
                        dx = cx - ell.centroid_mm[0]
                        dy = cy - ell.centroid_mm[1]
                        if math.sqrt(dx * dx + dy * dy) <= damage.fiber_break_radius_mm:
                            damage_factors[elem_idx] = DAMAGE_STIFFNESS_FACTOR
                            break

                elem_idx += 1

    return FeMesh(
        node_coords=node_coords,
        element_connectivity=connectivity,
        element_dof_maps=element_dof_maps,
        ply_indices=ply_indices,
        ply_angles_deg=ply_angles,
        damage_factors=damage_factors,
    )
