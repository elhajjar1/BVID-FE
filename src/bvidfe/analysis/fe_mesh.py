"""Structured hex mesh builder for the 3D FE tier.

Generates a regular brick mesh of the full panel with per-element ply
assignment and damage-aware stiffness reduction. The stiffness-reduction
approach (per-element factors inside delamination footprints) is a
simplification of zero-thickness cohesive surfaces, used for the
first-ply-failure solve. True cohesive surfaces are not implemented.

Damage is represented by TWO per-element factors applied component-wise
to the 6x6 elasticity matrix in the global Voigt frame
[xx, yy, zz, yz, xz, xy] = indices [0, 1, 2, 3, 4, 5]:

  - `damage_factors[e]`            — out-of-plane (OOP) factor in (0, 1].
                                     Scales the OOP block (rows/cols 2,
                                     3, 4) and the in-plane / OOP
                                     Poisson cross-coupling. Reduced to
                                     DAMAGE_OOP_FACTOR in elements that
                                     span a delaminated interface with
                                     their centroid inside that
                                     delamination's footprint, and in
                                     the fiber-break core.

  - `in_plane_damage_factors[e]`   — in-plane factor in (0, 1]. Scales
                                     the in-plane sub-block on rows/cols
                                     {0, 1, 5}. Reduced to
                                     DAMAGE_ZONE_INPLANE_FACTOR through
                                     the thickness inside the union of
                                     the delamination footprints (the
                                     damage zone), interpolated by the
                                     fraction of the element's area
                                     inside it, and to
                                     DAMAGE_FIBER_BREAK_INPLANE_FACTOR
                                     inside the fiber-break core (within
                                     `fiber_break_radius_mm` of any
                                     delamination centroid).

  - `fiber_tension_factors[e]`     — factor in (0, 1] on the fibre-
                                     direction modulus E11 (ply frame)
                                     when fe3d loads the panel in
                                     tension. 1 in the damage zone,
                                     whose delaminations and matrix
                                     cracks cut no fibres, and
                                     DAMAGE_FIBER_BREAK_INPLANE_FACTOR in
                                     the fiber-break core. Compression
                                     ignores it: there the whole in-plane
                                     sub-block takes the in-plane factor.

The two factors separate the interlaminar and in-plane effects of impact
damage. A delaminated interface loses its bond, so through-thickness
coupling (E33, G13, G23, and the E1-E3 / E2-E3 Poisson terms) drops
sharply there (O'Brien 1982; Pavier & Clarke 1995). The damage zone as a
whole also carries matrix cracks and fibre damage in every ply, which
lowers its in-plane stiffness: it is modelled as a soft inclusion over
the projected damage area (Soutis & Curtis 1996), the same footprint the
empirical tier's Soutis model uses. Without it the in-plane stress field
under in-plane load was the same with and without damage, and fe3d
first-ply failure could not see the impact. Inside the fiber-break core
under the impact site, fibers are broken in addition to interlaminar
separation (Camanho & Davila 2007; Maimi et al. 2007).

The zone's in-plane factor depends on the sign of the fibre load. In
compression it scales every in-plane term: delaminated sublaminates
buckle and matrix-cracked plies lose the lateral support that keeps the
fibres from kinking, and the factor is fitted to measured CAI. In tension
the fibres of the zone are intact and carry load as before, so only its
matrix-dominated moduli E22 and G12 take the factor (applied to the ply
compliance, Poisson terms kept) and E11 keeps its pristine value, as in
continuum damage models that give fibre tension its own damage variable,
set only by fibre failure (Matzenmiller, Lubliner & Taylor 1995; Maimi et
al. 2007), and in ply-discount analysis (E11 kept, E22 and G12
discounted). See fe_tier._build_elements.

The previous unified `DAMAGE_STIFFNESS_FACTOR = 0.30` (applied uniformly
to every entry of `_C_global` in `fe_tier._build_elements`) over-penalised
in-plane stiffness in delaminated zones — its docstring already noted
that "the plies themselves are intact (so in-plane load-carrying is
mostly preserved)", but the code did not honour that. The two-factor
model below makes the model match the docstring intent.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

import numpy as np

from bvidfe.analysis.config import AnalysisConfig, MeshParams
from bvidfe.damage.state import DamageState, DelaminationEllipse

# Out-of-plane stiffness fraction at delaminated interfaces. Representative
# of the post-delamination through-thickness / transverse-shear modulus
# loss reported in Bolotin (2001) review and Sun & Tao (1998).
DAMAGE_OOP_FACTOR = 0.05

# In-plane stiffness fraction inside the fiber-break core (matrix-crushed,
# fibers broken). Representative of the fiber-direction damage-saturation
# values in Camanho & Davila 2007 / Maimi et al. 2007 for CFRP.
DAMAGE_FIBER_BREAK_INPLANE_FACTOR = 0.30

# In-plane stiffness fraction of the damage zone (soft inclusion over the
# projected damage area, every ply). Its strengths drop with it (strain
# equivalence, see fe_tier._failure_stress_material). Under tension it
# scales only the matrix-dominated terms: the fibre-direction stiffness,
# and so the fibre tension strength, stay intact (fiber_tension_factors).
# Fitted to three
# measured CAI coupons with the delaminations kept inside the ASTM D7136
# impact window: the first NCAMP AS4/8552 [45/0/-45/90]3s coupon, the first
# Lovejoy & Scotti IM7/8552 [45/0/-45/90]4s coupon and their 40-ply hard
# laminate. Error of the fe3d strength against the measured CAI (5 mm mesh):
#   0.13: -26.7%, -6.8%, -7.9%     0.15: -20.6%, +1.2%, +0.3%
#   0.17: -14.6%, +9.0%, +7.9%     0.20: -5.6%, +20.6%, +11.1%
# 0.15 has the lowest mean error (7.4%). No single factor fits both
# families: NCAMP wants about 0.22, Lovejoy & Scotti about 0.15. With the
# old template (ply mismatch not wrapped modulo 180 deg) 0.17 fitted all
# three within 5%.
DAMAGE_ZONE_INPLANE_FACTOR = 0.15

# Sub-sampling grid (per side) for the fraction of an element's area inside
# the damage zone. A binary centroid test moved whole elements in or out of
# the zone as the mesh changed, a staircase that alone shifted the cross-ply
# knockdown by ~25% between 5 and 2.5 mm meshes.
ZONE_FRACTION_SAMPLES = 6


def estimate_fe_mesh_size(config: AnalysisConfig) -> dict:
    """Return a dict with n_elements, n_nodes, n_dof for a would-be fe_mesh build.

    Also ``nx``, ``ny`` and ``nz`` (elements along x, y and through the
    thickness). Cheap (no mesh is built): use it to check a config against
    ``fe_tier.FE3D_MAX_DOF`` BEFORE running the analysis."""
    nz, ny, nx = element_grid_shape(config)
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
    damage_factors: np.ndarray  # (n_elements,) out-of-plane factor in (0, 1]
    in_plane_damage_factors: np.ndarray  # (n_elements,) in-plane factor in (0, 1]
    # (n_elements,) factor on the ply-frame fibre-direction stiffness under
    # tension, in (0, 1]: 1 except in the fiber-break core
    fiber_tension_factors: np.ndarray
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


def _points_in_ellipse(x: np.ndarray, y: np.ndarray, ellipse: DelaminationEllipse) -> np.ndarray:
    """Vectorised ``_point_in_ellipse`` over coordinate arrays.

    cos/sin of the orientation are evaluated once (not once per element).
    The per-point arithmetic is the exact same sequence of scalar
    operations as ``_point_in_ellipse``, applied elementwise by numpy, so
    the boolean result is identical element by element.
    """
    c = math.cos(math.radians(-ellipse.orientation_deg))
    s = math.sin(math.radians(-ellipse.orientation_deg))
    dx = x - ellipse.centroid_mm[0]
    dy = y - ellipse.centroid_mm[1]
    xr = c * dx - s * dy
    yr = s * dx + c * dy
    return (xr / ellipse.major_mm) ** 2 + (yr / ellipse.minor_mm) ** 2 <= 1.0


@dataclass
class FeMeshSkeleton:
    """Damage-independent part of a structured FE mesh.

    Built once per configuration by ``build_fe_mesh_skeleton`` and reused
    across every iteration of a parametric sweep, where only the damage
    state changes. ``apply_damage`` fills in the per-element damage factors
    cheaply and returns a full :class:`FeMesh`.
    """

    node_coords: np.ndarray
    element_connectivity: np.ndarray
    element_dof_maps: List[np.ndarray]
    ply_indices: np.ndarray
    ply_angles_deg: np.ndarray
    centroid_x: np.ndarray  # (n_elements,) element centroid x
    centroid_y: np.ndarray  # (n_elements,) element centroid y
    cz_bot: np.ndarray  # (n_elements,) element bottom z
    cz_top: np.ndarray  # (n_elements,) element top z
    # Cumulative z position of every ply boundary, length n_plies + 1.
    # ply_top_z[k] is the bottom face of ply k; ply_top_z[k+1] its top
    # face (== the z of interface k). Carries per-ply thickness so
    # non-uniform laminates lay out the through-thickness grid correctly.
    ply_top_z: np.ndarray


def element_grid_shape(config: AnalysisConfig) -> tuple[int, int, int]:
    """``(nz, ny, nx)`` of the structured mesh: element ``e`` sits at
    ``e = k*(ny*nx) + j*nx + i``, so a per-element array reshapes to it."""
    mesh = config.mesh if config.mesh is not None else MeshParams()
    nx = max(1, math.ceil(config.panel.Lx_mm / mesh.in_plane_size_mm))
    ny = max(1, math.ceil(config.panel.Ly_mm / mesh.in_plane_size_mm))
    return len(config.layup_deg) * mesh.elements_per_ply, ny, nx


def build_fe_mesh_skeleton(config: AnalysisConfig) -> FeMeshSkeleton:
    """Build the damage-independent skeleton of the structured brick mesh.

    Connectivity and DOF maps are produced with pure numpy index arithmetic
    (no per-element Python loop), matching the element ordering
    ``elem = k*(ny*nx) + j*nx + i`` and the Abaqus hex node convention used
    by the original triple-loop builder exactly.
    """
    panel = config.panel
    layup = config.layup_deg
    n_plies = len(layup)
    mesh = config.mesh if config.mesh is not None else MeshParams()

    # Resolve per-ply thicknesses: a scalar applies uniformly, a sequence
    # gives one thickness per ply (length must equal n_plies).
    raw_t = config.ply_thickness_mm
    if isinstance(raw_t, (list, tuple, np.ndarray)):
        ply_thicknesses = [float(t) for t in raw_t]
    else:
        ply_thicknesses = [float(raw_t)] * n_plies
    # Cumulative ply-boundary z positions (length n_plies + 1).
    ply_top_z = np.empty(n_plies + 1, dtype=float)
    ply_top_z[0] = 0.0
    for k, t in enumerate(ply_thicknesses):
        ply_top_z[k + 1] = ply_top_z[k] + t

    Lx, Ly = panel.Lx_mm, panel.Ly_mm
    nz, ny, nx = element_grid_shape(config)

    # Nodes on a regular x-y grid; the z-grid follows per-ply thicknesses so
    # ply boundaries (and therefore interfaces) align exactly with element
    # faces. Each ply is subdivided into ``elements_per_ply`` equal-height
    # elements, so element height varies between plies but is constant
    # within a ply. For a uniform laminate this reproduces the previous
    # ``np.linspace(0, n_plies*h, nz+1)`` grid exactly.
    x_nodes = np.linspace(0.0, Lx, nx + 1)
    y_nodes = np.linspace(0.0, Ly, ny + 1)
    z_segments = [np.array([0.0])]
    for ply_i in range(n_plies):
        z_segments.append(
            np.linspace(
                ply_top_z[ply_i],
                ply_top_z[ply_i + 1],
                mesh.elements_per_ply + 1,
            )[1:]
        )
    z_nodes = np.concatenate(z_segments)
    node_coords = np.array(
        [[x, y, z] for z in z_nodes for y in y_nodes for x in x_nodes],
        dtype=float,
    )

    nxp1 = nx + 1
    nyp1 = ny + 1
    layer = nxp1 * nyp1  # nodes per k-layer

    # Element (i, j, k) grid in the original loop order: k outer, j, i inner.
    kk, jj, ii = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    ii = ii.ravel()
    jj = jj.ravel()
    kk = kk.ravel()  # length n_elements, ordered elem = k*(ny*nx) + j*nx + i

    def _nid(i: np.ndarray, j: np.ndarray, k: np.ndarray) -> np.ndarray:
        return k * layer + j * nxp1 + i

    # Connectivity (Abaqus hex convention), same 8-node ordering as before.
    connectivity = np.empty((ii.shape[0], 8), dtype=int)
    connectivity[:, 0] = _nid(ii, jj, kk)
    connectivity[:, 1] = _nid(ii + 1, jj, kk)
    connectivity[:, 2] = _nid(ii + 1, jj + 1, kk)
    connectivity[:, 3] = _nid(ii, jj + 1, kk)
    connectivity[:, 4] = _nid(ii, jj, kk + 1)
    connectivity[:, 5] = _nid(ii + 1, jj, kk + 1)
    connectivity[:, 6] = _nid(ii + 1, jj + 1, kk + 1)
    connectivity[:, 7] = _nid(ii, jj + 1, kk + 1)

    ply_idx_all = kk // mesh.elements_per_ply
    ply_indices = ply_idx_all.astype(int)
    layup_arr = np.asarray(layup, dtype=float)
    ply_angles = layup_arr[ply_idx_all]

    # dof_map[e] = [3*n + d for n in connectivity[e] for d in (0,1,2)]
    dof_maps_arr = (3 * connectivity[:, :, None] + np.arange(3)[None, None, :]).reshape(
        connectivity.shape[0], 24
    )
    element_dof_maps: List[np.ndarray] = [np.array(row, dtype=int) for row in dof_maps_arr]

    # Element centroids / through-thickness extents (identical floats to
    # the original 0.5*(x_nodes[i]+x_nodes[i+1]) etc.).
    cx_grid = 0.5 * (x_nodes[:-1] + x_nodes[1:])  # (nx,)
    cy_grid = 0.5 * (y_nodes[:-1] + y_nodes[1:])  # (ny,)
    centroid_x = cx_grid[ii]
    centroid_y = cy_grid[jj]
    cz_bot = z_nodes[kk]
    cz_top = z_nodes[kk + 1]

    return FeMeshSkeleton(
        node_coords=node_coords,
        element_connectivity=connectivity,
        element_dof_maps=element_dof_maps,
        ply_indices=ply_indices,
        ply_angles_deg=ply_angles,
        centroid_x=centroid_x,
        centroid_y=centroid_y,
        cz_bot=cz_bot,
        cz_top=cz_top,
        ply_top_z=ply_top_z,
    )


def _zone_area_fraction(skeleton: FeMeshSkeleton, damage: DamageState) -> np.ndarray:
    """Fraction of each element's x-y area inside the union of the
    delamination footprints, sampled on a ``ZONE_FRACTION_SAMPLES``-square
    grid per element. The zone runs through the thickness, so it is
    evaluated on the bottom layer of elements and repeated up the stack."""
    n_elem = skeleton.element_connectivity.shape[0]
    n_layer = int(np.count_nonzero(skeleton.cz_bot == skeleton.cz_bot.min()))
    # Nodes 0 and 2 are the (x_min, y_min) and (x_max, y_max) corners.
    lo = skeleton.node_coords[skeleton.element_connectivity[:n_layer, 0], :2]
    hi = skeleton.node_coords[skeleton.element_connectivity[:n_layer, 2], :2]
    t = (np.arange(ZONE_FRACTION_SAMPLES) + 0.5) / ZONE_FRACTION_SAMPLES
    inside = np.zeros(n_layer)
    for a in t:
        for b in t:
            px = lo[:, 0] + a * (hi[:, 0] - lo[:, 0])
            py = lo[:, 1] + b * (hi[:, 1] - lo[:, 1])
            hit = np.zeros(n_layer, dtype=bool)
            for ell in damage.delaminations:
                hit |= _points_in_ellipse(px, py, ell)
            inside += hit
    return (inside / t.size**2)[np.arange(n_elem) % n_layer]


def fibre_damage_hole_elements(
    mesh: FeMeshSkeleton | FeMesh, damage: DamageState, radius_mm: float
) -> np.ndarray:
    """Boolean mask of the elements whose x-y footprint overlaps a disk of
    ``radius_mm`` about any delamination centroid, in every ply: the
    elements under fe3d TAI's fibre-damage hole (see
    ``fe_tier._tai_fibre_damage_notch``). None for ``radius_mm <= 0``."""
    # Nodes 0 and 2 are the (x_min, y_min) and (x_max, y_max) corners.
    lo = mesh.node_coords[mesh.element_connectivity[:, 0], :2]
    hi = mesh.node_coords[mesh.element_connectivity[:, 2], :2]
    hit = np.zeros(lo.shape[0], dtype=bool)
    if radius_mm <= 0:
        return hit
    for ell in damage.delaminations:
        c = np.asarray(ell.centroid_mm, dtype=float)
        # Distance from the centroid to the nearest point of each element.
        gap = np.clip(c, lo, hi) - c
        hit |= np.hypot(gap[:, 0], gap[:, 1]) < radius_mm
    return hit


def apply_damage(skeleton: FeMeshSkeleton, damage: DamageState) -> FeMesh:
    """Fill per-element damage factors onto a skeleton, returning a FeMesh.

    Out-of-plane reduction where an element straddles a delaminated
    interface and its centroid lies in the ellipse; in-plane reduction of
    every element overlapping any ellipse (the damage zone, all plies),
    interpolated by the fraction of its area inside the zone; then
    fiber-break-core reduction (both factors, and the fibre-direction
    tension factor) within ``fiber_break_radius_mm`` of any delamination
    centroid.
    """
    n_elem = skeleton.element_connectivity.shape[0]
    cx = skeleton.centroid_x
    cy = skeleton.centroid_y
    cz_bot = skeleton.cz_bot
    cz_top = skeleton.cz_top
    ply_top_z = skeleton.ply_top_z

    damage_factors = np.ones(n_elem, dtype=float)
    in_plane_damage_factors = np.ones(n_elem, dtype=float)
    # Delaminations and matrix cracks cut no fibres: under tension the
    # zone keeps its fibre-direction stiffness. Only broken fibres lower it.
    fiber_tension_factors = np.ones(n_elem, dtype=float)

    # OOP: first ellipse (in order) whose interface the element straddles
    # and whose footprint contains the centroid sets the OOP factor. The
    # scalar code breaks on the first match; assigning a constant makes the
    # vectorised "any match" result identical.
    oop_hit = np.zeros(n_elem, dtype=bool)
    for ell in damage.delaminations:
        inside = _points_in_ellipse(cx, cy, ell)
        z_iface = ply_top_z[ell.interface_index + 1]
        straddles = (cz_bot <= z_iface) & (z_iface <= cz_top)
        oop_hit |= straddles & inside
    damage_factors[oop_hit] = DAMAGE_OOP_FACTOR
    # Damage zone: soft inclusion over the delamination footprints,
    # through the thickness. An element partly inside takes the
    # area-weighted factor.
    if damage.delaminations:
        frac = _zone_area_fraction(skeleton, damage)
        in_plane_damage_factors = 1.0 - frac * (1.0 - DAMAGE_ZONE_INPLANE_FACTOR)

    # Fiber-break core: within fiber_break_radius of any delamination
    # centroid → both factors reduced (OOP also forced to DAMAGE_OOP_FACTOR).
    if damage.fiber_break_radius_mm > 0:
        fb_hit = np.zeros(n_elem, dtype=bool)
        r = damage.fiber_break_radius_mm
        for ell in damage.delaminations:
            dx = cx - ell.centroid_mm[0]
            dy = cy - ell.centroid_mm[1]
            fb_hit |= np.sqrt(dx * dx + dy * dy) <= r
        damage_factors[fb_hit] = DAMAGE_OOP_FACTOR
        in_plane_damage_factors[fb_hit] = DAMAGE_FIBER_BREAK_INPLANE_FACTOR
        fiber_tension_factors[fb_hit] = DAMAGE_FIBER_BREAK_INPLANE_FACTOR

    return FeMesh(
        node_coords=skeleton.node_coords,
        element_connectivity=skeleton.element_connectivity,
        element_dof_maps=skeleton.element_dof_maps,
        ply_indices=skeleton.ply_indices,
        ply_angles_deg=skeleton.ply_angles_deg,
        damage_factors=damage_factors,
        in_plane_damage_factors=in_plane_damage_factors,
        fiber_tension_factors=fiber_tension_factors,
    )


def _skeleton_signature(config: AnalysisConfig) -> tuple:
    """Hashable signature of the mesh-defining (damage-independent) inputs.

    Two configs with equal signatures produce an identical skeleton, so the
    skeleton can be reused (e.g. across a ``sweep_energies`` run where only
    the damage state changes between iterations).
    """
    mesh = config.mesh if config.mesh is not None else MeshParams()
    raw_t = config.ply_thickness_mm
    if isinstance(raw_t, (list, tuple, np.ndarray)):
        t_sig: tuple = tuple(float(t) for t in raw_t)
    else:
        t_sig = (float(raw_t),)
    return (
        float(config.panel.Lx_mm),
        float(config.panel.Ly_mm),
        tuple(float(a) for a in config.layup_deg),
        t_sig,
        int(mesh.elements_per_ply),
        float(mesh.in_plane_size_mm),
    )


# Single-entry skeleton cache. A parametric sweep that varies only the
# damage state (sweep_energies) calls build_fe_mesh repeatedly with an
# identical mesh signature; caching the last skeleton makes those calls
# rebuild only the cheap per-element damage factors via apply_damage.
_SKELETON_CACHE: dict = {}


def build_fe_mesh(config: AnalysisConfig, damage: DamageState) -> FeMesh:
    """Build a structured brick mesh for the panel with per-element metadata.

    The mesh spans ``0 <= x <= Lx_mm``, ``0 <= y <= Ly_mm`` and the laminate
    thickness from ``z = 0`` (ply 0) upwards, at the resolution of
    ``config.mesh`` (``MeshParams()`` when ``None``); ``damage`` sets the
    per-element damage factors (see the module docstring). The
    damage-independent skeleton of the last mesh is cached, so a sweep
    that changes only the damage rebuilds only the factors.

    Examples
    --------
    >>> from bvidfe import AnalysisConfig, DamageState, ImpactEvent, PanelGeometry
    >>> cfg = AnalysisConfig(material="IM7/8552", layup_deg=[0, 90, 90, 0],
    ...                      ply_thickness_mm=0.152, panel=PanelGeometry(150, 100),
    ...                      impact=ImpactEvent(energy_J=10.0))
    >>> mesh = build_fe_mesh(cfg, DamageState())
    >>> mesh.n_elements  # 30 x 20 elements in-plane, one per ply
    2400
    """
    sig = _skeleton_signature(config)
    skeleton = _SKELETON_CACHE.get(sig)
    if skeleton is None:
        skeleton = build_fe_mesh_skeleton(config)
        # Keep only the most recent skeleton (sweeps fix the mesh signature
        # for the whole run; this bounds memory to one mesh).
        _SKELETON_CACHE.clear()
        _SKELETON_CACHE[sig] = skeleton
    return apply_damage(skeleton, damage)
