"""3D finite-element tier for BVID-FE.

The CAI residual is the minimum of two channels:
- first-ply failure on the damaged 3D mesh (this module's distinctive
  value: 3D stress states with per-element damage factors), LaRC05;
- the far-field stress at which the weakest buckled sublaminate grows its
  delamination (delegated to :mod:`bvidfe.analysis.semi_analytical`).

The TAI residual is first-ply failure alone, Tsai-Wu.

First-ply failure: one linear static solve of the structured hex mesh
(:mod:`bvidfe.analysis.fe_mesh`) at a reference end strain along x, scaled
analytically to the load at which the failure index reaches 1. The damage
zone (every ply under the union of the delamination footprints) is softened
in-plane by ``fe_mesh.DAMAGE_ZONE_INPLANE_FACTOR``, each element weighted by
the fraction of its area inside the zone, and its elements see the effective
stress of strain equivalence. The failure index, in each ply's material
frame, is averaged over each element's Gauss points and then over the
elements of the same through-thickness layer within
``FPF_AVERAGING_RADIUS_MM`` (2.5 mm) of it before it is checked against 1;
at the default 5 mm mesh that disk holds only the element itself. The
failure stress is the reaction force on the loaded edge over the gross
section. ``BvidAnalysis`` divides the damaged result by the same quantity
on the undamaged mesh and applies the ratio to the pristine strength.

Buckling onset (whole panel and sublaminate, Rayleigh-Ritz closed form) is
reported by :func:`fe3d_cai_buckling` but does not enter the residual.

Problems above ``FE3D_MAX_DOF`` degrees of freedom (150,000, or the
``BVIDFE_FE3D_MAX_DOF`` environment variable) raise :class:`FE3DSizeError`.
Stage timings go to stderr on the ``bvidfe.fe3d`` logger, at the level set
by ``BVIDFE_LOG_LEVEL`` (default ``INFO``).
"""

from __future__ import annotations

import logging
import os
import sys
import time
from typing import List

import numpy as np
from scipy.ndimage import convolve

from bvidfe._types import _CRITERION_NAMES, CriterionName
from bvidfe.analysis.config import AnalysisConfig
from bvidfe.analysis.fe_mesh import (
    FeMesh,
    build_fe_mesh,
    element_grid_shape,
    estimate_fe_mesh_size,
)
from bvidfe.core.laminate import Laminate
from bvidfe.damage.state import DamageState
from bvidfe.elements.hex8 import Hex8Element, build_geometry_table
from bvidfe.failure.larc05 import larc05_index, larc05_index_batch
from bvidfe.failure.puck import puck_index, puck_index_batch
from bvidfe.failure.tsai_wu import tsai_wu_index, tsai_wu_index_batch
from bvidfe.analysis.semi_analytical import (
    find_critical_interface,
    panel_buckling_load,
    sublaminate_buckling_load,
    weakest_sublaminate_growth,
)
from bvidfe.solver.boundary import uniaxial_x_bcs
from bvidfe.solver.static import solve_linear_static

# Module-level logger, one line per FE stage. As a library, bvidfe only adds
# a NullHandler; the CLI and the Streamlit app call enable_stderr_logging()
# so long runs stay observable from a terminal. Python API users opt in with
# that call or their own logging configuration.
_log = logging.getLogger("bvidfe.fe3d")
_log.addHandler(logging.NullHandler())

# Radius (mm) of the in-plane disk over which the failure index is averaged
# before first-ply failure is checked (average-stress criterion, Whitney &
# Nuismer 1974). A pointwise check has no length scale: at the re-entrant
# corners of the damage-zone footprint the stress is singular, so the
# pointwise failure load kept moving as the mesh was refined (cross-ply
# knockdown 0.175 / 0.223 / 0.211 at 5 / 2.5 / 2 mm, old damage template).
# Averaged over 2.5 mm, with the zone weighted by element area fraction and
# the disk by element overlap (_disk_kernel), the 8-ply cross-ply 15 J CAI
# knockdown is 0.217 / 0.185 / 0.170 / 0.174 at 5 / 2.5 / 2 / 1.67 mm: the
# fine meshes converge, but 5 mm under-resolves a zone this small (about
# 710 mm^2, a few elements across). 2.5 mm is about the characteristic
# distance of CFRP notched strength (1-4 mm); at a 5 mm mesh the disk holds
# only the element itself.
FPF_AVERAGING_RADIUS_MM = 2.5


def _resolve_log_level(default: str = "INFO") -> int:
    """Resolve BVIDFE_LOG_LEVEL → numeric level, falling back to ``default``
    with a stderr warning when the value is not a recognised level name.

    Without this guard a typo (``BVIDFE_LOG_LEVEL=DUBG``) raised ValueError
    at module import time and prevented the package from being imported.
    """
    raw = os.environ.get("BVIDFE_LOG_LEVEL", default).upper()
    level = logging.getLevelName(raw)
    if isinstance(level, int):
        return level
    sys.stderr.write(
        f"[bvidfe] BVIDFE_LOG_LEVEL={raw!r} is not a valid log level "
        f"(use DEBUG/INFO/WARNING/ERROR/CRITICAL); defaulting to {default}.\n"
    )
    return logging.getLevelName(default)


def _resolve_max_dof(default: int = 150000) -> int:
    """Resolve BVIDFE_FE3D_MAX_DOF → positive int, falling back to ``default``
    with a stderr warning when the value is not a positive integer."""
    raw = os.environ.get("BVIDFE_FE3D_MAX_DOF")
    if raw is None:
        return default
    try:
        v = int(raw)
        if v <= 0:
            raise ValueError(f"must be > 0 (got {v})")
        return v
    except ValueError as exc:
        sys.stderr.write(
            f"[bvidfe] BVIDFE_FE3D_MAX_DOF={raw!r} is not a positive int "
            f"({exc}); defaulting to {default}.\n"
        )
        return default


class _StderrHandler(logging.StreamHandler):
    """StreamHandler that writes to whatever ``sys.stderr`` is at emit time
    (so a replaced or captured stderr never leaves it holding a closed
    stream)."""

    def __init__(self) -> None:
        super().__init__(sys.stderr)

    @property  # type: ignore[override]
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, value) -> None:
        pass


def enable_stderr_logging() -> None:
    """Send ``bvidfe.fe3d`` stage logs to stderr at ``BVIDFE_LOG_LEVEL``
    (default INFO). Idempotent: repeated calls add one handler only."""
    if not any(isinstance(h, _StderrHandler) for h in _log.handlers):
        h = _StderrHandler()
        h.setFormatter(logging.Formatter("[bvidfe.fe3d %(asctime)s] %(message)s", "%H:%M:%S"))
        _log.addHandler(h)
    _log.setLevel(_resolve_log_level())


# An explicit BVIDFE_LOG_LEVEL still sets the logger's level for library use.
if "BVIDFE_LOG_LEVEL" in os.environ:
    _log.setLevel(_resolve_log_level())


def _t(msg: str, t0: float) -> None:
    """Log a timing line: message + seconds since t0."""
    _log.info("%s (%.2fs)", msg, time.time() - t0)


# Hard cap on fe3d problem size, so an oversized mesh raises FE3DSizeError
# instead of exhausting memory inside SuperLU (native code: an OOM there
# kills the process rather than raising). Peak RSS of a whole fe3d CAI run
# (two solves) on an 8-ply 150 x 100 mm panel, OMP_NUM_THREADS=1:
#   mesh    DOFs   ILU-CG (default)  direct LU fallback  old penalty + spsolve
#   5 mm    17.6k  0.20 GB           0.33 GB             0.55 GB
#   2.5 mm  67.5k  0.46 GB           1.3 GB              3.0 GB
#   1.25mm  265k   -                 -                   killed (OOM, 15 GB)
# ILU-CG memory grows about linearly (~1 GB at 150k DOFs); the direct LU,
# used only when ILU-CG and ILU-GMRES both fail, grows as ~n^1.4 in fill
# (113M factor entries at 63k free DOFs on a cross-ply panel), so ~4-8 GB
# at 150k. The cap keeps that worst case within a workstation's memory; it
# admits 8-ply 150 x 100 mm panels down to ~1.7 mm and 16-ply ones at
# 2.5 mm. Override via the BVIDFE_FE3D_MAX_DOF env var at your own risk; an
# invalid value falls back to the default with a stderr warning instead of
# raising at import time.
FE3D_MAX_DOF: int = _resolve_max_dof()


class FE3DSizeError(RuntimeError):
    """Raised when an fe3d problem exceeds the safe-size cap."""


def _guard_problem_size(cfg: AnalysisConfig) -> None:
    """Raise FE3DSizeError if the mesh would exceed the safe-size cap.

    Called at the top of fe3d_cai_buckling / fe3d_cai / fe3d_tai so callers
    get a clean Python exception instead of a native-code crash.
    """
    stats = estimate_fe_mesh_size(cfg)
    if stats["n_dof"] > FE3D_MAX_DOF:
        raise FE3DSizeError(
            f"fe3d problem too large: {stats['n_elements']:,} elements / "
            f"{stats['n_dof']:,} DOFs exceeds the safe-size cap of "
            f"{FE3D_MAX_DOF:,} DOFs. Increase MeshParams.in_plane_size_mm, "
            f"decrease elements_per_ply, or switch to tier='empirical' / "
            f"tier='semi_analytical'. Override the cap via the "
            f"BVIDFE_FE3D_MAX_DOF environment variable at your own risk."
        )


# Voigt indices: [xx, yy, zz, yz, xz, xy] = [0, 1, 2, 3, 4, 5]
# In-plane components live on rows/cols {0, 1, 5}; out-of-plane on {2, 3, 4}.
# After rotation about z (Hex8Element.ply_angle_deg) these blocks remain
# decoupled in C_global, so the mask below applies cleanly in the global
# frame. See fe_mesh.DAMAGE_OOP_FACTOR / DAMAGE_FIBER_BREAK_INPLANE_FACTOR
# docstrings for the physical model.
_INPLANE_VOIGT = np.array([0, 1, 5])


def _build_elements(mesh: FeMesh, lam: Laminate) -> List[Hex8Element]:
    """Build a Hex8Element for each mesh element with component-wise damage scaling.

    Two factors are applied to the 6x6 elasticity matrix:
      - `mesh.damage_factors[e]` scales OOP and in-plane/OOP cross terms
      - `mesh.in_plane_damage_factors[e]` scales the in-plane sub-block
    """
    elements: List[Hex8Element] = []
    # On the regular cuboid grid produced by build_fe_mesh, every element of a
    # given ply is a pure translate of every other, so its 2x2x2 Gauss-point
    # (B, detJ, gradN) table is identical. Compute it once per ply (keyed by
    # ply index AND element geometry — the latter keeps mixed-thickness plies,
    # issue #5/PR#13, correct) and share it across that ply's elements.
    geom_cache: dict[tuple, object] = {}
    for eidx in range(mesh.n_elements):
        node_ids = mesh.element_connectivity[eidx]
        node_coords = mesh.node_coords[node_ids]
        mat = lam.material
        ply_angle = float(mesh.ply_angles_deg[eidx])
        ply_i = int(mesh.ply_indices[eidx])
        dims = node_coords.max(axis=0) - node_coords.min(axis=0)
        cache_key = (
            ply_i,
            round(float(dims[0]), 9),
            round(float(dims[1]), 9),
            round(float(dims[2]), 9),
        )
        b_table = geom_cache.get(cache_key)
        if b_table is None:
            b_table = build_geometry_table(node_coords)
            geom_cache[cache_key] = b_table
        elem = Hex8Element(node_coords, mat, ply_angle_deg=ply_angle, b_table=b_table)
        f_oop = float(mesh.damage_factors[eidx])
        f_ip = float(mesh.in_plane_damage_factors[eidx])
        # Build a 6x6 element-wise scaling mask: start from f_oop everywhere
        # (OOP block + Poisson cross-coupling), then overwrite the 3x3
        # in-plane sub-block on rows/cols {0, 1, 5} with f_ip.
        mask = np.full((6, 6), f_oop)
        mask[np.ix_(_INPLANE_VOIGT, _INPLANE_VOIGT)] = f_ip
        elem._C_global = elem._C_global * mask
        elements.append(elem)
    return elements


def _resolve_material(cfg: AnalysisConfig):
    if isinstance(cfg.material, str):
        from bvidfe.core.material import MATERIAL_LIBRARY

        return MATERIAL_LIBRARY[cfg.material]
    return cfg.material


def _fe3d_preflight(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    *,
    label: str = "fe3d",
) -> tuple[FeMesh, List[Hex8Element], float]:
    """Run the size guard, build the damaged mesh, and build elements.

    Factors the duplicated pre-flight prologue used by every fe3d entry
    point (``fe3d_cai_buckling``, ``_fe3d_cai_first_ply_failure``,
    ``fe3d_tai``). The ``label`` differentiates the per-stage log line so
    a stderr trace still distinguishes "fe3d FPF" from "fe3d buckling".

    Returns
    -------
    (mesh, elements, t0)
        ``t0`` is ``time.time()`` captured before the mesh build, suitable
        for relative timing logs in the caller.
    """
    _guard_problem_size(cfg)
    t0 = time.time()
    _log.info("%s: mesh build start", label)
    mesh = build_fe_mesh(cfg, damage)
    _t(f"mesh build done: {mesh.n_elements} elements, {mesh.n_dof} DOFs", t0)
    elements = _build_elements(mesh, lam)
    _t("elements built", t0)
    return mesh, elements, t0


def _require_known_criterion(criterion: str) -> None:
    """Raise ValueError for a criterion name the FPF solve cannot dispatch.

    Without this guard any unrecognised string (a typo, or a criterion added
    to ``CriterionName`` without an analytic-scaling branch here) fell through
    to the Tsai-Wu ``else`` branch and silently returned Tsai-Wu results.
    """
    if criterion not in _CRITERION_NAMES:
        raise ValueError(
            f"unknown failure criterion {criterion!r}; "
            f"valid options are {sorted(_CRITERION_NAMES)}"
        )


def _solve_failure_strain_analytic(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    strain_sign: int,
    criterion: CriterionName,
    strain_cap: float = 0.05,
    averaging_radius_mm: float | None = FPF_AVERAGING_RADIUS_MM,
) -> float:
    """Find the applied strain magnitude at which max failure index hits 1,
    using a single FE solve + analytic scaling.

    Linear elasticity: stress field scales linearly with applied strain. So if
    we solve once at a reference strain `eps_ref`, the stress at any strain
    multiplier `c` is just `c * sigma_ref`. Tsai-Wu is quadratic in stress:
        idx(c) = c * (F . sigma_ref) + c^2 * (sigma_ref . F_ij . sigma_ref)
    At the critical strain, idx = 1; positive root of the quadratic.

    LaRC05 is piecewise-quadratic (tensile vs. compressive modes branch on
    the sign of s1 / s2) but each branch is still quadratic in c for a fixed
    stress direction, so we evaluate all four mode branches at a reference
    strain and pick the binding one analytically.

    Puck's exposure factor is positively homogeneous of degree 1 in stress
    (the fiber term ``|s1|/X`` and every action-plane term
    ``sqrt(quadratic) + linear * sigma_n`` scale linearly for c > 0, and
    scaling by c > 0 never flips a sign-selected branch), so
    idx(c) = c * idx(1) and c_crit = 1 / idx_ref exactly.

    Stresses are rotated into each ply's material frame
    (``Hex8Element.stress_at_gauss_points_material``) before the criterion
    sees them, as effective stresses in the damage zone
    (``_failure_stress_material``); both maps are linear, so the scaling
    arguments above hold.

    Raises ``ValueError`` for a criterion outside ``CriterionName``.

    This replaces the prior 10-12 iteration bisection (each iteration
    reassembled and re-solved the full FE system), giving a ~10x speedup
    for the FPF path with identical results on the quadratic branches.

    The per-element inner Gauss-point loop is vectorised via the batch
    criterion helpers (``larc05_index_batch`` / ``tsai_wu_index_batch``).
    Numerical equivalence with the prior scalar form is locked by
    ``tests/analysis/test_fpf_strain_solve_equivalence.py`` against the
    ``_solve_failure_strain_analytic_scalar_ref`` reference implementation
    kept in this module purely for the equivalence test, which passes
    ``averaging_radius_mm=None`` (pointwise, as the reference is); see
    ``_critical_strain_multiplier``.
    """
    _require_known_criterion(criterion)
    u_ref = _reference_displacement(cfg, mesh, elements, strain_sign, strain_cap)
    c_crit_min = _critical_strain_multiplier(
        cfg, mesh, elements, u_ref, criterion, averaging_radius_mm
    )
    if not np.isfinite(c_crit_min) or c_crit_min <= 0:
        return strain_cap  # nothing failed up to strain_cap
    return min(strain_cap, c_crit_min * strain_cap)


def _reference_displacement(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    strain_sign: int,
    strain_cap: float,
) -> np.ndarray:
    """Displacement field at the reference strain ``strain_sign * strain_cap``."""
    bcs = uniaxial_x_bcs(mesh.node_coords, strain_sign * strain_cap, boundary=cfg.panel.boundary)
    return solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs)


def _failure_stress_material(
    elem: Hex8Element, u_elem: np.ndarray, in_plane_factor: float
) -> np.ndarray:
    """Ply-frame stress (n_gp, 6) the failure criterion sees.

    For an element in the damage zone (``in_plane_factor < 1``) this is the
    effective stress of strain equivalence (Lemaitre): damaged material
    fails at the strain intact material would, so its strength drops with
    its stiffness. The in-plane rows come from the pristine stiffness at the
    element's actual strain; the out-of-plane rows keep the actual (damaged)
    tractions. Dividing the in-plane stress by the factor instead would also
    scale the through-thickness Poisson term in those rows.
    """
    if in_plane_factor >= 1.0:
        return elem.stress_at_gauss_points_material(u_elem)
    C_eff = elem._C_global.copy()
    C_eff[_INPLANE_VOIGT, :] = elem._compute_global_stiffness()[_INPLANE_VOIGT, :]
    return elem.stress_at_gauss_points_material(u_elem, C=C_eff)


def _critical_strain_multiplier(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    u_ref: np.ndarray,
    criterion: CriterionName,
    averaging_radius_mm: float | None = FPF_AVERAGING_RADIUS_MM,
) -> float:
    """Smallest multiplier ``c`` on ``u_ref`` at which the failure index
    reaches 1 (``inf`` when nothing is loaded towards failure). See
    ``_solve_failure_strain_analytic`` for the scaling.

    With ``averaging_radius_mm`` (default ``FPF_AVERAGING_RADIUS_MM``) the
    index is averaged over each element's Gauss points and then over a disk
    of that radius in the same through-thickness layer, each element
    weighted by the fraction of its area inside the disk
    (``_disk_kernel``); ``None`` checks every Gauss point on its own."""
    if averaging_radius_mm is None:
        return _critical_strain_multiplier_pointwise(cfg, mesh, elements, u_ref, criterion)
    a, b = _index_coefficients(cfg, mesh, elements, u_ref, criterion)
    nz, ny, nx = element_grid_shape(cfg)
    if a.size != nz * ny * nx:
        raise ValueError(f"mesh has {a.size} elements, expected a {nz}x{ny}x{nx} grid")
    kernel = _disk_kernel(cfg.panel.Lx_mm / nx, cfg.panel.Ly_mm / ny, averaging_radius_mm)
    # Near the panel edge the disk is clipped: average over what is inside.
    weight = convolve(np.ones((ny, nx)), kernel, mode="constant")
    a_avg = np.stack(
        [convolve(layer, kernel, mode="constant") / weight for layer in a.reshape(nz, ny, nx)]
    )
    b_avg = np.stack(
        [convolve(layer, kernel, mode="constant") / weight for layer in b.reshape(nz, ny, nx)]
    )
    return float(_smallest_positive_root(a_avg.ravel(), b_avg.ravel()).min())


def _disk_kernel(hx: float, hy: float, radius_mm: float, samples: int = 16) -> np.ndarray:
    """Weights of the elements of an ``hx`` x ``hy`` grid around one at the
    origin: the fraction of each element's area inside a disk of
    ``radius_mm`` centred on it, from ``samples`` x ``samples`` points per
    element. Counting whole elements whose centres lie inside the disk
    averaged over only the element and its four neighbours at 2-2.5 mm
    meshes, an area that shrank with the mesh, and the cross-ply 15 J
    knockdown kept falling with it (0.186 / 0.165 at 2.5 / 2 mm); weighted
    by overlap the disk keeps its area and the knockdown settles (0.185 /
    0.170 / 0.174 at 2.5 / 2 / 1.67 mm). At 5 mm the disk fits inside the
    element, so the kernel is the element alone."""
    rx = int(np.ceil(radius_mm / hx + 0.5))
    ry = int(np.ceil(radius_mm / hy + 0.5))
    t = (np.arange(samples) + 0.5) / samples - 0.5
    x = (np.arange(-rx, rx + 1)[:, None] + t[None, :]) * hx  # (2rx+1, samples)
    y = (np.arange(-ry, ry + 1)[:, None] + t[None, :]) * hy  # (2ry+1, samples)
    inside = x[None, :, None, :] ** 2 + y[:, None, :, None] ** 2 <= radius_mm**2
    return inside.mean(axis=(2, 3))  # (2ry+1, 2rx+1)


def _index_coefficients(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    u_ref: np.ndarray,
    criterion: CriterionName,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-element ``(a, b)`` with the Gauss-point mean of the failure index
    at multiplier ``c`` on ``u_ref`` equal to ``a c + b c^2``."""
    material = _resolve_material(cfg)
    in_plane = mesh.in_plane_damage_factors
    a = np.zeros(len(elements))
    b = np.zeros(len(elements))
    for eidx, elem in enumerate(elements):
        dof_map = mesh.element_dof_maps[eidx]
        sigma_ref = _failure_stress_material(elem, u_ref[dof_map], in_plane[eidx])
        if criterion == "larc05":
            b[eidx] = larc05_index_batch(material, sigma_ref).mean()
        elif criterion == "puck":
            a[eidx] = puck_index_batch(material, sigma_ref).mean()
        else:
            idx_ref = tsai_wu_index_batch(material, sigma_ref)
            idx_2 = tsai_wu_index_batch(material, 2.0 * sigma_ref)
            b_gp = (idx_2 - 2.0 * idx_ref) / 2.0
            b[eidx] = b_gp.mean()
            a[eidx] = (idx_ref - b_gp).mean()
    return a, b


def _smallest_positive_root(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Smallest ``c > 0`` with ``a c + b c^2 = 1`` elementwise (``inf``
    where there is none)."""
    tiny = 1e-14
    c = np.full(a.shape, np.inf)
    lin = (np.abs(b) < tiny) & (a >= tiny)
    c = np.where(lin, 1.0 / np.where(lin, a, 1.0), c)
    disc = a * a + 4.0 * b
    quad = (np.abs(b) >= tiny) & (disc >= 0)
    root = (-a + np.sqrt(np.where(quad, disc, 0.0))) / np.where(quad, 2.0 * b, 1.0)
    c = np.where(quad, root, c)
    return np.where(c > 0, c, np.inf)


def _critical_strain_multiplier_pointwise(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    u_ref: np.ndarray,
    criterion: CriterionName,
) -> float:
    """``_critical_strain_multiplier`` checked at every Gauss point on its
    own (``inf`` when no Gauss point is loaded towards failure)."""
    material = _resolve_material(cfg)
    in_plane = mesh.in_plane_damage_factors
    c_crit_min = np.inf
    for eidx, elem in enumerate(elements):
        dof_map = mesh.element_dof_maps[eidx]
        sigma_ref = _failure_stress_material(elem, u_ref[dof_map], in_plane[eidx])  # (n_gp, 6)

        if criterion == "larc05":
            idx_ref = larc05_index_batch(material, sigma_ref)  # (n_gp,)
            valid = idx_ref > 0
            if not bool(valid.any()):
                continue
            # LaRC05 modes are sums of squared normalised stresses, so
            # idx(c) = c^2 * idx(1)  ->  c_crit = 1 / sqrt(idx_ref).
            c_crit = np.where(valid, 1.0 / np.sqrt(np.where(valid, idx_ref, 1.0)), np.inf)
            c_crit_elem_min = float(c_crit.min())
        elif criterion == "puck":
            idx_ref = puck_index_batch(material, sigma_ref)  # (n_gp,)
            valid = idx_ref > 0
            if not bool(valid.any()):
                continue
            # Puck is degree-1 homogeneous: idx(c) = c * idx(1)  ->  c_crit = 1 / idx_ref.
            c_crit = np.where(valid, 1.0 / np.where(valid, idx_ref, 1.0), np.inf)
            c_crit_elem_min = float(c_crit.min())
        else:
            # Tsai-Wu: idx(c) = a*c + b*c^2. Solve a, b from two samples
            # (c=1 and c=2). System: a + b = idx_ref;  2a + 4b = idx_2.
            idx_ref = tsai_wu_index_batch(material, sigma_ref)  # (n_gp,)
            sigma_2 = 2.0 * sigma_ref
            idx_2 = tsai_wu_index_batch(material, sigma_2)  # (n_gp,)
            b = (idx_2 - 2.0 * idx_ref) / 2.0
            a = idx_ref - b
            valid = idx_ref > 0
            tiny = 1e-14
            # Initialise c_crit = +inf so masked-out points never drive the min.
            c_crit = np.full_like(idx_ref, np.inf, dtype=float)
            # Branch 1: |b| < tiny and |a| > tiny -> c = 1/a.
            mask_lin = valid & (np.abs(b) < tiny) & (np.abs(a) >= tiny)
            # Use np.where to keep divisor non-zero for inactive lanes.
            c_lin = 1.0 / np.where(mask_lin, a, 1.0)
            c_crit = np.where(mask_lin, c_lin, c_crit)
            # Branch 2: |b| >= tiny and disc >= 0 -> c = (-a + sqrt(disc)) / (2b).
            mask_quad = valid & (np.abs(b) >= tiny)
            disc = a * a + 4.0 * b
            mask_quad_real = mask_quad & (disc >= 0)
            sqrt_disc = np.sqrt(np.where(mask_quad_real, disc, 0.0))
            c_quad = (-a + sqrt_disc) / np.where(mask_quad_real, 2.0 * b, 1.0)
            c_crit = np.where(mask_quad_real, c_quad, c_crit)
            # Filter c_crit > 0 (Tsai-Wu can produce a non-physical
            # negative root when the linear and quadratic terms have
            # opposite signs in a way that doesn't bracket failure).
            c_crit = np.where(c_crit > 0, c_crit, np.inf)
            c_crit_elem_min = float(c_crit.min())

        if c_crit_elem_min < c_crit_min:
            c_crit_min = c_crit_elem_min
    return c_crit_min


def _loaded_edge_force(mesh: FeMesh, elements: List[Hex8Element], u: np.ndarray) -> float:
    """Magnitude of the x reaction force (N) on the loaded x_max face.

    The reaction at a constrained DOF is the internal force ``(K u)_i``,
    summed here from the element internal forces ``K_e u_e`` of the elements
    touching the face, so no global matrix is assembled.
    """
    x = mesh.node_coords[:, 0]
    on_face = np.abs(x - x.max()) < 1e-9
    force = 0.0
    for eidx in np.flatnonzero(on_face[mesh.element_connectivity].any(axis=1)):
        dof_map = mesh.element_dof_maps[eidx]
        f_elem = elements[eidx].stiffness_matrix() @ u[dof_map]
        local = np.flatnonzero(on_face[mesh.element_connectivity[eidx]])
        force += float(f_elem[3 * local].sum())
    return abs(force)


def _first_ply_failure_stress(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    lam: Laminate,
    strain_sign: int,
    criterion: CriterionName,
    strain_cap: float = 0.05,
) -> float:
    """Far-field stress (MPa) at first-ply failure: the reaction force on
    the loaded edge at the failure strain over the gross section
    ``Ly * h``.

    Multiplying the failure strain by the pristine CLT modulus credited a
    softened (damaged) panel with load it cannot carry: a panel softened
    uniformly x0.3 fails at 3.3x the strain and so reported 3.3x its
    intact strength. The failure multiplier and the reaction force come
    from the same FE solve.
    """
    _require_known_criterion(criterion)
    u_ref = _reference_displacement(cfg, mesh, elements, strain_sign, strain_cap)
    c_crit = _critical_strain_multiplier(cfg, mesh, elements, u_ref, criterion)
    if not np.isfinite(c_crit) or c_crit <= 0:
        c_crit = 1.0  # nothing failed up to strain_cap
    c_crit = min(c_crit, 1.0)
    area = cfg.panel.Ly_mm * lam.thickness_mm
    return c_crit * _loaded_edge_force(mesh, elements, u_ref) / area


def _solve_failure_strain_analytic_scalar_ref(
    cfg: AnalysisConfig,
    mesh: FeMesh,
    elements: List[Hex8Element],
    strain_sign: int,
    criterion: CriterionName,
    strain_cap: float = 0.05,
) -> float:
    """Reference (pre-vectorisation) scalar implementation of
    ``_solve_failure_strain_analytic``.

    Kept intact for the equivalence test in
    ``tests/analysis/test_fpf_strain_solve_equivalence.py`` so any future
    edit to the production routine immediately surfaces a numerical
    drift. NOT called from production code; do not use directly.
    """
    _require_known_criterion(criterion)
    material = _resolve_material(cfg)

    bcs = uniaxial_x_bcs(mesh.node_coords, strain_sign * strain_cap, boundary=cfg.panel.boundary)
    u_ref = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs)

    c_crit_min = np.inf
    for eidx, elem in enumerate(elements):
        dof_map = mesh.element_dof_maps[eidx]
        sigma_field_ref = _failure_stress_material(
            elem, u_ref[dof_map], mesh.in_plane_damage_factors[eidx]
        )
        for gp in range(sigma_field_ref.shape[0]):
            sigma_ref = sigma_field_ref[gp]
            if criterion == "larc05":
                idx_ref = larc05_index(material, sigma_ref)
            elif criterion == "puck":
                idx_ref = puck_index(material, sigma_ref)
            else:
                idx_ref = tsai_wu_index(material, sigma_ref)
            if idx_ref <= 0:
                continue
            if criterion == "larc05":
                c_crit = 1.0 / np.sqrt(idx_ref)
            elif criterion == "puck":
                c_crit = 1.0 / idx_ref
            else:
                sigma_2 = 2.0 * sigma_ref
                idx_2 = tsai_wu_index(material, sigma_2)
                b = (idx_2 - 2.0 * idx_ref) / 2.0
                a = idx_ref - b
                if abs(b) < 1e-14:
                    if abs(a) < 1e-14:
                        continue
                    c_crit = 1.0 / a
                else:
                    disc = a * a + 4.0 * b
                    if disc < 0:
                        continue
                    c_crit = (-a + np.sqrt(disc)) / (2.0 * b)
                if c_crit <= 0:
                    continue
            if c_crit < c_crit_min:
                c_crit_min = c_crit

    if not np.isfinite(c_crit_min) or c_crit_min <= 0:
        return strain_cap
    return min(strain_cap, c_crit_min * strain_cap)


def _fe3d_cai_first_ply_failure(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
    criterion: CriterionName = "larc05",
) -> float:
    """3D FE compression-after-impact residual strength via first-ply-failure (MPa).

    Scales one elastic solve of the damaged mesh to the load at which the
    selected failure criterion's index first reaches 1, and returns the
    loaded-edge reaction force at that load over the gross area (see
    ``_first_ply_failure_stress``). Default ``criterion="larc05"``; pass
    ``"tsai_wu"`` or ``"puck"`` to use those criteria instead. Unknown names
    raise ``ValueError``.
    """
    mesh, elements, t0 = _fe3d_preflight(cfg, damage, lam, label="fe3d FPF")
    sigma = _first_ply_failure_stress(cfg, mesh, elements, lam, strain_sign=-1, criterion=criterion)
    _t("FPF analytic solve done", t0)
    _log.info(
        "fe3d FPF done: residual = %.1f MPa (total %.2fs)",
        min(sigma, sigma_pristine_MPa),
        time.time() - t0,
    )
    return min(sigma, sigma_pristine_MPa)


def fe3d_cai_buckling(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
    sigma_ref_MPa: float = 1.0,
) -> tuple[float, float, List[str]]:
    """3D FE compression-after-impact buckling channel — closed-form delegation.

    Issue #129: the previous 3D K_g eigensolve under-predicted by ~100x on
    realistic panels because the constant-prestress + 3-2-1 rigid-body BC
    formulation didn't actually enforce a compressive state, and plain Hex8
    locks too aggressively for thin laminates at practical mesh sizes.
    Investigation also showed that even a proper static-prestress
    reformulation did not converge to closed-form within an order of
    magnitude. The Rayleigh-Ritz closed form (Timoshenko & Gere §9.2;
    Reddy 4.4.4) is exact for SSSS orthotropic rectangles, so the buckling
    channel now delegates to it via
    :func:`bvidfe.analysis.semi_analytical.panel_buckling_load`.

    For damaged panels the sublaminate-over-delamination check from
    :func:`bvidfe.analysis.semi_analytical.sublaminate_buckling_load` also
    runs; the minimum of the full-panel and worst-sublaminate buckling
    stress is returned.

    This is buckling *onset* and does not set the fe3d residual strength
    (see :func:`fe3d_cai`): a buckled sublaminate keeps carrying load until
    its delamination grows, and whole-panel buckling is suppressed by the
    anti-buckling guides of a CAI fixture. ``BvidAnalysis`` reports the
    value in ``AnalysisResults.buckling_eigenvalues``.
    """
    boundary = cfg.panel.boundary
    notes: List[str] = []

    t_total = float(sum(lam.ply_thicknesses_mm))
    N_cr_panel = panel_buckling_load(lam, cfg.panel.Lx_mm, cfg.panel.Ly_mm, boundary)
    sigma_critical = N_cr_panel / t_total if t_total > 0 else float("inf")

    if damage.delaminations:
        crit_idx = find_critical_interface(damage, lam)
        if crit_idx is not None:
            ellipses = [e for e in damage.delaminations if e.interface_index == crit_idx]
            critical_ellipse = max(ellipses, key=lambda e: e.area_mm2)
            N_cr_sub = sublaminate_buckling_load(lam, critical_ellipse, boundary=boundary)
            thicknesses = lam.ply_thicknesses_mm
            upper_t = sum(thicknesses[: crit_idx + 1])
            lower_t = sum(thicknesses[crit_idx + 1 :])
            sub_t = min(upper_t, lower_t)
            if sub_t > 0:
                sigma_sublam = N_cr_sub / sub_t
                sigma_critical = min(sigma_critical, sigma_sublam)

    if not np.isfinite(sigma_critical) or sigma_critical <= 0:
        note = (
            "fe3d buckling: closed-form Rayleigh-Ritz returned a degenerate "
            "result, so no buckling onset is reported"
        )
        return sigma_pristine_MPa, 0.0, [note]

    lambda_crit = sigma_critical / sigma_ref_MPa
    _log.info(
        "fe3d buckling: sigma_crit=%.1f MPa (closed-form Rayleigh-Ritz)",
        sigma_critical,
    )
    return min(sigma_critical, sigma_pristine_MPa), lambda_crit, notes


def fe3d_cai(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
) -> float:
    """3D FE compression-after-impact residual strength (MPa).

    The smaller of the far-field stress at which the weakest delamination
    grows by sublaminate buckling
    (:func:`bvidfe.analysis.semi_analytical.weakest_sublaminate_growth`)
    and first-ply failure on the damaged mesh
    (``_fe3d_cai_first_ply_failure``, capped at ``sigma_pristine_MPa``).

    Whole-panel buckling (:func:`fe3d_cai_buckling`) is left out: a CAI
    fixture's anti-buckling guides suppress it, so it is a property of the
    test setup rather than of the damage.
    """
    growth = weakest_sublaminate_growth(lam, damage, boundary=cfg.panel.boundary)
    sigma_growth = growth.growth_stress_MPa if growth is not None else float("inf")
    sigma_fpf = _fe3d_cai_first_ply_failure(cfg, damage, lam, sigma_pristine_MPa)
    return min(sigma_growth, sigma_fpf)


def fe3d_tai(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    sigma_pristine_MPa: float,
) -> float:
    """3D FE tension-after-impact residual strength (MPa).

    First-ply failure (Tsai-Wu) on the damaged mesh in tension, as the
    loaded-edge reaction force over the gross section, capped at
    ``sigma_pristine_MPa``. Raises :class:`FE3DSizeError` above the size cap.
    """
    mesh, elements, _t0 = _fe3d_preflight(cfg, damage, lam, label="fe3d TAI")
    sigma = _first_ply_failure_stress(cfg, mesh, elements, lam, strain_sign=+1, criterion="tsai_wu")
    return min(sigma, sigma_pristine_MPa)
