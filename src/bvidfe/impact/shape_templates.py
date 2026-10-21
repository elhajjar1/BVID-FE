"""Layup-dependent 'peanut' template distributing DPA into per-interface delamination ellipses.

Each interface's ellipse is parametrised by:
- aspect ratio AR = 1 + 0.025 * |delta_theta|, clipped to [1, 4]
- orientation = bisector of neighbour-ply angles
- relative size grows from impact face (small) toward back face (large)

Total DPA is enforced by a Brent root-find on a single scalar multiplier so the
polygon-union footprint equals the target DPA within 1%.

With ``max_half_extent_mm`` each ellipse is shrunk (aspect ratio kept) until
its bounding box fits inside the given half-extents about the centroid: the
delaminations stay inside the window the plate was supported on when it was
impacted. The union of differently oriented ellipses spreads well beyond the
circle of the same area (on the 40-ply hard Lovejoy & Scotti coupon, 108 mm
across a 102 mm coupon, which left fe3d a damage band through the whole
section). If the window cannot hold the target area, the saturated footprint
is returned and a ``DPACapClipWarning``-style ``UserWarning`` is raised.

Performance notes (issue #87):
- Per-interface template coefficients (``ar``, ``rel``, ``orient``) are
  precomputed once as NumPy arrays at ``distribute_damage`` entry.
- ``_build`` is vectorized over all interfaces for a given scalar.
- The expensive ``shapely.union_all`` call is invoked only at ~10 anchor
  scalars across the Brent bracket. Brent then operates on a cheap
  ``np.interp`` log-log interpolant; a single true union evaluation
  validates / refines the final scalar.
"""

from __future__ import annotations

import warnings
from typing import List, Optional, Tuple

import numpy as np
from scipy.optimize import brentq

from bvidfe.damage.state import DamageState, DelaminationEllipse


def _aspect_ratio(delta_theta_deg: float) -> float:
    return min(4.0, max(1.0, 1.0 + 0.025 * abs(delta_theta_deg)))


def _orientation_deg(lower_ply_deg: float, upper_ply_deg: float) -> float:
    return 0.5 * (lower_ply_deg + upper_ply_deg)


def _relative_size(interface_index: int, n_interfaces: int) -> float:
    """Weight grows from 0.3 near impact face to 1.0 near back face."""
    z = (interface_index + 1) / n_interfaces  # 0 < z <= 1
    return 0.3 + 0.7 * z


def _template_arrays(
    layup_deg: List[float], n_interfaces: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute per-interface (ar, rel, orient) as NumPy arrays."""
    layup = np.asarray(layup_deg, dtype=float)
    lower = layup[:-1]
    upper = layup[1:]
    dtheta = np.abs(upper - lower)
    ar = np.clip(1.0 + 0.025 * dtheta, 1.0, 4.0)
    orient = 0.5 * (lower + upper)
    idx = np.arange(n_interfaces, dtype=float)
    z = (idx + 1.0) / n_interfaces
    rel = 0.3 + 0.7 * z
    return ar, rel, orient


def distribute_damage(
    layup_deg: List[float],
    target_dpa_mm2: float,
    dent_depth_mm: float,
    fiber_break_radius_mm: float,
    centroid_mm: Tuple[float, float] = (0.0, 0.0),
    max_half_extent_mm: Optional[Tuple[float, float]] = None,
) -> List[DelaminationEllipse]:
    n_plies = len(layup_deg)
    n_interfaces = n_plies - 1
    if n_interfaces <= 0 or target_dpa_mm2 <= 0:
        return []

    ar_arr, rel_arr, orient_arr = _template_arrays(layup_deg, n_interfaces)
    # Per-interface "major / scalar" and "minor / scalar" factors.
    major_factor = ar_arr * rel_arr
    minor_factor = rel_arr
    # Native-Python orientation list (the dataclass stores it as a float).
    orient_list = orient_arr.tolist()
    theta = np.radians(orient_arr)
    cos_t, sin_t = np.abs(np.cos(theta)), np.abs(np.sin(theta))

    def _build(scalar: float) -> List[DelaminationEllipse]:
        # Vectorized ellipse-parameter computation. All semi-axes for a
        # given scalar come out in one NumPy multiplication; we then
        # zip into per-interface dataclasses.
        majors = scalar * major_factor
        minors = scalar * minor_factor
        if max_half_extent_mm is not None:
            # Half-widths of each rotated ellipse's bounding box.
            ex = np.hypot(majors * cos_t, minors * sin_t)
            ey = np.hypot(majors * sin_t, minors * cos_t)
            hx, hy = max_half_extent_mm
            k = np.minimum(1.0, np.minimum(hx / ex, hy / ey))
            majors, minors = majors * k, minors * k
        majors, minors = majors.tolist(), minors.tolist()
        return [
            DelaminationEllipse(
                interface_index=i,
                centroid_mm=centroid_mm,
                major_mm=majors[i],
                minor_mm=minors[i],
                orientation_deg=orient_list[i],
            )
            for i in range(n_interfaces)
        ]

    def _union_area(scalar: float) -> float:
        return DamageState(
            _build(scalar), dent_depth_mm, fiber_break_radius_mm
        ).projected_damage_area_mm2

    lo, hi = 0.1, 50.0
    area_hi = _union_area(hi)
    while area_hi < target_dpa_mm2:
        hi *= 2
        area_prev, area_hi = area_hi, _union_area(hi)
        if max_half_extent_mm is not None and area_hi <= area_prev * (1.0 + 1e-9):
            # Every ellipse is pinned to the window: the area cannot grow.
            warnings.warn(
                f"Target DPA ({target_dpa_mm2:.0f} mm^2) does not fit inside the "
                f"impact window (half-extents {max_half_extent_mm[0]:.1f} x "
                f"{max_half_extent_mm[1]:.1f} mm); using the largest footprint that "
                f"does ({area_hi:.0f} mm^2).",
                UserWarning,
                stacklevel=2,
            )
            return _build(hi)
        if hi > 1e6:
            raise RuntimeError("shape_templates: cannot bracket target DPA")

    # Sample the true union at ~10 log-spaced anchor scalars across the
    # bracket and build a cheap interpolant. Union area scales like
    # scalar**2 to leading order (with overlap corrections), so log-log
    # interpolation is near-linear and stable.
    n_anchors = 10
    anchors = np.geomspace(lo, hi, n_anchors)
    anchor_areas = np.array([_union_area(float(s)) for s in anchors])
    log_anchors = np.log(anchors)
    log_areas = np.log(anchor_areas)
    log_target = np.log(target_dpa_mm2)

    def _interp_residual(scalar: float) -> float:
        return float(np.interp(np.log(scalar), log_anchors, log_areas) - log_target)

    # The log-log relation is monotone, so the residual changes sign
    # within the bracket. Brent on the cheap interpolant gives an
    # initial scalar; we then verify with one true union call. If the
    # interpolant drifts more than 0.5% (well inside the existing 1%
    # contract), we re-solve on the true union from a tightened
    # bracket around the candidate.
    scalar_cheap = brentq(_interp_residual, lo, hi, xtol=1e-3)
    true_area = _union_area(scalar_cheap)
    rel_err = abs(true_area - target_dpa_mm2) / target_dpa_mm2
    if rel_err > 5e-3:
        # Tighten bracket around the candidate and solve on the true union.
        # The interpolant is monotone in scalar, so a 2x neighbourhood
        # always brackets the true root.
        tight_lo = max(lo, scalar_cheap * 0.5)
        tight_hi = min(hi, scalar_cheap * 2.0)
        # Ensure sign change; widen if necessary.
        while _union_area(tight_lo) - target_dpa_mm2 > 0 and tight_lo > lo:
            tight_lo = max(lo, tight_lo * 0.5)
        while _union_area(tight_hi) - target_dpa_mm2 < 0 and tight_hi < hi:
            tight_hi = min(hi, tight_hi * 2.0)
        scalar = brentq(
            lambda s: _union_area(s) - target_dpa_mm2,
            tight_lo,
            tight_hi,
            xtol=1e-3,
        )
    else:
        scalar = scalar_cheap

    return _build(scalar)
