"""Layup-dependent 'peanut' template distributing DPA into per-interface delamination ellipses.

Ply ``k`` of ``layup_deg`` sits at ``z_k`` with ply 0 at ``z = -h/2``, and ply
0 is taken as the impacted face (as the size ramp below and the C-scan
schema's ``interface_index`` assume). Interface ``i`` lies between ply ``i``
(impact-face side) and ply ``i + 1`` (back-face side).

Fibre directions are axial (``theta`` and ``theta + 180`` are the same ply), so
every angle below is taken modulo 180 deg and ``-45`` / ``135`` or ``90`` /
``-90`` give the same footprint. Each interface's ellipse is parametrised by:

- ply mismatch ``dtheta`` = the acute angle between the two fibre directions,
  in [0, 90] deg;
- aspect ratio AR = 1 + 0.025 * dtheta (1 to 3.25; the [1, 4] clip never binds);
- orientation = fibre angle of the back-face-side ply (ply ``i + 1``),
  reported in (-90, 90]: the peanut-shaped delamination at an interface
  elongates along the fibres of the ply below it, on the side away from the
  impactor (Liu 1988, *J. Compos. Mater.* 22, 674-692; Abrate 1998, *Impact
  on Composite Structures*);
- relative size grows from impact face (small) toward back face (large).

Interfaces between plies of the same direction (``dtheta = 0``, e.g. the
doubled 90/90 ply at the mid-plane of a symmetric laminate) carry no
delamination: with no bending-stiffness mismatch there is nothing to drive
one (Liu 1988), and delaminations are observed only at interfaces between
plies of different orientation (Abrate 1991, *Appl. Mech. Rev.* 44, 155-190). A layup whose plies all share
one direction keeps an ellipse at every interface (circles), so it still
receives the target damage area.

Total DPA is enforced by a Brent root-find on a single scalar multiplier so the
polygon-union footprint equals the target DPA within 1%. A target below
0.01 mm^2 (``_MIN_TARGET_DPA_MM2``, under any C-scan's resolution) gives no
delamination, like an impact at or below the onset energy.

With ``max_half_extent_mm`` each ellipse is shrunk (aspect ratio kept) until
its bounding box fits inside the given half-extents about the centroid: the
delaminations stay inside the window the plate was supported on when it was
impacted. The union of differently oriented ellipses spreads well beyond the
circle of the same area (on the 40-ply hard Lovejoy & Scotti coupon, 125 mm
across a 102 mm coupon, which would leave fe3d a damage band through the
whole section). If the window cannot hold the target area, the saturated footprint
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


def _ply_mismatch_deg(theta_a_deg, theta_b_deg):
    """Acute angle between two fibre directions, in [0, 90] deg.

    Works on scalars and NumPy arrays; fibre directions are axial, so the
    difference is taken modulo 180 deg.
    """
    d = np.mod(np.abs(np.asarray(theta_b_deg, dtype=float) - theta_a_deg), 180.0)
    return np.minimum(d, 180.0 - d)


def _axial_angle_deg(theta_deg):
    """Fibre angle mapped into (-90, 90] deg (``theta`` and ``theta + 180``
    are the same direction). Works on scalars and NumPy arrays."""
    return 90.0 - np.mod(90.0 - np.asarray(theta_deg, dtype=float), 180.0)


def _aspect_ratio(delta_theta_deg: float) -> float:
    return min(4.0, max(1.0, 1.0 + 0.025 * float(_ply_mismatch_deg(0.0, delta_theta_deg))))


def _orientation_deg(impact_side_ply_deg: float, back_side_ply_deg: float) -> float:
    """Ellipse orientation of the interface between the two plies: the fibre
    angle of the back-face-side ply (Liu 1988)."""
    del impact_side_ply_deg  # the impact-face-side ply does not set the axis
    return float(_axial_angle_deg(back_side_ply_deg))


def _relative_size(interface_index: int, n_interfaces: int) -> float:
    """Weight grows from 0.3 near impact face to 1.0 near back face."""
    z = (interface_index + 1) / n_interfaces  # 0 < z <= 1
    return 0.3 + 0.7 * z


def _template_arrays(
    layup_deg: List[float], n_interfaces: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute per-interface (ar, rel, orient) as NumPy arrays."""
    layup = np.asarray(layup_deg, dtype=float)
    impact_side = layup[:-1]
    back_side = layup[1:]
    dtheta = _ply_mismatch_deg(impact_side, back_side)
    ar = np.clip(1.0 + 0.025 * dtheta, 1.0, 4.0)
    orient = _axial_angle_deg(back_side)
    idx = np.arange(n_interfaces, dtype=float)
    z = (idx + 1.0) / n_interfaces
    rel = 0.3 + 0.7 * z
    return ar, rel, orient


def _delaminated_interfaces(layup_deg: List[float]) -> np.ndarray:
    """Indices of the interfaces that carry a delamination: those between
    plies of different fibre direction, or every interface when all plies
    share one direction."""
    layup = np.asarray(layup_deg, dtype=float)
    mismatch = _ply_mismatch_deg(layup[:-1], layup[1:])
    active = np.flatnonzero(mismatch > _SAME_DIRECTION_TOL_DEG)
    return active if active.size else np.arange(len(layup) - 1)


# Plies whose fibre directions differ by less than this are the same direction.
_SAME_DIRECTION_TOL_DEG = 1e-6

# A target DPA below this (mm^2) is below the resolution of any C-scan and
# of the 72-sided polygons the union is built from: it is no damage, and
# ``distribute_damage`` returns no delamination for it (as impact_to_damage
# does for an impact at or below the onset energy).
_MIN_TARGET_DPA_MM2 = 0.01
# Lower limit (mm) of the scalar bracket's downward extension; at 1e-4 mm
# the largest ellipse of any layup is far below ``_MIN_TARGET_DPA_MM2``.
_MIN_SCALAR = 1e-4


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
    if n_interfaces <= 0 or target_dpa_mm2 < _MIN_TARGET_DPA_MM2:
        return []

    ar_arr, rel_arr, orient_arr = _template_arrays(layup_deg, n_interfaces)
    # Same-direction interfaces carry no delamination (module docstring).
    active = _delaminated_interfaces(layup_deg)
    ar_arr, rel_arr, orient_arr = ar_arr[active], rel_arr[active], orient_arr[active]
    index_list = active.tolist()
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
                major_mm=majors[k],
                minor_mm=minors[k],
                orientation_deg=orient_list[k],
            )
            for k, i in enumerate(index_list)
        ]

    def _union_area(scalar: float) -> float:
        return DamageState(
            _build(scalar), dent_depth_mm, fiber_break_radius_mm
        ).projected_damage_area_mm2

    lo, hi = 0.1, 50.0
    # A target just above the resolution floor (an impact a few mJ above the
    # onset energy) sits below the union area at the default lower end, so
    # the bracket extends downward too (#163).
    while _union_area(lo) > target_dpa_mm2:
        lo *= 0.5
        if lo < _MIN_SCALAR:
            return []
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
