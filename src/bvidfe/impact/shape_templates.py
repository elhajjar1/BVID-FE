"""Layup-dependent 'peanut' template distributing DPA into per-interface delamination ellipses.

Each interface's ellipse is parametrised by:
- aspect ratio AR = 1 + 0.025 * |delta_theta|, clipped to [1, 4]
- orientation = bisector of neighbour-ply angles
- relative size grows from impact face (small) toward back face (large)

Total DPA is enforced by a Brent root-find on a single scalar multiplier so the
polygon-union footprint equals the target DPA within 1%.
"""

from __future__ import annotations

from typing import List, Tuple

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


def distribute_damage(
    layup_deg: List[float],
    target_dpa_mm2: float,
    dent_depth_mm: float,
    fiber_break_radius_mm: float,
    centroid_mm: Tuple[float, float] = (0.0, 0.0),
) -> List[DelaminationEllipse]:
    n_plies = len(layup_deg)
    n_interfaces = n_plies - 1
    if n_interfaces <= 0 or target_dpa_mm2 <= 0:
        return []

    templates = []
    for i in range(n_interfaces):
        dtheta = layup_deg[i + 1] - layup_deg[i]
        templates.append(
            {
                "i": i,
                "ar": _aspect_ratio(dtheta),
                "orient": _orientation_deg(layup_deg[i], layup_deg[i + 1]),
                "rel": _relative_size(i, n_interfaces),
            }
        )

    def _build(scalar: float) -> List[DelaminationEllipse]:
        return [
            DelaminationEllipse(
                interface_index=t["i"],
                centroid_mm=centroid_mm,
                major_mm=scalar * t["ar"] * t["rel"],
                minor_mm=scalar * t["rel"],
                orientation_deg=t["orient"],
            )
            for t in templates
        ]

    def _union_area(scalar: float) -> float:
        return DamageState(
            _build(scalar), dent_depth_mm, fiber_break_radius_mm
        ).projected_damage_area_mm2

    lo, hi = 0.1, 50.0
    while _union_area(hi) < target_dpa_mm2:
        hi *= 2
        if hi > 1e6:
            raise RuntimeError("shape_templates: cannot bracket target DPA")

    scalar = brentq(lambda s: _union_area(s) - target_dpa_mm2, lo, hi, xtol=1e-3)
    return _build(scalar)
