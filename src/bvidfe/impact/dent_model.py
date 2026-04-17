"""Empirical dent depth and fiber-break core radius from impact energy."""

from __future__ import annotations

import math

from bvidfe.core.material import OrthotropicMaterial


def dent_depth_mm(
    m: OrthotropicMaterial,
    E_impact_J: float,
    E_onset_J: float,
    h_mm: float,
) -> float:
    """Residual dent depth on the impact face, mm.

    d / h = beta * ((E_impact - E_onset) / (G_Ic * h^2))^gamma

    Zero at or below the damage-onset energy. `G_Ic` is in N/mm; energy inputs
    are in J; multiply (E_impact - E_onset) by 1e3 to get N*mm before the division.
    """
    if E_impact_J <= E_onset_J:
        return 0.0
    numerator = (E_impact_J - E_onset_J) * 1e3  # N*mm
    denom = m.G_Ic * h_mm**2  # N*mm
    raw = h_mm * m.dent_beta * (numerator / denom) ** m.dent_gamma
    # Physical cap: dent cannot exceed 50% of laminate thickness
    # (beyond that it's perforation, not BVID)
    return min(raw, 0.5 * h_mm)


def fiber_break_radius_mm(m: OrthotropicMaterial, E_impact_J: float) -> float:
    """Radius of the fiber-breakage core on the back face (mm).

    r_fb = eta * sqrt(max(0, E_impact - E_fb_threshold))

    Defaults to 0 unless the material overrides `fiber_break_eta`.
    """
    excess = max(0.0, E_impact_J - m.fiber_break_E_threshold)
    return m.fiber_break_eta * math.sqrt(excess)
