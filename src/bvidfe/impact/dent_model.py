"""Empirical dent depth and fiber-break core radius from impact energy."""

from __future__ import annotations

import math

from bvidfe.core.material import OrthotropicMaterial

# Dent depth of a flat-ended or 60 deg conical impactor relative to a
# hemispherical one at the same energy above onset. Fitted to Hasebe et al.
# (2025), T800S/3900-2B: the cone's median measured / hemispherical-law ratio
# is 8.6 (92 specimens; with the membrane-stiffened onset energy the law as
# used predicts the cone dents at a median 0.91x); 34 of the 40 flat-ended
# dents are recorded as 0, and 0.03 matches their mean (0.0034 mm).
_DENT_SHAPE_FACTOR = {"hemispherical": 1.0, "conical": 8.5, "flat": 0.03}


def dent_depth_mm(
    m: OrthotropicMaterial,
    E_impact_J: float,
    E_onset_J: float,
    h_mm: float,
    shape: str = "hemispherical",
) -> float:
    """Residual dent depth on the impact face, mm.

    d = dent_beta * (E_impact - E_onset) ** dent_gamma * S_shape, E in J,

    capped at ``h / 2`` and zero at or below the damage-onset energy.
    ``S_shape`` is 1 for a hemispherical impactor, 8.5 for the 60 deg cone
    and 0.03 for a flat-ended one (``_DENT_SHAPE_FACTOR``). The default
    constants (0.044 mm/J^0.57, 0.57) are fitted to the measured dents of
    Lovejoy & Scotti (IM7/8552) and Hasebe et al. (T800S/3900-2B); see
    ``docs/physics_models.md`` and ``validation/README.md``.
    """
    if E_impact_J <= E_onset_J:
        return 0.0
    raw = m.dent_beta * (E_impact_J - E_onset_J) ** m.dent_gamma
    raw *= _DENT_SHAPE_FACTOR.get(shape, 1.0)
    # Safety bound: a dent deeper than half the laminate is perforation, not
    # BVID. Of the measured cases it binds only on seven 8-ply laminates
    # struck by the cone at 5.2 J, whose measured dents exceed h / 2 too.
    return min(raw, 0.5 * h_mm)


def fiber_break_radius_mm(m: OrthotropicMaterial, E_impact_J: float) -> float:
    """Radius of the fiber-breakage core on the back face (mm).

    r_fb = eta * sqrt(max(0, E_impact - E_fb_threshold))

    Defaults to 0 unless the material overrides `fiber_break_eta`.
    """
    excess = max(0.0, E_impact_J - m.fiber_break_E_threshold)
    return m.fiber_break_eta * math.sqrt(excess)
