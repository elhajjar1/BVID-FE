"""Forward orchestrator: ImpactEvent -> DamageState.

Pipeline:
  1. compute Olsson onset energy (boundary- and shape-aware)
  2. if below threshold -> empty DamageState
  3. compute target DPA = alpha * (E - Eonset) * 1e3 / G_IIc, then apply
     multiplicative corrections for:
       - impactor mass, (5.5 kg / mass_kg) ** 0.1 (a warning below the
         plate's mass)
       - impactor-tip footprint (hemispherical / flat / conical)
       - impactor-diameter spread
     and cap it at 80% of the panel area (``DPACapClipWarning``)
  4. compute dent depth and fiber-break radius
  5. distribute DPA across interfaces via peanut template with union-area
     scaling, every ellipse kept inside the impacted panel
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Tuple

from bvidfe._validation import require_non_negative, require_positive, require_xy
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.damage.state import DamageState
from bvidfe.impact.dent_model import dent_depth_mm, fiber_break_radius_mm
from bvidfe.impact.olsson import onset_energy
from bvidfe.impact.shape_templates import distribute_damage


class SmallMassQuasiStaticWarning(UserWarning):
    """The impactor is comparable to or lighter than the plate, so Olsson's
    quasi-static threshold model may underpredict damage."""


class DPACapClipWarning(UserWarning):
    """The Olsson-predicted DPA exceeded the panel-area cap and was clipped."""


# Impactor-shape multiplier on target DPA, relative to a 16 mm hemisphere.
# Fitted to Hasebe et al. (Data in Brief 60, 2025), the only open dataset
# that varies the impactor on one material (T800S/3900-2B): over same-layup,
# same-energy groups, a 15.9 mm flat-ended impactor left 0.20x the C-scan
# area of a 12.7 mm hemisphere (14 groups; 11 of 40 flat specimens had
# none) and a 60 deg cone 1.20x (37 groups). Relative to a 16 mm hemisphere
# that is 0.2 and 1.3, the values used while the onset energy was ~10x too
# low. With the corrected threshold load the onset energy depends on the
# tip through the contact stiffness (the linearised cone contact is the
# softest, so the cone's onset is the highest), and the factors are refitted
# so that, per same-layup, same-energy group, the model's predicted /
# measured area ratio for each tip matches the 12.7 mm hemisphere's
# (median over groups): flat 0.21 (kept at 0.2) and conical 1.54 (1.5).
# With them the per-impactor median predicted / measured area is 1.00
# (flat), 1.02 (cone), 0.99, 1.05 and 0.96 (12.7, 25.4 and 31.8 mm
# hemispheres). The previous, unsourced 1.4 (flat) and 0.7 (conical) ran
# the other way.
_SHAPE_DPA_FACTOR: dict[str, float] = {
    "hemispherical": 1.0,
    "flat": 0.2,
    "conical": 1.5,
}

# Reference impactor diameter (mm) for the footprint-spread correction.
# Matches ASTM D7136 default. The exponent is fitted to the same data:
# 25.4 and 31.8 mm hemispheres left 0.88x and 0.87x the C-scan area of the
# 12.7 mm one (30 and 6 groups), against 0.81x and 0.76x at the previous
# 0.3. The corrected onset energy moved the model's same-group ratios to
# the 12.7 mm hemisphere by under 5% (25.4 mm: 1.01 to 1.04; 31.8 mm: 1.25
# to 1.27, on 6 cross-ply groups only), so 0.2 is kept.
_DIAMETER_REF_MM: float = 16.0
_DIAMETER_DPA_EXPONENT: float = 0.2


def _shape_dpa_factor(shape: str) -> float:
    return _SHAPE_DPA_FACTOR.get(shape, 1.0)


def _diameter_dpa_factor(diameter_mm: float, shape: str = "hemispherical") -> float:
    """Smaller impactors concentrate damage (larger factor); larger
    impactors spread it (smaller factor). Mild (exponent 0.2). A cone's
    contact is set by its angle, not its shank diameter, so the factor is
    1 for ``conical``."""
    if diameter_mm <= 0 or shape == "conical":
        return 1.0
    return (_DIAMETER_REF_MM / diameter_mm) ** _DIAMETER_DPA_EXPONENT


# Reference impactor mass (kg) for the small-mass amplification correction.
# 5.5 kg is the ASTM D7136 drop-weight standard and the calibration baseline
# for the empirical knockdown formulas in this package — at this mass the
# DPA-mass correction is exactly unity, so historical test data and the
# validation gate remain unshifted.
_MASS_REF_KG: float = 5.5
_MASS_DPA_EXPONENT: float = 0.10  # very mild: (5.5/m)**0.1


def _dynamic_amplification_factor(
    mass_kg: float, lam: Laminate, panel: PanelGeometry
) -> tuple[float, float]:
    """Mass-ratio-dependent amplification of DPA.

    Returns ``(daf, m_ratio)`` where ``m_ratio = m_impactor / m_plate``
    (both in kg) and::

        daf = (m_ref / m_impactor) ** 0.1        # mild, centered on m_ref = 5.5 kg

    At the ASTM D7136 reference mass (5.5 kg) DAF = 1.000 exactly, so the
    calibration of the empirical knockdown formulas is unperturbed. Lighter
    impactors at the same energy carry higher velocity / shorter contact
    time → slightly more delamination per Joule; heavier impactors → slightly
    less. The exponent 0.1 keeps the effect bounded to about [0.85, 1.17]
    across the practical 1 kg .. 20 kg range.

    We also return the mass ratio ``m_impactor / m_plate_eff`` so the caller
    can emit a "quasi-static validity" warning if the impactor is lighter
    than the plate (m_ratio < 1 — the small-mass regime where the Olsson
    quasi-static threshold can underpredict damage by 30%+).

    ``OrthotropicMaterial.rho`` is in kg/mm^3 (e.g. 1.57e-6 for CFRP,
    matching 1.57 g/cm^3); the m_ratio below is therefore in consistent kg.
    """
    rho_kg_per_mm3 = lam.material.rho
    m_plate_eff_kg = rho_kg_per_mm3 * panel.Lx_mm * panel.Ly_mm * lam.thickness_mm
    if m_plate_eff_kg <= 0 or mass_kg <= 0:
        return 1.0, float("inf")
    m_ratio = mass_kg / m_plate_eff_kg
    daf = (_MASS_REF_KG / mass_kg) ** _MASS_DPA_EXPONENT
    return daf, m_ratio


@dataclass
class ImpactEvent:
    """Parameters describing a single impact event.

    Attributes
    ----------
    energy_J : float
        Impact energy in J, finite and ``>= 0``. At or below the onset
        energy the damage state is empty.
    impactor : ImpactorGeometry
        Impactor tip; default the ASTM D7136 16 mm hemisphere.
    mass_kg : float
        Impactor mass in kg, finite and ``> 0`` (default 5.5, the ASTM
        D7136 reference, at which the mass correction is 1). The damage area scales by
        ``(5.5 / mass_kg) ** 0.1``.
    location_xy_mm : tuple of float
        Impact point ``(x, y)`` in mm on the impacted panel
        (``AnalysisConfig.impact_panel`` if set, else ``panel``). The
        default ``(0.0, 0.0)`` means the panel centre. A list (as JSON
        gives) is stored as a tuple.
    """

    energy_J: float
    impactor: ImpactorGeometry = field(default_factory=ImpactorGeometry)
    mass_kg: float = 5.5
    location_xy_mm: Tuple[float, float] = (0.0, 0.0)

    def __post_init__(self) -> None:
        require_non_negative(self.energy_J, "ImpactEvent.energy_J")
        require_positive(self.mass_kg, "ImpactEvent.mass_kg")
        if not isinstance(self.impactor, ImpactorGeometry):
            raise ValueError(
                f"ImpactEvent.impactor must be an ImpactorGeometry "
                f"(got {type(self.impactor).__name__})"
            )
        self.location_xy_mm = require_xy(self.location_xy_mm, "ImpactEvent.location_xy_mm")


def impact_to_damage(event: ImpactEvent, lam: Laminate, panel: PanelGeometry) -> DamageState:
    """Map an impact event to a full BVID damage state.

    Follows the pipeline in the module docstring. ``panel`` is the impacted
    panel: its size and ``boundary`` set the onset energy, and the
    delaminations are kept inside it, measured from the impact point.

    Returns a :class:`DamageState` in ``panel``'s frame (empty at or below
    the onset energy).

    Warns
    -----
    SmallMassQuasiStaticWarning
        Impactor lighter than the plate (outside Olsson's quasi-static regime).
    DPACapClipWarning
        Predicted damage area above 80% of the panel area, clipped to it.
    UserWarning
        The damage area does not fit inside the panel around the impact
        point; the largest footprint that fits is used. Also for a
        ``"free"`` panel boundary.

    Examples
    --------
    >>> from bvidfe import MATERIAL_LIBRARY, PanelGeometry
    >>> from bvidfe.core.laminate import Laminate
    >>> lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    >>> d = impact_to_damage(ImpactEvent(energy_J=10.0), lam, PanelGeometry(150, 100))
    >>> round(d.projected_damage_area_mm2)
    308
    """
    material = lam.material

    # If location is (0, 0), interpret as panel center for Olsson bending stiffness
    loc = event.location_xy_mm
    if loc == (0.0, 0.0):
        loc = (panel.Lx_mm / 2, panel.Ly_mm / 2)

    E_onset = onset_energy(lam, panel, event.impactor, location_xy_mm=loc)
    if event.energy_J <= E_onset:
        return DamageState([], dent_depth_mm=0.0, fiber_break_radius_mm=0.0)

    # Energy balance: a fraction alpha of the impact energy above onset goes
    # into mode-II delamination, DPA = alpha (E - E_onset) / G_IIc (mm^2;
    # energies J -> N*mm via *1e3). It used to be divided by the laminate
    # thickness as well, which left DPA in mm and flat in thickness at
    # constant energy per thickness, where measured areas grow with it:
    # Hasebe et al. (2025), 8 to 24 plies at fixed J/mm, and Lovejoy &
    # Scotti (2019), whose 7.6 mm coupons at the same 8.9 J/mm as their
    # 6.1 mm ones had 1.21x the C-scan area.
    h = lam.thickness_mm
    dpa_target = material.olsson_alpha * (event.energy_J - E_onset) * 1e3 / material.G_IIc

    # ---- multiplicative physics-motivated corrections --------------------
    # Dynamic amplification (small impactor mass => more damage per Joule)
    daf, m_ratio = _dynamic_amplification_factor(event.mass_kg, lam, panel)
    # Olsson's quasi-static threshold assumes the impactor is much heavier
    # than the plate (canonical rule: m_ratio > 40). Below ~1, dynamic
    # effects become significant and DPA is likely underpredicted regardless
    # of our mild DAF correction.
    if m_ratio < 1.0:
        warnings.warn(
            f"Impactor mass ({event.mass_kg:.2f} kg) is comparable to or lighter "
            f"than the plate's effective mass (ratio = {m_ratio:.2f}); Olsson's "
            f"quasi-static threshold model may underpredict damage in this "
            f"regime. Predictions should be interpreted qualitatively.",
            SmallMassQuasiStaticWarning,
            stacklevel=2,
        )
    # Impactor tip shape (flat-ended: less delamination; conical: more)
    shape_factor = _shape_dpa_factor(event.impactor.shape)
    # Impactor diameter spread (smaller D -> more concentrated -> larger DPA factor)
    diam_factor = _diameter_dpa_factor(event.impactor.diameter_mm, event.impactor.shape)

    dpa_target *= daf * shape_factor * diam_factor

    if dpa_target <= 0:
        return DamageState([], dent_depth_mm=0.0, fiber_break_radius_mm=0.0)

    # Cap DPA at 80% of panel area to avoid physically unreasonable footprints
    A_panel = panel.Lx_mm * panel.Ly_mm
    A_cap = 0.8 * A_panel
    if dpa_target > A_cap:
        warnings.warn(
            f"Olsson-predicted DPA ({dpa_target:.0f} mm^2) exceeds 80% of panel "
            f"area ({A_cap:.0f} mm^2). Clipping to {A_cap:.0f} mm^2. The impact "
            f"energy may exceed the panel's capacity to contain damage without "
            f"edge effects; consider a larger panel or lower energy.",
            DPACapClipWarning,
            stacklevel=2,
        )
        dpa_target = A_cap

    dent = dent_depth_mm(material, event.energy_J, E_onset, h, shape=event.impactor.shape)
    r_fb = fiber_break_radius_mm(material, event.energy_J)

    # Centroid for distribute_damage is the event location on the panel. The
    # delaminations stay inside the panel (support window) that was impacted.
    half_extent = (
        min(loc[0], panel.Lx_mm - loc[0]),
        min(loc[1], panel.Ly_mm - loc[1]),
    )
    ellipses = distribute_damage(
        layup_deg=lam.layup_deg,
        target_dpa_mm2=dpa_target,
        dent_depth_mm=dent,
        fiber_break_radius_mm=r_fb,
        centroid_mm=loc,
        max_half_extent_mm=half_extent,
    )

    return DamageState(delaminations=ellipses, dent_depth_mm=dent, fiber_break_radius_mm=r_fb)
