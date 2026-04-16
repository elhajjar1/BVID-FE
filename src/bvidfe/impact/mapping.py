"""Forward orchestrator: ImpactEvent -> DamageState.

Pipeline:
  1. compute Olsson onset energy
  2. if below threshold -> empty DamageState
  3. compute target DPA = alpha * (E - Eonset) * 1e3 / (G_IIc * h)
  4. compute dent depth and fiber-break radius
  5. distribute DPA across interfaces via peanut template with union-area scaling
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.damage.state import DamageState
from bvidfe.impact.dent_model import dent_depth_mm, fiber_break_radius_mm
from bvidfe.impact.olsson import onset_energy
from bvidfe.impact.shape_templates import distribute_damage


@dataclass
class ImpactEvent:
    """Parameters describing a single impact event."""

    energy_J: float
    impactor: ImpactorGeometry = field(default_factory=ImpactorGeometry)
    mass_kg: float = 5.5
    location_xy_mm: Tuple[float, float] = (0.0, 0.0)


def impact_to_damage(event: ImpactEvent, lam: Laminate, panel: PanelGeometry) -> DamageState:
    """Map an impact event to a full BVID damage state."""
    material = lam.material

    # If location is (0, 0), interpret as panel center for Olsson bending stiffness
    loc = event.location_xy_mm
    if loc == (0.0, 0.0):
        loc = (panel.Lx_mm / 2, panel.Ly_mm / 2)

    E_onset = onset_energy(lam, panel, event.impactor, location_xy_mm=loc)
    if event.energy_J <= E_onset:
        return DamageState([], dent_depth_mm=0.0, fiber_break_radius_mm=0.0)

    # Olsson: DPA = alpha * (E - Eonset) / (G_IIc * h)  (SI mm/N)
    # Energies in J -> N*mm via *1e3
    h = lam.thickness_mm
    dpa_target = material.olsson_alpha * (event.energy_J - E_onset) * 1e3 / (material.G_IIc * h)
    if dpa_target <= 0:
        return DamageState([], dent_depth_mm=0.0, fiber_break_radius_mm=0.0)

    dent = dent_depth_mm(material, event.energy_J, E_onset, h)
    r_fb = fiber_break_radius_mm(material, event.energy_J)

    # Centroid for distribute_damage is the event location on the panel
    ellipses = distribute_damage(
        layup_deg=lam.layup_deg,
        target_dpa_mm2=dpa_target,
        dent_depth_mm=dent,
        fiber_break_radius_mm=r_fb,
        centroid_mm=loc,
    )

    return DamageState(delaminations=ellipses, dent_depth_mm=dent, fiber_break_radius_mm=r_fb)
