from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage


def test_below_threshold_returns_empty_damage():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    pan = PanelGeometry(150, 100)
    ev = ImpactEvent(energy_J=0.01, impactor=ImpactorGeometry(), mass_kg=5.5)
    ds = impact_to_damage(ev, lam, pan)
    assert ds.delaminations == []
    assert ds.dent_depth_mm == 0.0
    assert ds.projected_damage_area_mm2 == 0.0


def test_above_threshold_produces_n_minus_1_ellipses():
    m = MATERIAL_LIBRARY["IM7/8552"]
    layup = [0, 45, -45, 90] * 4
    lam = Laminate(m, layup, 0.152)
    pan = PanelGeometry(150, 100)
    ev = ImpactEvent(energy_J=30.0, impactor=ImpactorGeometry(), mass_kg=5.5)
    ds = impact_to_damage(ev, lam, pan)
    assert len({e.interface_index for e in ds.delaminations}) == len(layup) - 1
    assert ds.dent_depth_mm > 0
    assert ds.projected_damage_area_mm2 > 0


def test_impact_event_location_defaults_to_panel_center():
    ev = ImpactEvent(energy_J=10.0, impactor=ImpactorGeometry(), mass_kg=5.5)
    assert ev.location_xy_mm == (0.0, 0.0)


def test_dpa_conserved_after_distribution():
    """The full path: Olsson → DPA_target → distribute_damage → union area matches target."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    pan = PanelGeometry(150, 100)
    ev = ImpactEvent(energy_J=30.0, impactor=ImpactorGeometry(), mass_kg=5.5)
    ds = impact_to_damage(ev, lam, pan)
    # impact_to_damage uses Olsson alpha * (E - Eonset) / (G_IIc * h) for DPA_target;
    # distribute_damage enforces union ≈ DPA_target within 1%.
    # Check union equals the target DPA implied by Olsson (recompute independently):
    from bvidfe.impact.olsson import onset_energy

    E_onset = onset_energy(lam, pan, ev.impactor)
    h = lam.thickness_mm
    dpa_target = m.olsson_alpha * (ev.energy_J - E_onset) * 1e3 / (m.G_IIc * h)
    assert dpa_target > 0
    assert abs(ds.projected_damage_area_mm2 - dpa_target) / dpa_target < 0.01
