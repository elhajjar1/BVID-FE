import warnings

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
    import warnings

    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    pan = PanelGeometry(150, 100)
    ev = ImpactEvent(energy_J=30.0, impactor=ImpactorGeometry(), mass_kg=5.5)
    # At 30 J the raw Olsson DPA exceeds 80% of the panel; capture the warning
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        ds = impact_to_damage(ev, lam, pan)
    # impact_to_damage caps DPA at 80% of panel area when Olsson would exceed it;
    # distribute_damage enforces union ≈ DPA_target (possibly capped) within 1%.
    from bvidfe.impact.olsson import onset_energy

    E_onset = onset_energy(lam, pan, ev.impactor)
    h = lam.thickness_mm
    dpa_raw = m.olsson_alpha * (ev.energy_J - E_onset) * 1e3 / (m.G_IIc * h)
    A_cap = 0.8 * pan.Lx_mm * pan.Ly_mm
    dpa_target = min(dpa_raw, A_cap)
    assert dpa_target > 0
    assert abs(ds.projected_damage_area_mm2 - dpa_target) / dpa_target < 0.01


def test_small_mass_emits_quasi_static_warning():
    """A light impactor on a heavy plate (m_ratio < 1) must emit the Olsson
    quasi-static-validity UserWarning.

    This safety-rail has no other test: the existing DAF / SummaryTab tests
    pin DAF == 1.0 at the 5.5 kg reference and a derived GUI note, neither of
    which depends on this specific warning. Flipping the ``m_ratio < 1.0``
    guard (or narrowing it) would slip through silently and hide the 30%+
    underprediction documented in the docstring.
    """
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0] * 16, 0.2)
    pan = PanelGeometry(500, 500)
    ev = ImpactEvent(energy_J=20.0, impactor=ImpactorGeometry(), mass_kg=0.5)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        impact_to_damage(ev, lam, pan)
    messages = [str(msg.message) for msg in w]
    assert any(
        "quasi-static" in msg and "mass" in msg for msg in messages
    ), messages


def test_impact_to_damage_clips_large_dpa():
    """DPA is clipped to 80% of panel area for high-energy impacts on small panels."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 2, 0.152)
    pan = PanelGeometry(40, 30)  # small panel -> easy to exceed
    ev = ImpactEvent(100.0, ImpactorGeometry(), mass_kg=5.5)  # lots of energy

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ds = impact_to_damage(ev, lam, pan)

    # DPA should be clipped to 0.8 * panel area
    A_panel = 40 * 30
    assert ds.projected_damage_area_mm2 <= 0.8 * A_panel * 1.02  # 2% slack for union
    # And we should have emitted at least one UserWarning
    assert any("exceeds 80%" in str(msg.message) for msg in w)
