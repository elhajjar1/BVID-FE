import warnings

import pytest

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.mapping import (
    ImpactEvent,
    SmallMassQuasiStaticWarning,
    _diameter_dpa_factor,
    _dynamic_amplification_factor,
    _shape_dpa_factor,
    impact_to_damage,
)
from bvidfe.impact.olsson import onset_energy


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


@pytest.mark.parametrize("dE_J", [1e-2, 1e-3, 1e-4, 1e-6])
def test_impact_just_above_onset_does_not_raise(dE_J):
    """Issue #163: an impact a few mJ above the onset energy targets a DPA
    under the union area at distribute_damage's default lower bracket end,
    and brentq's ``f(a) and f(b) must have different signs`` escaped from
    impact_to_damage. The damage is now either met within 1% or, below the
    0.01 mm^2 resolution floor, empty."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    E_onset = onset_energy(lam, pan, imp)
    ds = impact_to_damage(ImpactEvent(energy_J=E_onset + dE_J, impactor=imp, mass_kg=5.5), lam, pan)
    assert ds.dent_depth_mm >= 0.0
    dpa = ds.projected_damage_area_mm2
    assert dpa >= 0.0
    # The energy-balance target with the same mass, shape and diameter factors
    # impact_to_damage applies.
    dpa_target = m.olsson_alpha * dE_J * 1e3 / m.G_IIc
    dpa_target *= _dynamic_amplification_factor(5.5, lam, pan)[0]
    dpa_target *= _shape_dpa_factor(imp.shape) * _diameter_dpa_factor(imp.diameter_mm, imp.shape)
    if dpa_target >= 0.01:
        assert ds.delaminations
        assert abs(dpa - dpa_target) / dpa_target < 0.01  # distribute_damage's contract
    else:
        assert ds.delaminations == []
        assert dpa == 0.0


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
    E_onset = onset_energy(lam, pan, ev.impactor)
    dpa_raw = m.olsson_alpha * (ev.energy_J - E_onset) * 1e3 / m.G_IIc
    A_cap = 0.8 * pan.Lx_mm * pan.Ly_mm
    dpa_target = min(dpa_raw, A_cap)
    assert dpa_target > 0
    assert abs(ds.projected_damage_area_mm2 - dpa_target) / dpa_target < 0.01


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


def test_small_mass_emits_quasi_static_warning():
    """Issue #17: a sub-unity impactor/plate mass ratio must emit the
    Olsson quasi-static-validity UserWarning. Regressing the m_ratio < 1.0
    guard (or its message) would silently hide a documented 30%+ under-
    prediction in drop-tower / light-impactor simulations."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0] * 16, 0.2)
    pan = PanelGeometry(500, 500)
    # Above the onset energy (about 6.8 J on this panel), where the damage
    # area and the warning are computed.
    ev = ImpactEvent(energy_J=20.0, impactor=ImpactorGeometry(), mass_kg=0.5)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        impact_to_damage(ev, lam, pan)
    matching = [
        m
        for m in w
        if issubclass(m.category, SmallMassQuasiStaticWarning)
        and "quasi-static" in str(m.message)
        and "mass" in str(m.message)
    ]
    assert matching, [(m.category.__name__, str(m.message)) for m in w]
