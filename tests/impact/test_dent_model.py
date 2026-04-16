from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.dent_model import dent_depth_mm, fiber_break_radius_mm


def test_zero_below_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert dent_depth_mm(m, E_impact_J=1.0, E_onset_J=5.0, h_mm=2.0) == 0.0


def test_zero_at_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert dent_depth_mm(m, E_impact_J=5.0, E_onset_J=5.0, h_mm=2.0) == 0.0


def test_monotonic_in_energy():
    m = MATERIAL_LIBRARY["IM7/8552"]
    d1 = dent_depth_mm(m, E_impact_J=10.0, E_onset_J=5.0, h_mm=2.0)
    d2 = dent_depth_mm(m, E_impact_J=30.0, E_onset_J=5.0, h_mm=2.0)
    assert d2 > d1 > 0


def test_fiber_break_radius_zero_when_eta_zero():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert m.fiber_break_eta == 0.0
    assert fiber_break_radius_mm(m, E_impact_J=100.0) == 0.0


def test_fiber_break_radius_zero_below_threshold():
    import dataclasses

    m = MATERIAL_LIBRARY["IM7/8552"]
    m2 = dataclasses.replace(m, fiber_break_eta=0.1, fiber_break_E_threshold=50.0)
    # Below threshold
    assert fiber_break_radius_mm(m2, E_impact_J=10.0) == 0.0
    # Above threshold: positive
    assert fiber_break_radius_mm(m2, E_impact_J=100.0) > 0
