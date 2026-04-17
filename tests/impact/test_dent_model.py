from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.dent_model import dent_depth_mm, fiber_break_radius_mm


def test_zero_below_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert dent_depth_mm(m, E_impact_J=1.0, E_onset_J=5.0, h_mm=2.0) == 0.0


def test_zero_at_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert dent_depth_mm(m, E_impact_J=5.0, E_onset_J=5.0, h_mm=2.0) == 0.0


def test_monotonic_in_energy():
    # Use energies just above onset so the raw dent is well below the 0.5*h cap.
    # With h_mm=2.0, G_Ic=0.28, beta=0.05, gamma=0.5:
    #   raw(E=5.001) << 1.0 = 0.5*h  => no cap kicks in
    m = MATERIAL_LIBRARY["IM7/8552"]
    d1 = dent_depth_mm(m, E_impact_J=5.001, E_onset_J=5.0, h_mm=2.0)
    d2 = dent_depth_mm(m, E_impact_J=5.01, E_onset_J=5.0, h_mm=2.0)
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


def test_dent_depth_capped_at_half_thickness():
    m = MATERIAL_LIBRARY["IM7/8552"]
    h = 1.0  # mm
    # Very high energy that would otherwise produce a huge dent
    d = dent_depth_mm(m, E_impact_J=10_000.0, E_onset_J=1.0, h_mm=h)
    assert d <= 0.5 * h + 1e-9
