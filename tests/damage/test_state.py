import math

import pytest

from bvidfe.damage.state import DamageState, DelaminationEllipse


def test_ellipse_area():
    e = DelaminationEllipse(
        interface_index=0,
        centroid_mm=(0, 0),
        major_mm=10,
        minor_mm=5,
        orientation_deg=0,
    )
    assert math.isclose(e.area_mm2, math.pi * 10 * 5)


def test_ellipse_rejects_nonpositive_axes():
    with pytest.raises(ValueError):
        DelaminationEllipse(0, (0, 0), 0, 5, 0)
    with pytest.raises(ValueError):
        DelaminationEllipse(0, (0, 0), 10, -1, 0)


def test_ellipse_rejects_negative_interface_index():
    with pytest.raises(ValueError):
        DelaminationEllipse(-1, (0, 0), 10, 5, 0)


def test_damage_state_projected_area_two_identical_ellipses_union():
    e = DelaminationEllipse(0, (0, 0), 10, 5, 0)
    e2 = DelaminationEllipse(1, (0, 0), 10, 5, 0)  # same footprint, different interface
    ds = DamageState([e, e2], dent_depth_mm=0.3)
    # Union should equal single-ellipse area (they overlap 100% in plan view)
    assert math.isclose(ds.projected_damage_area_mm2, math.pi * 10 * 5, rel_tol=1e-2)


def test_damage_state_projected_area_two_non_overlapping():
    e1 = DelaminationEllipse(0, (0, 0), 10, 5, 0)
    e2 = DelaminationEllipse(1, (100, 100), 10, 5, 0)  # far apart
    ds = DamageState([e1, e2], dent_depth_mm=0.3)
    assert math.isclose(ds.projected_damage_area_mm2, 2 * math.pi * 10 * 5, rel_tol=1e-2)


def test_damage_state_empty_has_zero_dpa():
    ds = DamageState([], dent_depth_mm=0.0)
    assert ds.projected_damage_area_mm2 == 0.0


def test_damage_state_per_interface_area():
    e1 = DelaminationEllipse(0, (0, 0), 10, 5, 0)
    e2 = DelaminationEllipse(2, (0, 0), 8, 4, 0)
    ds = DamageState([e1, e2], dent_depth_mm=0.3)
    d = ds.per_interface_area
    assert math.isclose(d[0], math.pi * 10 * 5)
    assert math.isclose(d[2], math.pi * 8 * 4)
    assert 1 not in d
