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


@pytest.mark.parametrize("theta", [0.0, 30.0, 45.0, 90.0, -45.0, 135.0, 180.0])
def test_half_extents_match_the_footprint_bounding_box(theta):
    """``half_extents_mm`` is the half-width of the rotated ellipse's
    axis-aligned bounding box, which is what the footprint polygon spans
    (#165: the same extents feed the TAI width, the peanut template and the
    sublaminate buckling rectangle)."""
    from bvidfe.damage.state import _ellipse_polygon

    e = DelaminationEllipse(0, (75.0, 50.0), 20.0, 8.0, theta)
    ex, ey = e.half_extents_mm()
    assert isinstance(ex, float) and isinstance(ey, float)
    minx, miny, maxx, maxy = _ellipse_polygon(e, n_pts=720).bounds
    # The 720-gon is inscribed, so its box is at most ~1e-5 relative inside.
    assert 0.5 * (maxx - minx) == pytest.approx(ex, rel=2e-5)
    assert 0.5 * (maxy - miny) == pytest.approx(ey, rel=2e-5)
    assert ex == pytest.approx(
        math.hypot(20.0 * math.cos(math.radians(theta)), 8.0 * math.sin(math.radians(theta)))
    )
    assert ey == pytest.approx(
        math.hypot(20.0 * math.sin(math.radians(theta)), 8.0 * math.cos(math.radians(theta)))
    )


def test_half_extents_swap_at_90_degrees_and_are_bounded_by_the_axes():
    assert DelaminationEllipse(0, (0.0, 0.0), 20.0, 8.0, 0.0).half_extents_mm() == (20.0, 8.0)
    ex, ey = DelaminationEllipse(0, (0.0, 0.0), 20.0, 8.0, 90.0).half_extents_mm()
    assert (ex, ey) == pytest.approx((8.0, 20.0))
    ex, ey = DelaminationEllipse(0, (0.0, 0.0), 20.0, 8.0, 45.0).half_extents_mm()
    assert ex == pytest.approx(ey)
    assert 8.0 < ex < 20.0


def test_ellipse_half_extents_vectorised():
    import numpy as np

    from bvidfe.damage.state import ellipse_half_extents_mm

    ex, ey = ellipse_half_extents_mm(
        np.array([20.0, 20.0]), np.array([8.0, 8.0]), np.array([0.0, 90.0])
    )
    assert ex == pytest.approx([20.0, 8.0])
    assert ey == pytest.approx([8.0, 20.0])
