"""Delaminations stay inside the panel (support window) that was impacted.

The union of differently oriented template ellipses spreads well beyond the
circle of the same area. Unconfined, the 40-ply hard Lovejoy & Scotti coupon
got a 108 mm wide footprint on a 102 mm coupon, and fe3d saw a damage band
through the whole section (-47% on that coupon).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from shapely.ops import unary_union

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, _ellipse_polygon
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage
from bvidfe.impact.shape_templates import distribute_damage

_HARD = [45, 0, 90, 0, -45, 0, 45, 0, -45, 0] * 2
_HARD = _HARD + _HARD[::-1]


def _bounds(ellipses):
    return unary_union([_ellipse_polygon(e) for e in ellipses]).bounds


def test_every_ellipse_fits_the_window_and_the_area_is_kept():
    hx, hy = 62.5, 37.5
    ells = distribute_damage(_HARD, 6589.0, 0.3, 0.0, (76.2, 50.8), max_half_extent_mm=(hx, hy))
    x0, y0, x1, y1 = _bounds(ells)
    assert x1 - x0 <= 2 * hx + 0.5 and y1 - y0 <= 2 * hy + 0.5
    assert DamageState(ells).projected_damage_area_mm2 == pytest.approx(6589.0, rel=0.01)


def test_unconfined_hard_footprint_is_wider_than_the_coupon():
    """The reason for the confinement."""
    x0, y0, x1, y1 = _bounds(distribute_damage(_HARD, 6589.0, 0.3, 0.0, (76.2, 50.8)))
    assert y1 - y0 > 101.6


def test_a_window_too_small_for_the_area_warns_and_saturates():
    with pytest.warns(UserWarning, match="impact window"):
        ells = distribute_damage(_HARD, 6589.0, 0.3, 0.0, (20.0, 20.0), max_half_extent_mm=(20, 20))
    x0, y0, x1, y1 = _bounds(ells)
    assert x1 - x0 <= 40.5 and y1 - y0 <= 40.5
    assert DamageState(ells).projected_damage_area_mm2 < 6589.0


def test_no_window_is_the_unconfined_template():
    a = distribute_damage(_HARD, 3000.0, 0.3, 0.0, (0.0, 0.0))
    b = distribute_damage(_HARD, 3000.0, 0.3, 0.0, (0.0, 0.0), max_half_extent_mm=(1e6, 1e6))
    # Flattened: pytest >= 9.2 rejects nested sequences in approx().
    assert [v for e in a for v in (e.major_mm, e.minor_mm)] == pytest.approx(
        [v for e in b for v in (e.major_mm, e.minor_mm)]
    )


def test_impact_to_damage_keeps_the_damage_inside_the_impacted_panel():
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], _HARD, 0.19)
    window = PanelGeometry(125.0, 75.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = impact_to_damage(ImpactEvent(67.8, ImpactorGeometry(15.875), 5.5), lam, window)
    x0, y0, x1, y1 = _bounds(ds.delaminations)
    assert x0 >= -0.5 and y0 >= -0.5 and x1 <= 125.5 and y1 <= 75.5
    assert ds.projected_damage_area_mm2 > 5000.0


def test_off_centre_impact_is_confined_by_the_nearer_edges():
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.19)
    panel = PanelGeometry(150.0, 100.0)
    ev = ImpactEvent(20.0, ImpactorGeometry(16.0), 5.5, location_xy_mm=(30.0, 25.0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = impact_to_damage(ev, lam, panel)
    x0, y0, x1, y1 = _bounds(ds.delaminations)
    assert x0 >= -0.5 and y0 >= -0.5
    assert np.isclose(np.mean([e.centroid_mm[0] for e in ds.delaminations]), 30.0)
