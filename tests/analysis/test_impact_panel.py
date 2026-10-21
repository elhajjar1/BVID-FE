"""AnalysisConfig.impact_panel: impact on one panel, test on another.

Hasebe et al. (2025) impact 80 x 80 mm plates on a 60 x 60 mm clamped
window, then trim them to an 80 x 50 mm compression coupon.
"""

from __future__ import annotations

import warnings

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage

_LAYUP = [45, 0, -45, 90] * 2 + [90, -45, 0, 45] * 2
_COUPON = PanelGeometry(80.0, 50.0, "simply_supported")
_WINDOW = PanelGeometry(60.0, 60.0, "clamped")
_EVENT = ImpactEvent(energy_J=10.0, impactor=ImpactorGeometry(12.7), mass_kg=5.5)


def _cfg(**kw):
    base = dict(
        material="T800S/3900-2B",
        layup_deg=_LAYUP,
        ply_thickness_mm=0.19,
        panel=_COUPON,
        impact=_EVENT,
    )
    return AnalysisConfig(**{**base, **kw})


def test_damage_is_mapped_on_the_impact_panel():
    lam = Laminate(MATERIAL_LIBRARY["T800S/3900-2B"], _LAYUP, 0.19)
    expected = impact_to_damage(_EVENT, lam, _WINDOW)
    r = BvidAnalysis(_cfg(impact_panel=_WINDOW)).run()
    assert r.dpa_mm2 == pytest.approx(expected.projected_damage_area_mm2, rel=1e-6)
    assert r.damage.dent_depth_mm == pytest.approx(expected.dent_depth_mm)


def test_damage_is_centred_on_the_test_panel():
    r = BvidAnalysis(_cfg(impact_panel=_WINDOW)).run()
    assert r.damage.delaminations
    for e in r.damage.delaminations:
        assert e.centroid_mm == pytest.approx((40.0, 25.0))


def test_impact_panel_changes_the_damage_cap():
    """The 80% DPA cap is set by the panel that was impacted."""
    big = ImpactEvent(energy_J=200.0, impactor=ImpactorGeometry(12.7), mass_kg=5.5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        on_window = BvidAnalysis(_cfg(impact=big, impact_panel=PanelGeometry(30.0, 30.0))).run()
        on_coupon = BvidAnalysis(_cfg(impact=big)).run()
    assert on_window.dpa_mm2 == pytest.approx(0.8 * 30.0 * 30.0, rel=0.02)
    assert on_coupon.dpa_mm2 > on_window.dpa_mm2


def test_impact_panel_defaults_to_the_test_panel():
    a = BvidAnalysis(_cfg()).run()
    b = BvidAnalysis(_cfg(impact_panel=_COUPON)).run()
    assert a.residual_strength_MPa == pytest.approx(b.residual_strength_MPa)
    assert a.dpa_mm2 == pytest.approx(b.dpa_mm2)


def test_impact_panel_needs_an_impact():
    with pytest.raises(ValueError, match="impact_panel"):
        _cfg(impact=None, damage=DamageState([], dent_depth_mm=0.0), impact_panel=_WINDOW)
