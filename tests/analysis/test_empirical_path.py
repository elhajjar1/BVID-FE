import math

import pytest

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import ImpactEvent


def _base_cfg(**overrides):
    kw = dict(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90] * 4,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier="empirical",
    )
    kw.update(overrides)
    return AnalysisConfig(**kw)


def test_config_requires_exactly_one_of_impact_or_damage():
    with pytest.raises((AssertionError, ValueError)):
        _base_cfg(impact=None, damage=None)
    ds = DamageState([DelaminationEllipse(0, (0, 0), 10, 5, 0)], 0.3)
    ev = ImpactEvent(10.0, ImpactorGeometry(), mass_kg=5.5)
    with pytest.raises((AssertionError, ValueError)):
        _base_cfg(impact=ev, damage=ds)


def test_empirical_cai_impact_path():
    cfg = _base_cfg(impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5))
    r = BvidAnalysis(cfg).run()
    assert 0.3 < r.knockdown < 1.0
    assert r.residual_strength_MPa < r.pristine_strength_MPa
    assert r.damage.projected_damage_area_mm2 > 0
    assert r.tier_used == "empirical"


def test_empirical_damage_path_equivalent_to_impact_path():
    cfg_imp = _base_cfg(impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5))
    r_imp = BvidAnalysis(cfg_imp).run()
    cfg_dmg = _base_cfg(damage=r_imp.damage)
    r_dmg = BvidAnalysis(cfg_dmg).run()
    assert abs(r_imp.knockdown - r_dmg.knockdown) < 1e-6


def test_tension_path_runs():
    ds = DamageState(
        [DelaminationEllipse(3, (75, 50), 20, 12, 45)],
        dent_depth_mm=0.4,
    )
    cfg = _base_cfg(loading="tension", damage=ds)
    r = BvidAnalysis(cfg).run()
    assert r.knockdown < 1.0
    assert r.residual_strength_MPa > 0


def test_below_threshold_returns_pristine_strength():
    # Tiny energy => impact_to_damage returns empty damage => empirical returns pristine
    cfg = _base_cfg(impact=ImpactEvent(0.01, ImpactorGeometry(), mass_kg=5.5))
    r = BvidAnalysis(cfg).run()
    assert math.isclose(r.knockdown, 1.0, rel_tol=1e-6)
    assert r.damage.projected_damage_area_mm2 == 0.0


def test_analysis_results_summary_runs():
    cfg = _base_cfg(impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5))
    r = BvidAnalysis(cfg).run()
    s = r.summary()
    assert isinstance(s, str)
    assert "knockdown" in s.lower() or "residual" in s.lower()


def test_to_dict_round_trip():
    cfg = _base_cfg(impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5))
    r = BvidAnalysis(cfg).run()
    d = r.to_dict()
    assert "residual_strength_MPa" in d
    assert "knockdown" in d
    assert "damage" in d
