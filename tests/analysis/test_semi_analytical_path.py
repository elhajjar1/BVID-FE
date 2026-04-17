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
        tier="semi_analytical",
    )
    kw.update(overrides)
    return AnalysisConfig(**kw)


def test_semi_analytical_cai_runs_end_to_end():
    cfg = _base_cfg(impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5))
    r = BvidAnalysis(cfg).run()
    assert r.tier_used == "semi_analytical"
    assert 0 < r.knockdown < 1.0
    assert r.buckling_eigenvalues is not None
    assert r.critical_sublaminate is not None


def test_semi_analytical_tai_runs():
    ds = DamageState(
        [DelaminationEllipse(3, (75, 50), 20, 12, 45)],
        dent_depth_mm=0.4,
    )
    cfg = _base_cfg(loading="tension", damage=ds)
    r = BvidAnalysis(cfg).run()
    assert r.tier_used == "semi_analytical"
    assert r.knockdown < 1.0
    assert r.residual_strength_MPa > 0


def test_semi_analytical_vs_empirical_cai_both_positive():
    impact = ImpactEvent(25.0, ImpactorGeometry(), mass_kg=5.5)
    cfg_e = _base_cfg(tier="empirical", impact=impact)
    cfg_sa = _base_cfg(tier="semi_analytical", impact=impact)
    r_e = BvidAnalysis(cfg_e).run()
    r_sa = BvidAnalysis(cfg_sa).run()
    assert r_e.residual_strength_MPa > 0
    assert r_sa.residual_strength_MPa > 0
