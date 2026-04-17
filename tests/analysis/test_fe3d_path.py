from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import ImpactEvent


def _small_cfg(**overrides):
    """Small mesh to keep FE tests fast."""
    kw = dict(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(10, 5),
        loading="compression",
        tier="fe3d",
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.5),
    )
    kw.update(overrides)
    return AnalysisConfig(**kw)


def test_fe3d_cai_end_to_end():
    ds = DamageState([DelaminationEllipse(1, (5, 2.5), 3, 1.5, 0)], dent_depth_mm=0.2)
    cfg = _small_cfg(damage=ds)
    r = BvidAnalysis(cfg).run()
    assert r.tier_used == "fe3d"
    assert 0 < r.residual_strength_MPa
    assert 0 < r.knockdown <= 1.0


def test_fe3d_tai_end_to_end():
    ds = DamageState([DelaminationEllipse(1, (5, 2.5), 3, 1.5, 0)], dent_depth_mm=0.2)
    cfg = _small_cfg(damage=ds, loading="tension")
    r = BvidAnalysis(cfg).run()
    assert r.tier_used == "fe3d"
    assert r.residual_strength_MPa > 0


def test_fe3d_impact_path_runs():
    cfg = _small_cfg(impact=ImpactEvent(3.0, ImpactorGeometry(), mass_kg=5.5))
    r = BvidAnalysis(cfg).run()
    assert r.tier_used == "fe3d"
    assert r.residual_strength_MPa > 0
