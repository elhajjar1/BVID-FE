import pytest

from bvidfe.analysis.config import AnalysisConfig, MeshParams
from bvidfe.analysis.fe_tier import fe3d_cai, fe3d_tai
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import ImpactEvent


@pytest.fixture
def small_cfg():
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(10.0, 5.0),
        loading="compression",
        tier="fe3d",
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.5),
    )


def test_fe3d_cai_pristine_returns_positive_strength(small_cfg):
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], small_cfg.layup_deg, small_cfg.ply_thickness_mm)
    damage = DamageState([], dent_depth_mm=0.0)
    sigma = fe3d_cai(small_cfg, damage, lam, sigma_pristine_MPa=500.0)
    assert sigma > 0


def test_fe3d_cai_damaged_less_than_pristine(small_cfg):
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], small_cfg.layup_deg, small_cfg.ply_thickness_mm)
    # Damage with a modest ellipse at interface 1
    ds = DamageState([DelaminationEllipse(1, (5, 2.5), 3, 1.5, 0)], dent_depth_mm=0.2)
    sigma_pristine = fe3d_cai(small_cfg, DamageState([], 0.0), lam, sigma_pristine_MPa=500.0)
    sigma_damaged = fe3d_cai(small_cfg, ds, lam, sigma_pristine_MPa=500.0)
    assert 0 < sigma_damaged <= sigma_pristine


def test_fe3d_tai_pristine_positive(small_cfg):
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], small_cfg.layup_deg, small_cfg.ply_thickness_mm)
    damage = DamageState([], 0.0)
    sigma = fe3d_tai(small_cfg, damage, lam, sigma_pristine_MPa=800.0)
    assert sigma > 0
