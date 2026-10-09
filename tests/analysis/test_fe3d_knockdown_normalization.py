"""fe3d knockdown is fe3d's damaged strength over its own undamaged strength.

fe3d's failure stresses (first-ply failure on the 3D mesh, delamination
growth) are a different quantity from the shared
``pristine_strength_MPa`` (then a thickness-weighted ply-strength average), so
dividing one by the other reported a knockdown well below 1 for an
undamaged panel (0.37 in tension, 0.03 in compression for an 8-ply QI
150x100 mm panel). The knockdown is now normalised by an undamaged fe3d run,
and applied to the shared pristine strength so that, as for every tier,
``knockdown == residual_strength_MPa / pristine_strength_MPa``.
"""

from __future__ import annotations

import math

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, DamageState, DelaminationEllipse, MeshParams
from bvidfe.analysis.fe_tier import _fe3d_cai_first_ply_failure, fe3d_tai
from bvidfe.analysis.semi_analytical import weakest_sublaminate_growth
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY

_LAYUP = [0, 45, -45, 90, 90, -45, 45, 0]
_UNDAMAGED = DamageState(delaminations=[], dent_depth_mm=0.0)
_DELAMINATED = DamageState(
    delaminations=[DelaminationEllipse(3, (20.0, 15.0), 16.0, 12.0, 0.0)],
    dent_depth_mm=0.0,
)


def _config(loading: str, damage: DamageState) -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_LAYUP,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(40.0, 30.0),
        loading=loading,
        tier="fe3d",
        damage=damage,
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=4.0),
    )


def _raw_fe3d_strength(cfg: AnalysisConfig, damage: DamageState) -> float:
    """Uncapped fe3d strength in MPa, computed from the fe_tier channels."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], _LAYUP, 0.152)
    if cfg.loading == "compression":
        growth = weakest_sublaminate_growth(lam, damage, boundary=cfg.panel.boundary)
        sigma_growth = growth.growth_stress_MPa if growth is not None else math.inf
        return min(sigma_growth, _fe3d_cai_first_ply_failure(cfg, damage, lam, math.inf))
    return fe3d_tai(cfg, damage, lam, math.inf)


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_undamaged_fe3d_knockdown_is_one(loading):
    r = BvidAnalysis(_config(loading, _UNDAMAGED)).run()
    assert r.knockdown == pytest.approx(1.0, rel=1e-12)
    assert r.residual_strength_MPa == pytest.approx(r.pristine_strength_MPa, rel=1e-12)


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_fe3d_knockdown_is_ratio_of_damaged_to_undamaged(loading):
    cfg = _config(loading, _DELAMINATED)
    r = BvidAnalysis(cfg).run()
    expected = min(1.0, _raw_fe3d_strength(cfg, _DELAMINATED) / _raw_fe3d_strength(cfg, _UNDAMAGED))
    assert r.knockdown == pytest.approx(expected, rel=1e-9)
    assert r.residual_strength_MPa == pytest.approx(r.knockdown * r.pristine_strength_MPa)
    assert any("undamaged" in n for n in r.notes), r.notes


def test_delamination_lowers_fe3d_compression_knockdown():
    r = BvidAnalysis(_config("compression", _DELAMINATED)).run()
    assert 0.0 < r.knockdown < 1.0
