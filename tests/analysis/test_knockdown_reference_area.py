"""The CAI knockdown scales damage by a fixed area, not the panel's.

Dividing by the panel area sent the knockdown to 1 as the panel grew,
although impact damage is a local notch; it also made Hasebe et al.'s
80 x 50 mm coupons look far more damage-sensitive than ASTM ones.
"""

from __future__ import annotations

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.soutis_openhole import SOUTIS_REFERENCE_AREA_MM2

_LAYUP = [45, 0, -45, 90] * 2 + [90, -45, 0, 45] * 2


def _knockdown(Lx, Ly, tier="empirical"):
    damage = DamageState(
        [DelaminationEllipse(7, (Lx / 2, Ly / 2), 20.0, 15.0, 0.0)], dent_depth_mm=0.2
    )
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=_LAYUP,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(Lx, Ly),
        tier=tier,
        damage=damage,
    )
    return BvidAnalysis(cfg).run()


def test_reference_area_is_the_astm_d7137_coupon():
    assert SOUTIS_REFERENCE_AREA_MM2 == 150.0 * 100.0


@pytest.mark.parametrize("tier", ["empirical", "semi_analytical"])
def test_same_damage_gives_the_same_soutis_knockdown_on_any_panel(tier):
    small = _knockdown(80.0, 50.0, tier)
    astm = _knockdown(150.0, 100.0, tier)
    large = _knockdown(600.0, 600.0, tier)
    assert small.dpa_mm2 == pytest.approx(astm.dpa_mm2, rel=1e-6)
    if tier == "empirical":
        assert small.knockdown == pytest.approx(astm.knockdown, rel=1e-12)
        assert large.knockdown == pytest.approx(astm.knockdown, rel=1e-12)
    # Panel size still enters the semi-analytical tier through the boundary
    # stiffness of the growth check, never by relaxing the Soutis bound.
    assert large.knockdown <= astm.knockdown + 1e-12


def test_knockdown_formula():
    m = MATERIAL_LIBRARY["IM7/8552"]
    r = _knockdown(150.0, 100.0)
    expected = 1.0 / (1.0 + m.soutis_k_s * (r.dpa_mm2 / SOUTIS_REFERENCE_AREA_MM2) ** m.soutis_m)
    assert r.knockdown == pytest.approx(expected, rel=1e-12)
