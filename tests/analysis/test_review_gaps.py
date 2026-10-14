"""Gaps found while bringing the documentation up to date."""

from __future__ import annotations

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis
from bvidfe.analysis import bvid
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.analysis.semi_analytical import semi_analytical_cai
from bvidfe.damage.state import DamageState, DelaminationEllipse

_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]


def _cfg(damage: DamageState, tier: str = "empirical") -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        loading="compression",
        tier=tier,
        damage=damage,
    )


@pytest.mark.parametrize("index", [7, 12])
def test_delamination_interface_outside_the_layup_is_rejected(index):
    """An 8-ply laminate has interfaces 0-6. Index 7+ used to be silently
    ignored by the closed-form tiers and raised IndexError in fe3d."""
    damage = DamageState([DelaminationEllipse(index, (75.0, 50.0), 20.0, 10.0, 0.0)], 0.2)
    with pytest.raises(ValueError, match="interface_index"):
        _cfg(damage)


def test_last_interface_is_accepted():
    damage = DamageState([DelaminationEllipse(6, (75.0, 50.0), 20.0, 10.0, 0.0)], 0.2)
    assert BvidAnalysis(_cfg(damage)).run().knockdown < 1.0


def test_fe3d_buckling_notes_reach_the_results(monkeypatch):
    note = "fe3d buckling: degenerate"
    monkeypatch.setattr(bvid.BvidAnalysis, "_fe3d_raw_strength", lambda self, lam, dmg: 100.0)
    monkeypatch.setattr(bvid, "fe3d_cai_buckling", lambda cfg, dmg, lam, s: (s, 0.0, [note]))
    damage = DamageState([DelaminationEllipse(3, (75.0, 50.0), 20.0, 10.0, 0.0)], 0.2)
    cfg = _cfg(damage, tier="fe3d")
    r = BvidAnalysis(cfg).run()
    assert note in r.notes
    assert r.buckling_eigenvalues is None


def test_semi_analytical_residual_is_a_python_float():
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], _QI8, 0.152)
    damage = DamageState([DelaminationEllipse(3, (75.0, 50.0), 25.0, 15.0, 0.0)], 0.2)
    r = semi_analytical_cai(lam, damage, 600.0, 15000.0)
    assert type(r.residual_strength_MPa) is float
