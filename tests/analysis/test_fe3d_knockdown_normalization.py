"""fe3d knockdown is fe3d's damaged strength over its own undamaged strength.

fe3d's failure stresses (first-ply failure on the 3D mesh, delamination
growth) are a different quantity from the shared
``pristine_strength_MPa`` (then a thickness-weighted ply-strength average), so
dividing one by the other reported a knockdown well below 1 for an
undamaged panel (0.37 in tension, 0.03 in compression for an 8-ply QI
150x100 mm panel). The knockdown is now normalised by an undamaged fe3d run,
and applied to the shared pristine strength so that, as for every tier,
``knockdown == residual_strength_MPa / pristine_strength_MPa``.

Only fe3d's own first-ply-failure stress is normalised. In compression the
delamination growth stress is the closed-form far-field stress that
``semi_analytical`` uses, and it used to be rescaled by
``sigma0 / FPF_undamaged`` along with first-ply failure, which put fe3d above
``semi_analytical`` wherever growth governs (issue #166). It is now compared
with the normalised first-ply failure as it is.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, DamageState, DelaminationEllipse, MeshParams
from bvidfe.analysis import bvid
from bvidfe.analysis.fe_tier import _fe3d_cai_first_ply_failure, fe3d_tai
from bvidfe.analysis.semi_analytical import weakest_sublaminate_growth
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.mapping import ImpactEvent

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


def _expected_residual(cfg: AnalysisConfig, damage: DamageState, sigma0: float) -> float:
    """fe3d residual strength in MPa, computed from the fe_tier channels:
    first-ply failure normalised by the undamaged run and, in compression,
    the delamination growth stress as it is."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], _LAYUP, 0.152)
    if cfg.loading == "compression":
        fpf = _fe3d_cai_first_ply_failure
        ratio = fpf(cfg, damage, lam, math.inf) / fpf(cfg, _UNDAMAGED, lam, math.inf)
        growth = weakest_sublaminate_growth(lam, damage, boundary=cfg.panel.boundary)
        sigma_growth = growth.growth_stress_MPa if growth is not None else math.inf
        return min(sigma_growth, sigma0 * min(1.0, ratio))
    ratio = fe3d_tai(cfg, damage, lam, math.inf) / fe3d_tai(cfg, _UNDAMAGED, lam, math.inf)
    return sigma0 * min(1.0, ratio)


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_undamaged_fe3d_knockdown_is_one(loading):
    r = BvidAnalysis(_config(loading, _UNDAMAGED)).run()
    assert r.knockdown == pytest.approx(1.0, rel=1e-12)
    assert r.residual_strength_MPa == pytest.approx(r.pristine_strength_MPa, rel=1e-12)


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_fe3d_knockdown_is_ratio_of_damaged_to_undamaged(loading):
    cfg = _config(loading, _DELAMINATED)
    r = BvidAnalysis(cfg).run()
    expected = _expected_residual(cfg, _DELAMINATED, r.pristine_strength_MPa)
    assert r.residual_strength_MPa == pytest.approx(expected, rel=1e-9)
    assert r.residual_strength_MPa == pytest.approx(r.knockdown * r.pristine_strength_MPa)
    assert any("undamaged" in n for n in r.notes), r.notes


def test_delamination_lowers_fe3d_compression_knockdown():
    r = BvidAnalysis(_config("compression", _DELAMINATED)).run()
    assert 0.0 < r.knockdown < 1.0


def _growth_at(stress_MPa: float):
    """A real SublaminateGrowth record with its growth stress replaced."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], _LAYUP, 0.152)
    growth = weakest_sublaminate_growth(lam, _DELAMINATED)
    assert growth is not None
    return dataclasses.replace(growth, growth_stress_MPa=stress_MPa)


@pytest.mark.parametrize(
    "growth_MPa, fpf_damaged, governs",
    [(100.0, 300.0, "delamination growth"), (400.0, 150.0, "first-ply failure")],
)
def test_fe3d_cai_compares_growth_and_normalised_fpf_on_one_scale(
    monkeypatch, growth_MPa, fpf_damaged, governs
):
    """With fe3d's undamaged first-ply failure at 400 MPa (below sigma0 =
    601.3 MPa), a 100 MPa growth stress used to come out at
    601.3 x 100 / 400 = 150 MPa. It is now 100 MPa, the semi_analytical
    value; first-ply failure is still normalised (601.3 x 150 / 400)."""

    def channels(cfg, damage, lam, cap):
        if damage.delaminations:
            return _growth_at(growth_MPa), fpf_damaged, 1.0
        return None, 400.0, 1.0

    monkeypatch.setattr(bvid, "fe3d_cai_channels_and_modulus", channels)
    r = BvidAnalysis(_config("compression", _DELAMINATED)).run()
    sigma0 = r.pristine_strength_MPa
    assert r.residual_strength_MPa == pytest.approx(min(growth_MPa, sigma0 * fpf_damaged / 400.0))
    assert r.knockdown == pytest.approx(r.residual_strength_MPa / sigma0)
    assert any(n.startswith(f"fe3d CAI: {governs} governs") for n in r.notes), r.notes


def test_fe3d_cai_growth_stress_is_capped_at_pristine(monkeypatch):
    monkeypatch.setattr(
        bvid,
        "fe3d_cai_channels_and_modulus",
        lambda cfg, damage, lam, cap: (
            (_growth_at(1e4), 1e4, 1.0) if damage.delaminations else (None, 1e3, 1.0)
        ),
    )
    r = BvidAnalysis(_config("compression", _DELAMINATED)).run()
    assert r.knockdown == 1.0


def _impact_config(material: str, layup: list, tier: str) -> AnalysisConfig:
    return AnalysisConfig(
        material=material,
        layup_deg=layup,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        loading="compression",
        tier=tier,
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0) if tier == "fe3d" else None,
    )


# Delamination growth governs both tiers at 5 J on a 150 x 100 mm panel
# (10 mm fe3d mesh), and fe3d's undamaged LaRC05 first-ply failure is below
# sigma0 (FPF_undamaged / sigma0 = 0.64 and 0.86), so rescaling the growth
# stress put fe3d 27% and 3% above semi_analytical (282.9 against 221.9 MPa
# and 298.1 against 288.7 MPa). [+-45]2s is not here: its damaged
# first-ply failure governs (129 MPa against a 217 MPa growth stress).
_GROWTH_GOVERNS = [
    pytest.param("IM7/8552", [45, -45, 90, 90, -45, 45], id="IM7/8552 [45/-45/90]s"),
    pytest.param("AS4/3501-6", [0, 45, -45, 90, 90, -45, 45, 0], id="AS4/3501-6 QI8"),
]


@pytest.mark.parametrize("material, layup", _GROWTH_GOVERNS)
def test_fe3d_cai_is_semi_analytical_where_growth_governs(material, layup):
    """Where growth governs, fe3d is the growth stress itself, as in
    semi_analytical, not the growth stress rescaled by sigma0 / FPF_undamaged."""
    fe3d = BvidAnalysis(_impact_config(material, layup, "fe3d")).run()
    semi = BvidAnalysis(_impact_config(material, layup, "semi_analytical")).run()
    assert any(
        n.startswith("fe3d CAI: delamination growth governs") for n in fe3d.notes
    ), fe3d.notes
    assert fe3d.residual_strength_MPa <= semi.residual_strength_MPa
    # Both are the same closed-form growth stress.
    assert fe3d.residual_strength_MPa == pytest.approx(semi.residual_strength_MPa, rel=1e-12)
