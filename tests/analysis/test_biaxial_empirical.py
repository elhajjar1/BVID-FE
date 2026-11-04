"""Biaxial loading in the empirical tier (issue #115, PR 1), and the
``NotImplementedError`` of the tiers that do not model it yet."""

from __future__ import annotations

import math
import time

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, MeshParams, PanelGeometry
from bvidfe.analysis import fe_tier
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.laminate_strength import unnotched_strength
from bvidfe.failure.lekhnitskii import biaxial_point_stress_factor
from bvidfe.failure.soutis_openhole import fibre_damage_tai, tai_hole_radius_mm

QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
CROSS = [0, 90, 0, 90, 90, 0, 90, 0]
_IM7 = MATERIAL_LIBRARY["IM7/8552"]
_NEW_FIELDS = (
    "pristine_transverse_stress_MPa",
    "residual_transverse_stress_MPa",
    "biaxial_interaction_factor",
)


def _cfg(loading="compression", beta=0.0, energy=15.0, layup=QI8, **kw):
    args = dict(
        material="IM7/8552",
        layup_deg=list(layup),
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        loading=loading,
        impact=ImpactEvent(energy_J=energy),
        transverse_stress_ratio=beta,
    )
    args.update(kw)
    return AnalysisConfig(**args)


def _run(cfg):
    return BvidAnalysis(cfg).run()


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_uniaxial_run_has_no_biaxial_fields_or_notes(loading):
    r = _run(_cfg(loading))
    for name in _NEW_FIELDS:
        assert getattr(r, name) is None
        assert r.to_dict()[name] is None
    assert not any("Biaxial" in n for n in r.notes)
    assert "biaxial" not in r.summary()


@pytest.mark.parametrize("layup", [QI8, CROSS])
@pytest.mark.parametrize("beta", [-1.0, -0.5, 0.5, 1.0])
@pytest.mark.parametrize("energy", [5.0, 30.0])
def test_cai_is_homothetic(layup, beta, energy):
    """D6: the Soutis knockdown does not depend on the ratio; the residual
    is that knockdown times the pristine strength along the ray."""
    uni = _run(_cfg("compression", 0.0, energy, layup))
    bi = _run(_cfg("compression", beta, energy, layup))
    lam = Laminate(_IM7, layup, 0.152)
    assert bi.pristine_strength_MPa == unnotched_strength(lam, "compression", beta)
    assert bi.knockdown == pytest.approx(uni.knockdown, rel=1e-12)
    assert bi.residual_strength_MPa == pytest.approx(
        uni.knockdown * bi.pristine_strength_MPa, rel=1e-12
    )
    assert bi.biaxial_interaction_factor == pytest.approx(
        bi.pristine_strength_MPa / uni.pristine_strength_MPa, rel=1e-12
    )


@pytest.mark.parametrize("layup", [QI8, CROSS])
@pytest.mark.parametrize("beta", [-1.0, -0.5, 0.5, 1.0])
@pytest.mark.parametrize("energy", [5.0, 30.0])
def test_tai_knockdown_carries_the_lekhnitskii_factor(layup, beta, energy):
    """D7: kd(beta) = kd(0) * F(beta) at the run's equivalent hole."""
    uni = _run(_cfg("tension", 0.0, energy, layup))
    bi = _run(_cfg("tension", beta, energy, layup))
    lam = Laminate(_IM7, layup, 0.152)
    R = tai_hole_radius_mm(_IM7, bi.damage, 100.0)
    F = biaxial_point_stress_factor(lam, R, _IM7.wn_d0_mm, beta)
    assert F != 1.0
    assert bi.pristine_strength_MPa == unnotched_strength(lam, "tension", beta)
    assert bi.knockdown == pytest.approx(uni.knockdown * F, rel=1e-12)
    assert bi.biaxial_interaction_factor == pytest.approx(
        bi.residual_strength_MPa / uni.residual_strength_MPa, rel=1e-12
    )
    assert any(f"F = {F:.3f}" in n for n in bi.notes if n.startswith("TAI:"))


def test_tai_factor_is_small_at_bvid_hole_sizes():
    """The TAI hole is 1-4 mm across, so the lateral stress moves the TAI
    knockdown by a few percent only (probe 2.4)."""
    for beta in (-1.0, 1.0):
        uni = _run(_cfg("tension", 0.0, 15.0))
        bi = _run(_cfg("tension", beta, 15.0))
        assert bi.knockdown / uni.knockdown == pytest.approx(1.0, abs=0.05)


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", [-0.6, 0.8])
def test_transverse_stresses_are_signed_sigma_y_at_failure(loading, beta):
    r = _run(_cfg(loading, beta))
    s = 1.0 if loading == "tension" else -1.0
    assert r.pristine_transverse_stress_MPa == pytest.approx(s * beta * r.pristine_strength_MPa)
    assert r.residual_transverse_stress_MPa == pytest.approx(s * beta * r.residual_strength_MPa)
    d = r.to_dict()
    for name in _NEW_FIELDS:
        assert d[name] == getattr(r, name)
    text = r.summary()
    assert "biaxial (sigma_y/sigma_x" in text and "biaxial_interaction_factor" in text


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_interaction_factor_is_against_the_uniaxial_run(loading):
    uni = _run(_cfg(loading, 0.0))
    for beta in (-1.0, 0.3):
        bi = _run(_cfg(loading, beta))
        assert bi.biaxial_interaction_factor == pytest.approx(
            bi.residual_strength_MPa / uni.residual_strength_MPa, rel=1e-12
        )


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("energy", [5.0, 15.0, 30.0])
def test_a_tiny_ratio_is_continuous_with_the_uniaxial_run(loading, energy):
    uni = _run(_cfg(loading, 0.0, energy))
    for beta in (1e-9, -1e-9):
        bi = _run(_cfg(loading, beta, energy))
        for name in ("residual_strength_MPa", "pristine_strength_MPa", "knockdown", "dpa_mm2"):
            assert getattr(bi, name) == pytest.approx(getattr(uni, name), rel=1e-6)
        assert bi.biaxial_interaction_factor == pytest.approx(1.0, rel=1e-6)
        assert abs(bi.pristine_transverse_stress_MPa) < 1e-5


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", [-1.0, -0.5, 0.5, 1.0])
def test_every_biaxial_run_says_it_is_experimental(loading, beta):
    r = _run(_cfg(loading, beta))
    note = r.notes[0]
    assert note.startswith(f"Biaxial sigma_y/sigma_x = {beta:+.2f}: experimental, not validated")
    assert "no biaxial residual-strength data" in note
    if loading == "compression":
        assert "homothetic" in note
    else:
        assert "characteristic point" in note
    if loading == "tension" and beta > 0:
        assert "every ply has cracked" in note  # the D5 caveat
    else:
        assert "every ply has cracked" not in note
    assert sum("experimental" in n for n in r.notes) == 1


def test_undamaged_panel_keeps_the_pristine_strength():
    undamaged = DamageState([], dent_depth_mm=0.0)
    for loading in ("compression", "tension"):
        for beta in (-1.0, 0.5):
            r = _run(_cfg(loading, beta, impact=None, damage=undamaged))
            assert r.knockdown == 1.0
            assert r.residual_strength_MPa == r.pristine_strength_MPa


def test_inspection_driven_biaxial_tai():
    damage = DamageState([DelaminationEllipse(2, (75.0, 50.0), 20.0, 12.0, 0.0)])
    uni = _run(_cfg("tension", 0.0, impact=None, damage=damage))
    bi = _run(_cfg("tension", 0.5, impact=None, damage=damage))
    lam = Laminate(_IM7, QI8, 0.152)
    R = tai_hole_radius_mm(_IM7, damage, 100.0)
    F = biaxial_point_stress_factor(lam, R, _IM7.wn_d0_mm, 0.5)
    assert bi.knockdown == pytest.approx(uni.knockdown * F, rel=1e-12)


def test_fibre_damage_tai_needs_the_laminate_for_a_ratio():
    damage = DamageState([DelaminationEllipse(2, (75.0, 50.0), 20.0, 12.0, 0.0)])
    assert fibre_damage_tai(_IM7, damage, 1.0, 3.0, 100.0) == fibre_damage_tai(
        _IM7, damage, 1.0, 3.0, 100.0, transverse_stress_ratio=0.0
    )
    with pytest.raises(ValueError, match="laminate"):
        fibre_damage_tai(_IM7, damage, 1.0, 3.0, 100.0, transverse_stress_ratio=0.5)


@pytest.mark.parametrize("tier", ["semi_analytical", "fe3d"])
@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_tiers_without_a_biaxial_model_refuse_a_ratio(tier, loading, monkeypatch):
    """No silent uniaxial fallback; the message names the follow-up PR. It
    is raised before any fe3d solve."""

    def no_solve(*a, **k):
        raise AssertionError("fe3d solved")

    monkeypatch.setattr(fe_tier, "solve_linear_static", no_solve)
    cfg = _cfg(loading, 0.5, tier=tier, mesh=MeshParams(in_plane_size_mm=10.0))
    t0 = time.perf_counter()
    with pytest.raises(NotImplementedError) as exc:
        _run(cfg)
    assert time.perf_counter() - t0 < 5.0
    msg = str(exc.value)
    assert f"tier={tier!r}" in msg and "#115" in msg
    assert ("PR 2" if tier == "semi_analytical" else "PR 3") in msg
    assert "tier='empirical'" in msg


@pytest.mark.parametrize(
    "call",
    [
        lambda c, d, lam: fe_tier.fe3d_cai_buckling(c, d, lam, math.inf),
        lambda c, d, lam: fe_tier.fe3d_cai_channels_and_modulus(c, d, lam, math.inf),
        lambda c, d, lam: fe_tier.fe3d_tai_and_modulus(c, d, lam, math.inf),
        lambda c, d, lam: fe_tier.tai_fibre_damage_hole(c, d, lam),
    ],
)
def test_fe3d_entry_points_refuse_a_ratio_too(call):
    cfg = _cfg("compression", -0.3, tier="fe3d", mesh=MeshParams(in_plane_size_mm=10.0))
    damage = DamageState([DelaminationEllipse(2, (75.0, 50.0), 20.0, 12.0, 0.0)])
    with pytest.raises(NotImplementedError, match="#115"):
        call(cfg, damage, Laminate(_IM7, QI8, 0.152))


@pytest.mark.parametrize("tier", ["semi_analytical", "fe3d"])
def test_tiers_without_a_biaxial_model_run_at_zero_ratio(tier):
    cfg = _cfg("compression", 0.0, 5.0, tier=tier, mesh=MeshParams(in_plane_size_mm=10.0))
    r = _run(cfg)
    assert 0 < r.knockdown <= 1 and r.biaxial_interaction_factor is None
