"""Tests for the fe3d CAI buckling channel (Rayleigh-Ritz closed-form
delegation, #129)."""

import contextlib

import pytest

from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.analysis.fe_tier import fe3d_cai_buckling
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import DPACapClipWarning, ImpactEvent


@contextlib.contextmanager
def _tiny_panel_impact_warnings():
    """An impact on ``small_cfg``'s 10 x 5 mm panel predicts far more DPA
    than the panel holds: it is clipped to 80% of the panel area, and that
    area does not fit inside the impact window."""
    with (
        pytest.warns(DPACapClipWarning),
        pytest.warns(UserWarning, match="does not fit inside the impact window"),
    ):
        yield


@pytest.fixture
def small_cfg():
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(10, 5),
        loading="compression",
        tier="fe3d",
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.5),
    )


def test_fe3d_cai_buckling_returns_positive_strength(small_cfg):
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], small_cfg.layup_deg, small_cfg.ply_thickness_mm)
    damage = DamageState([], dent_depth_mm=0.0)
    sigma, lambda_crit, notes = fe3d_cai_buckling(small_cfg, damage, lam, sigma_pristine_MPa=500.0)
    assert sigma > 0
    assert lambda_crit > 0
    # A clean solve notes the whole-panel buckling stress, which the
    # knockdown leaves out, and that the onset is the closed form (#133).
    assert len(notes) == 2, notes
    assert notes[0].startswith("fe3d whole-panel buckling at"), notes
    assert notes[1].startswith(f"fe3d buckling onset {lambda_crit:.1f} MPa (whole panel;"), notes
    assert "closed-form rayleigh-ritz" in notes[1].lower(), notes


def test_fe3d_cai_buckling_damaged_less_than_pristine(small_cfg):
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], small_cfg.layup_deg, small_cfg.ply_thickness_mm)
    ds = DamageState([DelaminationEllipse(1, (5, 2.5), 3, 1.5, 0)], dent_depth_mm=0.2)
    sigma_pristine, _, _ = fe3d_cai_buckling(small_cfg, DamageState([], 0.0), lam, 500.0)
    sigma_damaged, _, _ = fe3d_cai_buckling(small_cfg, ds, lam, 500.0)
    assert 0 < sigma_damaged <= sigma_pristine


def test_fe3d_cai_buckling_sublaminate_path_binds_on_thin_panel():
    """The damaged-vs-pristine ordering on ``small_cfg`` saturates against
    ``sigma_pristine_MPa`` (the small coupon's plate buckling load is huge),
    so the test above can't verify the sublaminate path actually engages.
    This test uses a wider panel with a large delamination so the
    sublaminate buckling stress is strictly below pristine and strictly
    below the full-panel buckling stress — proving min(panel, sublam) picks
    the sublaminate when it should.
    """
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0, boundary="simply_supported"),
        loading="compression",
        tier="fe3d",
        damage=DamageState(
            [DelaminationEllipse(3, (75.0, 50.0), 40.0, 30.0, 0.0)], dent_depth_mm=0.5
        ),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0),
    )
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], cfg.layup_deg, cfg.ply_thickness_mm)
    sigma_pristine_panel, _, _ = fe3d_cai_buckling(
        cfg, DamageState([], 0.0), lam, sigma_pristine_MPa=600.0
    )
    sigma_damaged, lambda_crit, notes = fe3d_cai_buckling(
        cfg, cfg.damage, lam, sigma_pristine_MPa=600.0
    )
    assert sigma_damaged < sigma_pristine_panel
    assert sigma_damaged < 600.0  # neither path saturates
    # The onset note names the sublaminate it comes from (#133).
    onset = [n for n in notes if n.startswith("fe3d buckling onset")]
    assert onset == [
        f"fe3d buckling onset {lambda_crit:.1f} MPa (sublaminate over interface 3; "
        f"buckling_eigenvalues) is from the closed-form Rayleigh-Ritz kernel "
        f"semi_analytical uses, not a 3D eigensolve (#129)"
    ], notes


def test_fe3d_cai_buckling_canonical_panel_value():
    """Regression test pinning the buckling stress on the canonical
    200x150 mm 8-ply quasi-isotropic IM7/8552 pristine reference (issue #129).

    The simply supported Rayleigh-Ritz with D16 / D26 converges to 12.204
    MPa here (#132 PR0; the Navier closed form that ignores them is 12.582
    MPa, 3.0% higher). A regression to a different formulation (or a
    unit-conversion bug) would move this number by an order of magnitude.
    """
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(200.0, 150.0, boundary="simply_supported"),
        loading="compression",
        tier="fe3d",
        damage=DamageState([], dent_depth_mm=0.0),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0),
    )
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], cfg.layup_deg, cfg.ply_thickness_mm)
    sigma, _, _ = fe3d_cai_buckling(cfg, cfg.damage, lam, sigma_pristine_MPa=887.5)
    assert sigma == pytest.approx(12.204, abs=1e-3), f"expected ~12.204 MPa, got {sigma}"


def test_bvid_analysis_fe3d_cai_uses_buckling_path(small_cfg):
    """BvidAnalysis with tier=fe3d should populate buckling_eigenvalues."""
    with _tiny_panel_impact_warnings():
        r = BvidAnalysis(small_cfg).run()
    assert r.tier_used == "fe3d"
    assert r.buckling_eigenvalues is not None
    assert len(r.buckling_eigenvalues) >= 1
    assert r.buckling_eigenvalues[0] > 0
    # Issue #133: the success path says in-band that the buckling channel is
    # the closed form, and which CAI channel sets the strength.
    assert any(
        n.startswith("fe3d buckling onset") and "closed-form" in n.lower() for n in r.notes
    ), r.notes
    assert any(n.startswith("fe3d CAI: ") and " governs " in n for n in r.notes), r.notes


def test_fe3d_cai_buckling_degenerate_closed_form_emits_note(small_cfg, monkeypatch):
    """If the Rayleigh-Ritz closed-form returns a degenerate (non-finite or
    non-positive) value for both the full-panel and the worst-sublaminate
    paths, fe3d_cai_buckling falls back to pristine AND emits a runtime
    note. Without the note the user would see knockdown=1.0 with no
    indication that the buckling channel was inactive."""
    from bvidfe.analysis import fe_tier as ft
    from bvidfe.core.laminate import Laminate as Lam

    monkeypatch.setattr(ft, "panel_buckling_load", lambda *_a, **_kw: float("inf"))
    monkeypatch.setattr(ft, "sublaminate_buckling_load", lambda *_a, **_kw: float("inf"))
    lam = Lam(MATERIAL_LIBRARY["IM7/8552"], small_cfg.layup_deg, small_cfg.ply_thickness_mm)
    sigma, lambda_crit, notes = fe3d_cai_buckling(
        small_cfg, DamageState([], 0.0), lam, sigma_pristine_MPa=500.0
    )
    assert sigma == 500.0  # fell back to pristine
    assert lambda_crit == 0.0
    assert any("rayleigh-ritz" in n.lower() or "degenerate" in n.lower() for n in notes), notes


def test_bvid_analysis_fe3d_residual_ignores_buckling_onset(small_cfg, monkeypatch):
    """Buckling onset is reported in ``buckling_eigenvalues`` but does not
    set the fe3d residual: a buckled sublaminate keeps carrying load until
    its delamination grows, and whole-panel buckling is suppressed by a CAI
    fixture's anti-buckling guides. A degenerate closed form therefore
    leaves the residual unchanged and raises no fallback tag."""
    from bvidfe.analysis import fe_tier as ft

    with _tiny_panel_impact_warnings():
        baseline = BvidAnalysis(small_cfg).run()
    monkeypatch.setattr(ft, "panel_buckling_load", lambda *_a, **_kw: float("inf"))
    monkeypatch.setattr(ft, "sublaminate_buckling_load", lambda *_a, **_kw: float("inf"))
    with _tiny_panel_impact_warnings():
        r = BvidAnalysis(small_cfg).run()
    assert r.residual_strength_MPa == baseline.residual_strength_MPa
    assert r.buckling_eigenvalues is None
    assert r.warnings == baseline.warnings  # no extra fallback warning
