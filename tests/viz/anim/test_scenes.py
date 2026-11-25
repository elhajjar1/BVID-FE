"""The two scenes: timeline invariants (model outputs only, at the keyframes)
and frames (exact 1080p RGBA, not blank, with the citation footer)."""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from matplotlib.text import Text

from bvidfe import BvidAnalysis, ImpactEvent, PanelGeometry, load_cscan_json
from bvidfe.citation import frame_credit, short_credit
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.viz.anim import frame as fr
from bvidfe.viz.anim.scenes import damage_stack, energy_sweep
from bvidfe.viz.anim.scenes._common import damage_region, onset_energy_J
from bvidfe.viz.anim.spec import AnimationSpec
from tests.viz.anim._anim_helpers import LAYUP_8, small_config

TIMES = (0.0, 0.5, 1.0)
SAMPLE_CSCAN = Path(__file__).resolve().parents[3] / "examples" / "sample_cscan.json"


def _texts(fig) -> list[str]:
    return [t.get_text() for t in fig.findobj(Text)]


def _render(tl, t: float, spec: AnimationSpec | None = None) -> tuple[np.ndarray, list[str]]:
    fig = fr.new_figure(spec)
    tl.draw(fig, t)
    return fr.render_rgba(fig), [*_texts(fig), *fr.find_text(fig, fr.CITATION_GID)]


# ---- A1: energy sweep timeline ------------------------------------------------------


def test_sweep_keyframes_are_bvid_analysis_outputs(sweep_timeline):
    tl = sweep_timeline
    cfg = tl.config
    for k in (0, len(tl.keyframes) // 2, len(tl.keyframes) - 1):
        kf = tl.keyframes[k]
        r = BvidAnalysis(replace(cfg, impact=replace(cfg.impact, energy_J=kf.energy_J))).run()
        assert kf.dpa_mm2 == r.dpa_mm2
        assert kf.dent_mm == r.damage.dent_depth_mm
        assert kf.knockdown == r.knockdown
        assert kf.residual_MPa == r.residual_strength_MPa
        assert kf.pristine_MPa == r.pristine_strength_MPa
        assert list(kf.delaminations) == r.damage.delaminations


def test_sweep_physics_invariants(sweep_timeline):
    tl = sweep_timeline
    e = tl.energies
    dpa, dent, kd = tl.series("dpa_mm2"), tl.series("dent_mm"), tl.series("knockdown")
    assert np.all(np.diff(e) > 0)
    assert np.all(np.diff(dpa) >= 0) and np.all(np.diff(dent) >= 0)
    assert np.all(np.diff(kd) <= 1e-12)
    below = e <= tl.onset_J
    assert below.any() and (~below).any()
    assert np.all(dpa[below] == 0) and np.all(dent[below] == 0)
    assert np.all(kd[below] == 1.0)
    assert all(not tl.keyframes[i].delaminations for i in np.flatnonzero(below))
    # The default grid starts at 0.8 x onset, ends at the config's energy and
    # holds the onset energy itself.
    assert e[0] == pytest.approx(0.8 * tl.onset_J)
    assert e[-1] == pytest.approx(20.0)
    assert np.any(e == tl.onset_J)
    assert tl.onset_J == onset_energy_J(tl.config)


def test_sweep_energy_at_follows_the_schedule(sweep_timeline):
    tl = sweep_timeline
    e = tl.energies
    assert tl.energy_at(0.0) == e[0] and tl.energy_at(1.0) == e[-1]
    ts = np.linspace(0, 1, 50)
    assert np.all(np.diff([tl.energy_at(t) for t in ts]) >= 0)


def test_sweep_explicit_energies_get_the_onset():
    cfg = small_config()
    onset = onset_energy_J(cfg)
    tl = energy_sweep.build_timeline(cfg, energies=[10.0, 0.0, 5.0, 5.0])
    assert list(tl.energies) == sorted([0.0, onset, 5.0, 10.0])


@pytest.mark.parametrize(
    "kw",
    [
        {"energies": [1.0]},
        {"energies": [1.0, np.nan]},
        {"energies": [-1.0, 2.0]},
        {"energies": [2.0, 2.0]},
        {"n_energies": 1},
        {"n_energies": 2.5},
        {"energy_min": 10.0, "energy_max": 5.0},
    ],
)
def test_sweep_rejects_bad_energies(kw):
    with pytest.raises(ValueError):
        energy_sweep.build_timeline(small_config(), **kw)


def test_sweep_refuses_fe3d_and_inspection_configs():
    with pytest.raises(ValueError, match="too slow"):
        energy_sweep.build_timeline(small_config(tier="fe3d"))
    damage = DamageState([DelaminationEllipse(0, (75.0, 50.0), 10.0, 5.0, 0.0)])
    with pytest.raises(ValueError, match="impact-driven"):
        energy_sweep.build_timeline(small_config(impact=None, damage=damage))
    with pytest.raises(ValueError, match="AnalysisConfig"):
        energy_sweep.build_timeline("IM7/8552")  # type: ignore[arg-type]


def test_sweep_default_energy_max_above_onset():
    # An energy at or below the onset would sweep nothing: 4 x onset instead.
    cfg = small_config(impact=ImpactEvent(energy_J=0.5))
    grid = energy_sweep.energy_grid(cfg, n_energies=5)
    assert grid[-1] == pytest.approx(4.0 * onset_energy_J(cfg))


def test_semi_analytical_sweep():
    tl = energy_sweep.build_timeline(small_config(tier="semi_analytical"), n_energies=3)
    assert tl.tier == "semi_analytical" and len(tl.keyframes) == 4


def test_ellipses_fade_in_between_keyframes():
    e = DelaminationEllipse(2, (75.0, 50.0), 20.0, 10.0, 30.0)
    half = energy_sweep.interpolate_ellipses([], [e], 0.5)
    assert len(half) == 1
    d = half[0]
    assert (d.major_mm, d.minor_mm, d.orientation_deg) == (10.0, 5.0, 30.0)
    assert 0 < d.fade < 1
    assert energy_sweep.interpolate_ellipses([], [e], 0.0) == []
    full = energy_sweep.interpolate_ellipses([e], [e], 0.3)[0]
    assert (full.major_mm, full.fade) == (20.0, 1.0)


def test_zoom_only_for_small_damage():
    panel = PanelGeometry(150.0, 100.0)
    small = [DelaminationEllipse(0, (75.0, 50.0), 5.0, 2.0, 0.0)]
    big = [DelaminationEllipse(0, (75.0, 50.0), 40.0, 20.0, 0.0)]
    box, zoomed = damage_region(panel, small, (75.0, 50.0), 1.5)
    assert zoomed and box[0] < 70 < 80 < box[2]
    assert damage_region(panel, big, (75.0, 50.0), 1.5) == ((0.0, 0.0, 150.0, 100.0), False)
    assert damage_region(panel, [], None, 1.5)[1] is False


# ---- A2: damage stack timeline ---------------------------------------------------------


def test_stack_from_config_runs_the_impact_mapping(stack_timeline):
    tl = stack_timeline
    r = BvidAnalysis(tl.config).run()
    assert tl.damage.delaminations == r.damage.delaminations
    assert tl.damage.dent_depth_mm == r.damage.dent_depth_mm
    assert tl.n_plies == 8 and tl.layup_deg == tuple(float(a) for a in LAYUP_8)
    assert tl.onset_J == pytest.approx(onset_energy_J(tl.config))
    assert tl.impact_point_mm == (75.0, 50.0)


def test_stack_from_a_cscan():
    damage = load_cscan_json(SAMPLE_CSCAN)
    tl = damage_stack.build_timeline(damage, layup_deg=LAYUP_8, ply_thickness_mm=0.152)
    assert tl.n_plies == 8 and tl.config is None and tl.source == "supplied damage state"
    assert tl.panel.Lx_mm == 150.0 and tl.panel.Ly_mm == 100.0
    bare = damage_stack.build_timeline(damage, intro=False)
    assert bare.layup_deg is None and bare.n_plies == 8 and not bare.intro
    with pytest.raises(ValueError, match="past the last interface"):
        damage_stack.build_timeline(damage, layup_deg=[0, 90, 0])
    with pytest.raises(ValueError, match="come from the config"):
        damage_stack.build_timeline(small_config(), layup_deg=LAYUP_8)
    with pytest.raises(ValueError):
        damage_stack.build_timeline("damage")  # type: ignore[arg-type]


def test_fibre_segments_stay_in_the_box():
    box = (10.0, 20.0, 40.0, 30.0)
    for angle in (0.0, 45.0, -45.0, 90.0, 30.0):
        segs = damage_stack.fibre_segments(angle, box, 3.0)
        assert segs
        pts = np.array([p for s in segs for p in s])
        assert pts[:, 0].min() >= 10 - 1e-9 and pts[:, 0].max() <= 40 + 1e-9
        assert pts[:, 1].min() >= 20 - 1e-9 and pts[:, 1].max() <= 30 + 1e-9


def test_far_corner_carries_the_swatch():
    # Camera on -y (azimuth -90 deg): the far corners are the two at y = ly.
    w = damage_stack.corner_weights(150.0, 100.0, -90.0)
    assert w[2] == pytest.approx(1.0) and w[3] == pytest.approx(1.0)
    assert w[0] == 0.0 and w[1] == 0.0


# ---- Frames ----------------------------------------------------------------------------


@pytest.mark.parametrize("scene", ["sweep", "stack"])
@pytest.mark.parametrize("t", TIMES)
def test_frame_is_1080p_rgba_with_citation(scene, t, request):
    tl = request.getfixturevalue(f"{scene}_timeline")
    rgba, texts = _render(tl, t)
    assert rgba.shape == (1080, 1920, 4) and rgba.dtype == np.uint8
    assert (rgba[..., 3] == 255).all()
    assert rgba[..., :3].std() > 10  # not blank
    assert frame_credit() in texts
    credit = frame_credit()
    assert "please cite" in credit and "Elhajjar" in credit and "github.com" in credit


def test_sweep_frame_labels(sweep_timeline):
    _, texts = _render(sweep_timeline, 0.0)
    assert any(s.startswith("below damage onset (") for s in texts)
    assert "empirical tier" in " ".join(texts)
    _, texts = _render(sweep_timeline, 1.0)
    assert not any(s.startswith("below damage onset") for s in texts)


def test_biaxial_frames_say_experimental():
    cfg = small_config(transverse_stress_ratio=0.5)
    tl = energy_sweep.build_timeline(cfg, n_energies=3)
    fig = fr.new_figure()
    tl.draw(fig, 0.5)
    caveats = fr.find_text(fig, fr.CAVEAT_GID)
    assert any("experimental, not validated" in c for c in caveats)
    stack = damage_stack.build_timeline(cfg)
    stack.draw(fig, 0.5)
    assert any("experimental, not validated" in c for c in fr.find_text(fig, fr.CAVEAT_GID))


def test_below_onset_stack_says_so():
    tl = damage_stack.build_timeline(small_config(impact=ImpactEvent(energy_J=0.2)))
    assert not tl.damage.delaminations
    _, texts = _render(tl, 0.8)
    assert any("below damage onset" in s for s in texts)


@pytest.mark.parametrize("scene", ["sweep", "stack"])
def test_frames_are_deterministic(scene, request):
    tl = request.getfixturevalue(f"{scene}_timeline")
    a, _ = _render(tl, 0.37)
    b, _ = _render(tl, 0.37)
    assert hashlib.sha256(a.tobytes()).hexdigest() == hashlib.sha256(b.tobytes()).hexdigest()


def test_frame_at_another_16_9_size(sweep_timeline):
    rgba, texts = _render(sweep_timeline, 0.5, AnimationSpec(width=640))
    assert rgba.shape == (360, 640, 4)
    # A preview-size frame carries the short credit, in type it can show.
    assert short_credit() in texts and frame_credit() not in texts
    _, texts = _render(sweep_timeline, 0.5, AnimationSpec(width=1280))
    assert frame_credit() in texts


def test_preview_frames_writes_pngs(tmp_path, stack_timeline):
    from PIL import Image

    from bvidfe.viz.anim import preview_frames

    frames = preview_frames(stack_timeline, (0.0, 1.0), out_dir=tmp_path, stem="s")
    assert [f.shape for f in frames] == [(1080, 1920, 4)] * 2
    for name in ("s_t0.00.png", "s_t1.00.png"):
        with Image.open(tmp_path / name) as im:
            assert im.size == (1920, 1080)


@pytest.mark.parametrize("n", [1, 4])
def test_caveats_fit_the_frame(n):
    """Caveats that do not fit one row are all set smaller, and every one
    is drawn inside the margins; one that fits keeps the 18 px type."""
    fig = fr.new_figure()
    notes = [f"caveat number {k}: experimental, not validated, qualitative only" for k in range(n)]
    fr.draw_caveats(fig, notes)
    drawn = [t for t in fig.texts if t.get_gid() == fr.CAVEAT_GID]
    assert [t.get_text() for t in drawn] == notes
    fig.canvas.draw()
    right = max(t.get_window_extent().x1 for t in drawn) * 100.0 / fig.dpi
    assert right <= fr.W - fr.MARGIN + 1
    sizes = {round(t.get_fontsize(), 6) for t in drawn}
    assert sizes == ({round(fr.pt(18), 6)} if n == 1 else {max(sizes)})
    if n == 4:
        assert max(sizes) < fr.pt(18)


def test_layup_label():
    assert fr.layup_label([0, 45, -45, 90] * 4) == "[0/45/-45/90]4"
    assert fr.layup_label(LAYUP_8) == "[0/45/-45/90]s"
    assert fr.layup_label([0, 45, -45, 90] * 2 + [90, -45, 45, 0] * 2) == "[0/45/-45/90]2s"
    assert fr.layup_label([0, 90, 45]) == "[0/90/45]"
    long = list(range(0, 15 * 7, 7))
    assert fr.layup_label(long).endswith("(15 plies)")
