"""A3, the fe3d load ramp: timeline invariants (the ramp ends where
``run()`` puts the residual strength, the first element to reach a load
ratio of 1 is fe3d's first-ply failure), frames, the face cap, the CLI and
the encode. Small panels at 5 and 10 mm keep the fe3d solves under a
second; the full 10 s render of a 150 x 100 mm, 16-ply panel is slow."""

from __future__ import annotations

import hashlib
from dataclasses import replace
from unittest import mock

import numpy as np
import pytest
from matplotlib.text import Text
from PIL import Image

from bvidfe import (
    AnalysisConfig,
    BvidAnalysis,
    ImpactEvent,
    MeshParams,
    PanelGeometry,
    load_cscan_json,
)
from bvidfe.citation import frame_credit
from bvidfe.viz.anim import (
    AnimationSpec,
    encode_timeline,
    estimate_render_seconds,
    list_animations,
    render_cai_ramp,
)
from bvidfe.viz.anim import frame as fr
from bvidfe.viz.anim.cli import main
from bvidfe.viz.anim.scenes import cai_ramp
from tests.viz.anim._anim_helpers import LAYUP_8, needs_ffmpeg, probe
from tests.viz.anim.test_scenes import SAMPLE_CSCAN

TIMES = (0.0, 0.5, 1.0)


def _cfg(**kw) -> AnalysisConfig:
    """8-ply QI IM7/8552, 60 x 40 mm, 5 mm fe3d mesh, a 3 J impact."""
    base = dict(
        material="IM7/8552",
        layup_deg=LAYUP_8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(60.0, 40.0),
        tier="fe3d",
        impact=ImpactEvent(energy_J=3.0),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=5.0),
    )
    base.update(kw)
    return AnalysisConfig(**base)


# Delamination growth governs (as in tests/analysis/test_fe_fields.py).
GROWTH = dict(
    layup_deg=[45, -45, 90, 90, -45, 45],
    panel=PanelGeometry(150.0, 100.0),
    impact=ImpactEvent(energy_J=5.0),
    mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0),
)


@pytest.fixture(scope="module")
def ramp():
    return cai_ramp.build_timeline(_cfg())


@pytest.fixture(scope="module")
def growth_ramp():
    return cai_ramp.build_timeline(_cfg(**GROWTH))


def _texts(fig) -> list[str]:
    return [t.get_text() for t in fig.findobj(Text)]


def _render(tl, t: float) -> tuple[np.ndarray, list[str]]:
    fig = fr.new_figure()
    tl.draw(fig, t)
    return fr.render_rgba(fig), _texts(fig)


def _ramp_end(tl) -> float:
    return tl.schedule.bounds("ramp")[1]


# ---- Timeline --------------------------------------------------------------------------


def test_ramp_ends_at_first_ply_failure_where_run_puts_it(ramp):
    f = ramp.fields
    r = BvidAnalysis(ramp.config).run()
    assert f.governing == "first_ply_failure" and f.fpf_knockdown < 1
    assert f.residual_strength_MPa == r.residual_strength_MPa
    assert f.knockdown == r.knockdown
    end = _ramp_end(ramp)
    assert ramp.stress_at(end) == f.damaged.fpf_stress_MPa == ramp.end_stress_MPa
    assert ramp.stress_at(1.0) == ramp.end_stress_MPa
    assert ramp.reported_stress_at(1.0) == r.residual_strength_MPa
    # The ramp's end on the reported scale is the residual strength.
    assert f.reported_stress(ramp.end_stress_MPa) == pytest.approx(
        r.residual_strength_MPa, rel=1e-12
    )


def test_load_ratio_reaches_one_exactly_at_first_ply_failure(ramp):
    f = ramp.fields
    end = _ramp_end(ramp)
    ratio = ramp.load_ratio(ramp.stress_at(end))
    assert ratio.max() == pytest.approx(1.0, abs=1e-9)
    assert int(np.argmax(ratio)) == f.damaged.fpf_element
    assert ratio.max() >= cai_ramp.FAILED
    # Rising through the ramp, below 1 until its end.
    ts = np.linspace(0.0, end, 40)
    peaks = [ramp.max_ratio(ramp.stress_at(t)) for t in ts]
    assert peaks[0] == 0.0
    assert np.all(np.diff(peaks) >= 0)
    assert max(peaks[:-1]) < cai_ramp.FAILED
    # The undamaged panel at the same applied stress: the fpf knockdown.
    assert ramp.max_ratio(ramp.end_stress_MPa, damaged=False) == pytest.approx(
        f.fpf_knockdown, rel=1e-9
    )


def test_displacement_is_zero_at_the_start(ramp):
    d = ramp.fields.damaged
    assert ramp.stress_at(0.0) == 0.0
    assert not np.any(d.displacement(ramp.stress_at(0.0)))
    nodes = cai_ramp._drawn_points(ramp, d.mesh.node_coords, d.u_ref.reshape(-1, 3), 0.0)
    h = d.mesh.node_coords[:, 2].max()
    assert np.array_equal(nodes[:, :2], d.mesh.node_coords[:, :2])
    assert np.allclose(nodes[:, 2], (h - d.mesh.node_coords[:, 2]) * ramp.view.thickness_scale)
    assert d.end_shortening_mm(ramp.stress_at(0.0)) == 0.0


def test_growth_governed_ramp_ends_at_the_growth_stress(growth_ramp):
    tl = growth_ramp
    f = tl.fields
    r = BvidAnalysis(tl.config).run()
    assert f.governing == "delamination_growth"
    assert f.residual_strength_MPa == r.residual_strength_MPa == f.growth_stress_MPa
    assert tl.end_stress_MPa == pytest.approx(f.applied_stress(f.growth_stress_MPa), rel=1e-12)
    assert tl.max_ratio(tl.end_stress_MPa) < 1.0  # no element fails first
    assert tl.reported_stress_at(1.0) == r.residual_strength_MPa
    _, texts = _render(tl, 1.0)
    joined = " ".join(texts)
    assert "Delamination growth governs" in joined and "no FE field" in joined
    assert f"{r.residual_strength_MPa:.1f} MPa" in texts
    assert f"{r.knockdown:.3f}" in texts
    assert "first-ply failure (not reached)" in texts


def test_view_draws_the_block_on_one_side_of_the_cut(ramp):
    v = ramp.view
    d = ramp.fields.damaged
    nz, ny, nx = d.grid_shape
    k, j, i = cai_ramp._unravel(d.fpf_element, d.grid_shape)
    assert v.cut_row == j
    _, rows, cols = np.unravel_index(v.parent, (nz, ny, nx))
    assert np.all(rows >= j) if v.keep_high else np.all(rows <= j)
    # The cut face holds every layer of the cut row, the first element too.
    cut = v.parent[v.n_top : v.n_top + v.n_cut]
    assert len(cut) == nz * nx and d.fpf_element in cut
    assert v.n_faces <= cai_ramp.face_count_bound(ramp.config)
    assert v.faces.min() >= 0 and v.faces.max() < d.mesh.n_nodes
    assert not v.faces.flags.writeable
    assert v.deformation_scale >= 1 and v.thickness_scale >= 1


def test_pointwise_ratio_and_bad_ratio(ramp):
    pw = cai_ramp.CaiRampTimeline(
        config=ramp.config,
        fields=ramp.fields,
        view=cai_ramp._view(ramp.fields, "pointwise"),
        ratio="pointwise",
    )
    # Without a notch the pointwise check fails no later than the averaged one.
    assert pw.max_ratio(pw.end_stress_MPa) >= 1.0 - 1e-9
    with pytest.raises(ValueError, match="averaged"):
        cai_ramp.build_timeline(_cfg(), ratio="gauss")  # type: ignore[arg-type]


def test_face_cap_is_checked_before_solving():
    cfg = _cfg(panel=PanelGeometry(150.0, 100.0), mesh=MeshParams(1, 1.0))
    assert cai_ramp.face_count_bound(cfg) > cai_ramp.MAX_FACES
    with mock.patch.object(cai_ramp, "fe3d_field_history") as solve:
        with pytest.raises(ValueError, match="faces a frame"):
            cai_ramp.build_timeline(cfg)
        with pytest.raises(ValueError, match="in_plane_size_mm"):
            cai_ramp.build_timeline(_cfg(), max_faces=10)
    solve.assert_not_called()
    with pytest.raises(ValueError, match="AnalysisConfig"):
        cai_ramp.build_timeline("IM7/8552")  # type: ignore[arg-type]


def test_ratio_label_rounds_down_below_failure():
    assert cai_ramp.ratio_label(0.996) == "0.99"
    assert cai_ramp.ratio_label(0.0) == "0.00"
    assert cai_ramp.ratio_label(1.0 - 1e-12) == "1.00"
    assert cai_ramp.ratio_label(1.0) == "1.00"
    assert cai_ramp.ratio_label(1.237) == "1.24"


def test_render_estimate(ramp):
    spec = AnimationSpec()
    est = estimate_render_seconds(ramp, spec)
    assert est == pytest.approx(300 * (ramp.frame_seconds_estimate() + 0.04))
    assert estimate_render_seconds(ramp, spec, frame_seconds=1.0) == pytest.approx(300 * 1.04)
    short = AnimationSpec(duration_s=1.0, fps=10)
    assert estimate_render_seconds(ramp, short) == pytest.approx(est / 30)
    gif = AnimationSpec(format="gif")
    assert estimate_render_seconds(ramp, gif) < est
    assert "cai-ramp" in list_animations()


# ---- Frames ----------------------------------------------------------------------------


@pytest.mark.parametrize("t", TIMES)
def test_frame_is_1080p_with_citation(ramp, t):
    rgba, texts = _render(ramp, t)
    assert rgba.shape == (1080, 1920, 4) and rgba.dtype == np.uint8
    assert (rgba[..., 3] == 255).all()
    assert rgba[..., :3].std() > 10
    assert frame_credit() in texts


def test_numbers_on_screen_are_run_s(ramp):
    r = BvidAnalysis(ramp.config).run()
    _, texts = _render(ramp, 1.0)
    residual = f"{r.residual_strength_MPa:.1f} MPa"
    assert texts.count(residual) >= 2  # the title readout and the result block
    assert f"{r.knockdown:.3f}" in texts
    assert f"{r.pristine_strength_MPa:.1f} MPa" in texts
    assert "1.00" in texts  # the damaged panel's highest load ratio
    joined = " ".join(texts)
    assert "First-ply failure at " + residual in joined
    assert "first-ply failure: ply " in joined
    _, texts = _render(ramp, 0.0)
    assert "0.0 MPa" in texts and "0.00" in texts
    assert not any(s.startswith("First-ply failure") for s in texts)


def test_frame_labels_the_caveats_and_scales(ramp):
    fig = fr.new_figure()
    ramp.draw(fig, 0.5)
    caveats = fr.find_text(fig, fr.CAVEAT_GID)
    assert any("linear elastic to first-ply failure, no progressive damage" in c for c in caveats)
    assert any("fe3d: qualitative" in c for c in caveats)
    joined = " ".join(_texts(fig))
    v = ramp.view
    assert f"Displacements ×{v.deformation_scale:g}, thickness ×{v.thickness_scale:g}" in joined
    assert "fe3d mesh 5 mm, 1 element per ply, 768 elements" in " ".join(
        fr.find_text(fig, fr.SUMMARY_GID)
    )
    assert "fe3d tier" in joined


@pytest.mark.parametrize(
    "kw, caveat",
    [
        ({"loading": "tension"}, "TAI: not validated"),
        ({"transverse_stress_ratio": 0.5}, "experimental, not validated"),
    ],
)
def test_tai_and_biaxial_frames(kw, caveat):
    tl = cai_ramp.build_timeline(_cfg(**kw))
    r = BvidAnalysis(tl.config).run()
    assert tl.fields.residual_strength_MPa == r.residual_strength_MPa
    fig = fr.new_figure()
    tl.draw(fig, 1.0)
    assert any(caveat in c for c in fr.find_text(fig, fr.CAVEAT_GID))
    texts = _texts(fig)
    assert f"{r.residual_strength_MPa:.1f} MPa" in texts
    if kw.get("loading") == "tension":
        assert any(s.startswith("Tension after impact") for s in texts)
        assert "End extension (mm)" in texts


def test_below_onset_and_cscan_frames():
    tl = cai_ramp.build_timeline(_cfg(impact=ImpactEvent(energy_J=0.3)))
    assert not tl.fields.damage.delaminations
    assert tl.fields.knockdown == 1.0
    _, texts = _render(tl, 1.0)
    assert any(s.startswith("Below damage onset") for s in texts)
    cfg = _cfg(impact=None, damage=load_cscan_json(SAMPLE_CSCAN), panel=PanelGeometry(150, 100))
    tl = cai_ramp.build_timeline(cfg)
    assert tl.fields.residual_strength_MPa == BvidAnalysis(cfg).run().residual_strength_MPa
    _, texts = _render(tl, 0.5)
    assert any("damage as supplied" in s for s in texts)


def test_frames_are_deterministic(ramp):
    a, _ = _render(ramp, 0.37)
    b, _ = _render(ramp, 0.37)
    assert hashlib.sha256(a.tobytes()).hexdigest() == hashlib.sha256(b.tobytes()).hexdigest()


def test_tier_is_not_read(ramp):
    tl = cai_ramp.build_timeline(replace(ramp.config, tier="empirical"))
    assert tl.fields.residual_strength_MPa == ramp.fields.residual_strength_MPa
    assert tl.tier == "fe3d"


# ---- Encode and CLI --------------------------------------------------------------------


def test_mp4_smoke_encode(ramp, tmp_path):
    needs_ffmpeg()
    result = encode_timeline(ramp, tmp_path / "ramp.mp4", AnimationSpec(duration_s=1.0, fps=4))
    info = probe(result.path)
    assert info["codec"] == "h264" and (info["width"], info["height"]) == (1920, 1080)
    assert info["frames"] == 4 and info["pix_fmt"] == "yuv420p"


def test_render_cai_ramp_gif(tmp_path):
    progress = []
    tl = render_cai_ramp(
        _cfg(),
        tmp_path / "r.gif",
        AnimationSpec(format="gif", duration_s=1, fps=2, width=320),
        progress=lambda i, n: progress.append(i),
    )
    assert progress == [1, 2] and isinstance(tl, cai_ramp.CaiRampTimeline)
    with Image.open(tmp_path / "r.gif") as im:
        assert im.size == (320, 180)


CLI_ARGS = [
    "cai-ramp",
    "--material",
    "IM7/8552",
    "--layup",
    "0,45,-45,90,90,-45,45,0",
    "--thickness",
    "0.152",
    "--panel",
    "60x40",
    "--energy",
    "3",
]


def test_cli_preview(tmp_path, capsys):
    rc = main([*CLI_ARGS, "--in-plane-size", "10", "--preview", "--out", str(tmp_path / "p.mp4")])
    assert rc == 0
    printed = capsys.readouterr().out.split()
    assert printed == [str(tmp_path / f"p_t{t}.png") for t in ("0.00", "0.50", "1.00")]
    for p in printed:
        with Image.open(p) as im:
            assert im.size == (1920, 1080)


def test_cli_errors(tmp_path, capsys):
    rc = main([*CLI_ARGS[:-3], "150x100", "--energy", "3", "--in-plane-size", "0.5"])
    assert rc == 1
    assert "faces a frame" in capsys.readouterr().err
    with pytest.raises(SystemExit) as exc:
        main([*CLI_ARGS, "--tier", "fe3d"])
    assert exc.value.code == 2
    assert "unrecognized arguments: --tier" in capsys.readouterr().err


# ---- Full render (slow) ----------------------------------------------------------------


@pytest.mark.slow
def test_full_10s_1080p_render(tmp_path):
    needs_ffmpeg()
    cfg = _cfg(
        layup_deg=[0, 45, -45, 90] * 4,
        panel=PanelGeometry(150.0, 100.0),
        impact=ImpactEvent(energy_J=20.0),
    )
    path = tmp_path / "ramp.mp4"
    tl = render_cai_ramp(cfg, path)
    info = probe(path)
    assert info["codec"] == "h264" and (info["width"], info["height"]) == (1920, 1080)
    assert info["fps"] == pytest.approx(30.0) and info["frames"] == 300
    assert 0.2e6 < path.stat().st_size < 30e6
    assert tl.fields.residual_strength_MPa == BvidAnalysis(cfg).run().residual_strength_MPa
