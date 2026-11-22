"""AnimationSpec validation, the timeline schedule and easing, and the
shared style helpers."""

from __future__ import annotations

import math

import numpy as np
import pytest

from bvidfe.viz.anim.spec import AnimationSpec
from bvidfe.viz.anim.timeline import (
    EASINGS,
    Phase,
    Schedule,
    bracket,
    lerp,
    nearest,
)
from bvidfe.viz.style import CATEGORICAL, series_color, series_linestyle

# ---- AnimationSpec ---------------------------------------------------------------


def test_spec_defaults_are_1080p_h264():
    s = AnimationSpec()
    assert s.size == (1920, 1080)
    assert (s.fps, s.duration_s, s.codec, s.crf, s.pix_fmt, s.format) == (
        30,
        10.0,
        "libx264",
        18,
        "yuv420p",
        "mp4",
    )
    assert s.n_frames == 300
    t = s.frame_times()
    assert t[0] == 0.0 and t[-1] == 1.0 and len(t) == 300


def test_spec_gif_defaults_to_a_small_preview():
    s = AnimationSpec(format="gif")
    assert s.size == (640, 360) and s.fps == 12


def test_spec_one_side_gives_the_other():
    assert AnimationSpec(width=1280).size == (1280, 720)
    assert AnimationSpec(height=720).size == (1280, 720)


def test_spec_dpi_gives_exact_pixels():
    for w in (640, 960, 1280, 1920, 3840):
        s = AnimationSpec(width=w)
        assert int(19.2 * s.dpi) == w
        assert int(10.8 * s.dpi) == s.height


def test_spec_at_least_two_frames():
    assert AnimationSpec(duration_s=0.01, fps=1).n_frames == 2


@pytest.mark.parametrize(
    "kw",
    [
        {"width": 1920, "height": 1200},  # not 16:9
        {"width": 1922, "height": 1081},  # odd
        {"width": 0},
        {"width": 1920.0, "height": 1080},
        {"width": True},
        {"fps": 0},
        {"fps": 121},
        {"fps": 29.97},
        {"duration_s": 0},
        {"duration_s": -1},
        {"duration_s": math.nan},
        {"duration_s": math.inf},
        {"duration_s": True},
        {"duration_s": 1e6},
        {"crf": -1},
        {"crf": 52},
        {"crf": 18.0},
        {"codec": ""},
        {"codec": "-y"},
        {"pix_fmt": 3},
        {"format": "webm"},
    ],
)
def test_spec_rejects(kw):
    with pytest.raises(ValueError):
        AnimationSpec(**kw)


# ---- Timeline --------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(EASINGS))
def test_easings_fix_ends_and_rise(name):
    f = EASINGS[name]
    assert f(0.0) == pytest.approx(0.0) and f(1.0) == pytest.approx(1.0)
    u = np.linspace(0, 1, 101)
    assert np.all(np.diff([f(x) for x in u]) >= -1e-12)


def test_schedule_phases_and_progress():
    s = Schedule((Phase("a", 1.0), Phase("b", 2.0, "ease_in_out"), Phase("c", 1.0)))
    assert s.bounds("b") == pytest.approx((0.25, 0.75))
    assert s.at(0.0) == ("a", 0.0)
    assert s.at(0.5)[0] == "b" and s.at(0.5)[1] == pytest.approx(0.5)
    assert s.at(1.0) == ("c", 1.0)
    assert s.progress(0.1, "b") == 0.0 and s.progress(0.9, "b") == 1.0
    assert s.progress(-1.0, "a") == 0.0 and s.progress(2.0, "c") == 1.0
    with pytest.raises(KeyError):
        s.bounds("missing")
    with pytest.raises(ValueError):
        s.at(math.nan)


@pytest.mark.parametrize(
    "phases",
    [(), (Phase("a", 1.0), Phase("a", 1.0))],
)
def test_schedule_rejects(phases):
    with pytest.raises(ValueError):
        Schedule(phases)


@pytest.mark.parametrize("kw", [{"weight": 0.0}, {"weight": math.nan}, {"easing": "bounce"}])
def test_phase_rejects(kw):
    with pytest.raises(ValueError):
        Phase("a", **{"weight": 1.0, **kw})


def test_bracket_and_nearest():
    keys = [0.0, 1.0, 3.0]
    assert bracket(keys, -1.0) == (0, 0, 0.0)
    assert bracket(keys, 0.5) == (0, 1, 0.5)
    assert bracket(keys, 1.0) == (1, 2, 0.0)
    assert bracket(keys, 2.0) == (1, 2, 0.5)
    assert bracket(keys, 9.0) == (2, 2, 0.0)
    assert nearest(keys, 2.1) == 2 and nearest(keys, 1.9) == 1
    assert lerp(2.0, 4.0, 0.25) == 2.5
    with pytest.raises(ValueError):
        bracket([], 1.0)


# ---- Style -----------------------------------------------------------------------


def test_series_styles_are_unique_for_40_interfaces():
    pairs = {(series_color(k), str(series_linestyle(k))) for k in range(40)}
    assert len(pairs) == 40
    assert [series_color(k) for k in range(8)] == list(CATEGORICAL)
    assert all(series_linestyle(k) == "solid" for k in range(3))
    assert all(series_linestyle(k) != "solid" for k in range(3, 40))


def test_plot_validation_uses_the_shared_palette():
    import importlib.util
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    spec = importlib.util.spec_from_file_location(
        "plot_validation", root / "validation" / "plot_validation.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["plot_validation"] = mod
    try:
        spec.loader.exec_module(mod)
        from bvidfe.viz import style

        assert mod.SURFACE == style.SURFACE and mod.TEXT == style.TEXT
        assert {c for _, c, _ in mod.STYLE.values()} <= set(CATEGORICAL)
    finally:
        sys.modules.pop("plot_validation", None)
        sys.modules.pop("validate_bvid_public", None)
