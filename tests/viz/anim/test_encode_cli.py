"""Encoding smoke tests (1 s, low fps, per scene), the GIF and PNG outputs,
the ``bvidfe-animate`` CLI, and the full 10 s renders (slow)."""

from __future__ import annotations

import subprocess
import sys

import pytest
from PIL import Image

from bvidfe import __version__
from bvidfe.viz.anim import (
    AnimationSpec,
    encode_timeline,
    list_animations,
    render_damage_stack,
    render_energy_sweep,
)
from bvidfe.viz.anim.cli import main
from tests.viz.anim._anim_helpers import needs_ffmpeg, probe, small_config

CONFIG_ARGS = [
    "--material",
    "IM7/8552",
    "--layup",
    "0,45,-45,90,90,-45,45,0",
    "--thickness",
    "0.152",
    "--panel",
    "150x100",
]


@pytest.mark.parametrize("scene", ["sweep", "stack"])
def test_mp4_smoke_encode(scene, request, tmp_path):
    needs_ffmpeg()
    tl = request.getfixturevalue(f"{scene}_timeline")
    spec = AnimationSpec(duration_s=1.0, fps=4)
    result = encode_timeline(tl, tmp_path / "clip.mp4", spec)
    assert result.path == tmp_path / "clip.mp4" and result.n_frames == 4
    info = probe(result.path)
    assert info["codec"] == "h264"
    assert (info["width"], info["height"]) == (1920, 1080)
    assert info["fps"] == pytest.approx(4.0)
    assert info["frames"] == 4
    assert info["pix_fmt"] == "yuv420p"


def test_gif_and_png_outputs(sweep_timeline, tmp_path):
    gif = encode_timeline(
        sweep_timeline, tmp_path / "prev.mp4", AnimationSpec(format="gif", duration_s=1, fps=3)
    )
    assert gif.path == tmp_path / "prev.gif"
    with Image.open(gif.path) as im:
        assert im.size == (640, 360) and im.n_frames == 3
    png = encode_timeline(
        sweep_timeline, tmp_path / "seq", AnimationSpec(format="png", duration_s=1, fps=2)
    )
    files = sorted(png.path.glob("*.png"))
    assert len(files) == 2
    with Image.open(files[0]) as im:
        assert im.size == (1920, 1080)


def test_render_functions_return_the_timeline(tmp_path):
    progress = []
    spec = AnimationSpec(format="gif", duration_s=1, fps=2, width=320)
    tl = render_energy_sweep(
        small_config(),
        energies=[0.2, 10.0],  # the onset (0.84 J) is added
        path=tmp_path / "a.gif",
        spec=spec,
        progress=lambda i, n: progress.append((i, n)),
    )
    assert len(tl.keyframes) == 3 and progress == [(1, 2), (2, 2)]
    assert (tmp_path / "a.gif").exists()
    tl2 = render_damage_stack(small_config(), path=tmp_path / "b.gif", spec=spec, intro=False)
    assert tl2.damage.delaminations and (tmp_path / "b.gif").exists()
    assert set(list_animations()) == {"energy-sweep", "damage-stack", "cai-ramp"}


# ---- CLI -------------------------------------------------------------------------------


@pytest.mark.parametrize("scene", ["energy-sweep", "damage-stack"])
def test_cli_preview(scene, tmp_path, capsys):
    extra = ["--energy-max", "15", "--n-energies", "4"] if scene == "energy-sweep" else []
    extra += ["--energy", "15"] if scene == "damage-stack" else []
    rc = main([scene, *CONFIG_ARGS, *extra, "--preview", "--out", str(tmp_path / "p.mp4")])
    assert rc == 0
    printed = capsys.readouterr().out.split()
    assert printed == [str(tmp_path / f"p_t{t}.png") for t in ("0.00", "0.50", "1.00")]
    for p in printed:
        with Image.open(p) as im:
            assert im.size == (1920, 1080)


def test_cli_damage_stack_from_cscan_encodes(tmp_path, capsys):
    from tests.viz.anim.test_scenes import SAMPLE_CSCAN

    rc = main(
        [
            "damage-stack",
            *CONFIG_ARGS,
            "--cscan",
            str(SAMPLE_CSCAN),
            "--format",
            "gif",
            "--duration",
            "1",
            "--fps",
            "2",
            "--width",
            "320",
            "--no-intro",
            "--out",
            str(tmp_path / "s.gif"),
        ]
    )
    assert rc == 0
    out = capsys.readouterr()
    assert out.out.strip() == str(tmp_path / "s.gif")
    assert "320 x 180, 2 frames at 2 fps" in out.err


@pytest.mark.parametrize(
    "argv, message",
    [
        (["energy-sweep", *CONFIG_ARGS], "--energy or --energy-max"),
        (["energy-sweep", *CONFIG_ARGS, "--energy", "10", "--tier", "fe3d"], "empirical or"),
        (
            ["energy-sweep", *CONFIG_ARGS, "--energy", "10", "--tier", "empirical,fe3d"],
            "one tier",
        ),
        (["energy-sweep", *CONFIG_ARGS, "--cscan", "examples/sample_cscan.json"], "damage-stack"),
        (["damage-stack", *CONFIG_ARGS], "--energy"),
        (["damage-stack", *CONFIG_ARGS, "--energy", "10", "--width", "1000"], "16:9"),
        (["damage-stack", "--layup", "0,90", "--energy", "10"], "missing required"),
    ],
)
def test_cli_usage_errors(argv, message, capsys, monkeypatch):
    from tests.viz.anim.test_scenes import SAMPLE_CSCAN

    monkeypatch.chdir(SAMPLE_CSCAN.parents[1])
    with pytest.raises(SystemExit) as exc:
        main(argv)
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_cli_module_version():
    out = subprocess.run(
        [sys.executable, "-m", "bvidfe.viz.anim", "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    assert out.strip() == f"bvidfe-animate {__version__}"


# ---- Full renders (slow) -----------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("scene", ["sweep", "stack"])
def test_full_10s_1080p_render(scene, tmp_path):
    needs_ffmpeg()
    cfg = small_config(layup_deg=[0, 45, -45, 90] * 4, impact=small_config().impact)
    path = tmp_path / f"{scene}.mp4"
    if scene == "sweep":
        render_energy_sweep(cfg, path=path)
    else:
        render_damage_stack(cfg, path=path)
    info = probe(path)
    assert info["codec"] == "h264" and (info["width"], info["height"]) == (1920, 1080)
    assert info["fps"] == pytest.approx(30.0) and info["frames"] == 300
    assert 0.2e6 < path.stat().st_size < 30e6
