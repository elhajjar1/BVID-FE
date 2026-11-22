"""ffmpeg discovery order, output paths and the format fallbacks."""

from __future__ import annotations

import os
import sys

import matplotlib as mpl
import pytest

from bvidfe.viz.anim import writer
from bvidfe.viz.anim.spec import AnimationSpec
from bvidfe.viz.anim.writer import (
    FFmpegFallbackWarning,
    FFmpegNotFoundError,
    find_ffmpeg,
    output_path,
    require_ffmpeg,
    resolve_format,
)


def _exe(tmp_path, name):
    p = tmp_path / name
    p.write_text("#!/bin/sh\n")
    p.chmod(0o755)
    return str(p)


@pytest.fixture
def no_ffmpeg(monkeypatch):
    """No env var, default rcParams, nothing on PATH, no imageio-ffmpeg."""
    monkeypatch.delenv("BVIDFE_FFMPEG", raising=False)
    monkeypatch.setitem(mpl.rcParams, "animation.ffmpeg_path", "ffmpeg")
    monkeypatch.setattr(writer.shutil, "which", lambda name, *a, **k: None)
    monkeypatch.setitem(sys.modules, "imageio_ffmpeg", None)


class _FakeImageio:
    def __init__(self, exe):
        self.exe = exe

    def get_ffmpeg_exe(self):
        if self.exe is None:
            raise RuntimeError("no binary")
        return self.exe


def test_env_var_wins(monkeypatch, tmp_path, no_ffmpeg):
    env = _exe(tmp_path, "env-ffmpeg")
    monkeypatch.setenv("BVIDFE_FFMPEG", env)
    monkeypatch.setitem(mpl.rcParams, "animation.ffmpeg_path", _exe(tmp_path, "rc-ffmpeg"))
    monkeypatch.setattr(
        writer.shutil, "which", lambda name, *a, **k: "/opt/ffmpeg" if name == "ffmpeg" else None
    )
    assert find_ffmpeg() == env


def test_env_var_that_does_not_exist_raises(monkeypatch, tmp_path, no_ffmpeg):
    monkeypatch.setenv("BVIDFE_FFMPEG", str(tmp_path / "missing"))
    with pytest.raises(FFmpegNotFoundError, match="BVIDFE_FFMPEG"):
        find_ffmpeg()


def test_rcparams_path_second(monkeypatch, tmp_path, no_ffmpeg):
    rc = _exe(tmp_path, "rc-ffmpeg")
    monkeypatch.setitem(mpl.rcParams, "animation.ffmpeg_path", rc)
    monkeypatch.setattr(
        writer.shutil, "which", lambda name, *a, **k: "/opt/ffmpeg" if name == "ffmpeg" else None
    )
    assert find_ffmpeg() == rc


def test_rcparams_path_that_does_not_exist_is_skipped(monkeypatch, tmp_path, no_ffmpeg):
    monkeypatch.setitem(mpl.rcParams, "animation.ffmpeg_path", str(tmp_path / "gone"))
    monkeypatch.setattr(
        writer.shutil, "which", lambda name, *a, **k: "/opt/ffmpeg" if name == "ffmpeg" else None
    )
    assert find_ffmpeg() == "/opt/ffmpeg"


def test_path_third(monkeypatch, tmp_path, no_ffmpeg):
    monkeypatch.setattr(
        writer.shutil, "which", lambda name, *a, **k: "/opt/ffmpeg" if name == "ffmpeg" else None
    )
    monkeypatch.setitem(sys.modules, "imageio_ffmpeg", _FakeImageio(_exe(tmp_path, "bundled")))
    assert find_ffmpeg() == "/opt/ffmpeg"


def test_imageio_ffmpeg_last(monkeypatch, tmp_path, no_ffmpeg):
    bundled = _exe(tmp_path, "bundled")
    monkeypatch.setitem(sys.modules, "imageio_ffmpeg", _FakeImageio(bundled))
    assert find_ffmpeg() == bundled


def test_imageio_ffmpeg_without_binary(monkeypatch, no_ffmpeg):
    monkeypatch.setitem(sys.modules, "imageio_ffmpeg", _FakeImageio(None))
    assert find_ffmpeg() is None


def test_none_found_names_the_extra(no_ffmpeg):
    assert find_ffmpeg() is None
    with pytest.raises(FFmpegNotFoundError, match=r"bvidfe\[anim\]"):
        require_ffmpeg()
    assert resolve_format(AnimationSpec(format="auto")) == "png"
    assert resolve_format(AnimationSpec(format="gif")) == "gif"


def test_mp4_without_ffmpeg_raises(no_ffmpeg, tmp_path, sweep_timeline):
    from bvidfe.viz.anim import encode_timeline

    with pytest.raises(FFmpegNotFoundError):
        encode_timeline(sweep_timeline, tmp_path / "x.mp4", AnimationSpec(duration_s=1, fps=2))
    assert not (tmp_path / "x.mp4").exists()


def test_auto_without_ffmpeg_writes_pngs(no_ffmpeg, tmp_path, sweep_timeline):
    from bvidfe.viz.anim import encode_timeline

    spec = AnimationSpec(duration_s=1, fps=2, format="auto", width=640)
    with pytest.warns(FFmpegFallbackWarning, match="PNG sequence"):
        result = encode_timeline(sweep_timeline, tmp_path / "clip.mp4", spec)
    assert result.format == "png" and result.path == tmp_path / "clip_frames"
    files = sorted(os.listdir(result.path))
    assert files == ["frame_00000.png", "frame_00001.png"]
    assert result.size_bytes > 0


def test_output_paths(tmp_path):
    assert output_path(tmp_path / "a.mp4", "mp4") == tmp_path / "a.mp4"
    assert output_path(tmp_path / "a.mov", "mp4") == tmp_path / "a.mp4"
    assert output_path(tmp_path / "a.mp4", "gif") == tmp_path / "a.gif"
    assert output_path(tmp_path / "a.mp4", "png") == tmp_path / "a_frames"
    assert output_path(tmp_path / "frames", "png") == tmp_path / "frames"
