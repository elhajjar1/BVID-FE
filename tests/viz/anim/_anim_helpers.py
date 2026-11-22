"""Helpers of the animation tests: a small, fast configuration and a video
probe."""

from __future__ import annotations

import shutil
import subprocess

import pytest

from bvidfe import AnalysisConfig, ImpactEvent, PanelGeometry

LAYUP_8 = [0, 45, -45, 90, 90, -45, 45, 0]


def small_config(**kw) -> AnalysisConfig:
    """8-ply IM7/8552 on the 150 x 100 mm coupon, impact-driven, empirical."""
    base = dict(
        material="IM7/8552",
        layup_deg=LAYUP_8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        impact=ImpactEvent(energy_J=20.0),
    )
    base.update(kw)
    return AnalysisConfig(**base)


def probe(path) -> dict:
    """``codec``, ``width``, ``height``, ``fps``, ``frames`` and ``pix_fmt``
    of a video, from ffprobe when it is on PATH, else from imageio-ffmpeg."""
    ffprobe = shutil.which("ffprobe")
    if ffprobe:
        out = subprocess.run(
            [
                ffprobe,
                "-v",
                "error",
                "-select_streams",
                "v:0",
                "-count_frames",
                "-show_entries",
                "stream=codec_name,width,height,r_frame_rate,nb_read_frames,pix_fmt",
                "-of",
                "default=nw=1",
                str(path),
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        kv = dict(line.split("=", 1) for line in out.strip().splitlines())
        num, den = kv["r_frame_rate"].split("/")
        return {
            "codec": kv["codec_name"],
            "width": int(kv["width"]),
            "height": int(kv["height"]),
            "fps": float(num) / float(den),
            "frames": int(kv["nb_read_frames"]),
            "pix_fmt": kv["pix_fmt"],
        }
    imageio_ffmpeg = pytest.importorskip("imageio_ffmpeg")
    gen = imageio_ffmpeg.read_frames(str(path))
    meta = next(gen)
    gen.close()
    frames, _ = imageio_ffmpeg.count_frames_and_secs(str(path))
    return {
        "codec": meta["codec"],
        "width": meta["size"][0],
        "height": meta["size"][1],
        "fps": float(meta["fps"]),
        "frames": int(frames),
        "pix_fmt": meta["pix_fmt"].split("(")[0],
    }


def needs_ffmpeg() -> None:
    from bvidfe.viz.anim.writer import find_ffmpeg

    if find_ffmpeg() is None:
        pytest.skip("no ffmpeg (install the anim extra)")
