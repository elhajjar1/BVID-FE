"""Frame sinks: ffmpeg discovery and the MP4, GIF and PNG-sequence writers.

Frames are streamed one at a time (``FFMpegWriter`` pipes each to ffmpeg's
stdin), so memory does not grow with the clip length; the GIF writer holds
its reduced-size frames until it saves.
"""

from __future__ import annotations

import os
import shutil
import time
import warnings
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Final, Iterator, Optional, Union

import matplotlib as mpl
from matplotlib.animation import FFMpegWriter, PillowWriter
from matplotlib.figure import Figure

from bvidfe.viz.anim.spec import AnimationSpec

PathLike = Union[str, "os.PathLike[str]"]
ProgressCallback = Callable[[int, int], None]

FFMPEG_ENV = "BVIDFE_FFMPEG"
_MPL_FFMPEG_KEY: Final = "animation.ffmpeg_path"
_INSTALL_HINT = (
    'install the optional extra, pip install "bvidfe[anim]" (it bundles an ffmpeg '
    "through imageio-ffmpeg), or put ffmpeg on PATH, or set BVIDFE_FFMPEG to its path"
)


GIF_COLOURS = 160  # palette size of a GIF preview


class FFmpegNotFoundError(RuntimeError):
    """No ffmpeg executable was found for an MP4 render."""


class FFmpegFallbackWarning(UserWarning):
    """``format="auto"`` found no ffmpeg and wrote a PNG sequence instead."""


def _executable(candidate: str) -> Optional[str]:
    """``candidate`` resolved to an existing executable (a path or a name
    on PATH), else ``None``."""
    found = shutil.which(candidate)
    if found:
        return found
    p = Path(candidate)
    return str(p) if p.is_file() else None


def find_ffmpeg() -> Optional[str]:
    """The ffmpeg executable to use, or ``None`` when there is none.

    In order: the ``BVIDFE_FFMPEG`` environment variable; matplotlib's
    ``rcParams["animation.ffmpeg_path"]`` when set to something other than
    the plain ``"ffmpeg"`` and existing; ``ffmpeg`` on PATH; the binary of
    ``imageio-ffmpeg`` (the ``anim`` extra).

    Raises ``FFmpegNotFoundError`` when ``BVIDFE_FFMPEG`` is set but names
    no executable (an explicit setting is not silently skipped).
    """
    env = os.environ.get(FFMPEG_ENV, "").strip()
    if env:
        found = _executable(env)
        if found is None:
            raise FFmpegNotFoundError(f"{FFMPEG_ENV}={env!r} is not an ffmpeg executable")
        return found
    rc = str(mpl.rcParams[_MPL_FFMPEG_KEY] or "").strip()
    if rc and rc != "ffmpeg":
        found = _executable(rc)
        if found is not None:
            return found
    found = shutil.which("ffmpeg")
    if found:
        return found
    try:
        import imageio_ffmpeg
    except ImportError:
        return None
    try:
        exe = imageio_ffmpeg.get_ffmpeg_exe()
    except RuntimeError:  # the package is installed without its binary
        return None
    return _executable(str(exe))


def require_ffmpeg() -> str:
    """``find_ffmpeg()``, or ``FFmpegNotFoundError`` naming the ``anim`` extra."""
    exe = find_ffmpeg()
    if exe is None:
        raise FFmpegNotFoundError(f"MP4 output needs ffmpeg, which was not found: {_INSTALL_HINT}")
    return exe


def resolve_format(spec: AnimationSpec) -> str:
    """The concrete format of ``spec``: ``"auto"`` is ``"mp4"`` when ffmpeg is
    found, else ``"png"``."""
    if spec.format != "auto":
        return spec.format
    return "mp4" if find_ffmpeg() is not None else "png"


def output_path(path: PathLike, fmt: str) -> Path:
    """Where a render of format ``fmt`` to ``path`` goes. A PNG sequence
    goes into a directory: ``path`` itself when it has no suffix, else
    ``<stem>_frames`` next to it."""
    p = Path(path)
    if fmt == "png":
        return p if p.suffix == "" else p.with_name(p.stem + "_frames")
    suffix = {"mp4": ".mp4", "gif": ".gif"}[fmt]
    return p if p.suffix.lower() == suffix else p.with_suffix(suffix)


@dataclass(frozen=True)
class RenderResult:
    """What an encode wrote: the file (or PNG directory), its format and size."""

    path: Path
    format: str
    width: int
    height: int
    fps: int
    n_frames: int
    elapsed_s: float

    @property
    def size_bytes(self) -> int:
        """File size, or the total of the PNG files of a sequence."""
        if self.path.is_dir():
            return sum(f.stat().st_size for f in self.path.glob("*.png"))
        return self.path.stat().st_size


class _SharedPaletteGifWriter(PillowWriter):
    """matplotlib's ``PillowWriter`` with one palette for the whole clip and
    transparent deltas.

    Each frame is mapped, without dithering, onto a palette taken from the
    first, middle and last frames; every pixel equal to the previous frame's
    is then written as transparent, so a frame stores only what changed.
    The decoded frames equal the palette-mapped ones exactly. A GIF of a
    mostly static frame comes out several times smaller."""

    def finish(self) -> None:
        import numpy as np
        from PIL import Image

        grabbed = self._frames  # type: ignore[attr-defined]  # set by PillowWriter.setup
        frames = [im.convert("RGB") for im in grabbed]
        n = len(frames)
        picks = [frames[0], frames[n // 2], frames[-1]] if n > 2 else frames
        w, h = frames[0].size
        sheet = Image.new("RGB", (w, h * len(picks)))
        for k, im in enumerate(picks):
            sheet.paste(im, (0, k * h))
        palette_im = sheet.quantize(colors=GIF_COLOURS, method=Image.Quantize.FASTOCTREE)
        colours = (palette_im.getpalette() or [])[: 3 * GIF_COLOURS]
        palette = colours + [0, 0, 0] * (256 - len(colours) // 3)
        clear = len(colours) // 3  # the first unused index is transparent
        out = []
        prev = None
        for im in frames:
            idx = np.asarray(im.quantize(palette=palette_im, dither=Image.Dither.NONE))
            marked = idx.copy()
            if prev is not None:
                marked[idx == prev] = clear
            prev = idx
            frame = Image.fromarray(marked, mode="P")
            frame.putpalette(palette)
            out.append(frame)
        out[0].save(
            self.outfile,
            save_all=True,
            append_images=out[1:],
            duration=int(1000 / self.fps),
            loop=0,
            optimize=False,
            transparency=clear,
            disposal=1,
        )


class _PngSequenceWriter:
    """One PNG per frame, ``frame_00000.png`` ... in a directory."""

    def __init__(self, fig: Figure, directory: Path, dpi: float) -> None:
        self.fig, self.directory, self.dpi = fig, directory, dpi
        self.count = 0
        directory.mkdir(parents=True, exist_ok=True)
        for old in directory.glob("frame_*.png"):
            old.unlink()

    def grab_frame(self) -> None:
        self.fig.savefig(
            self.directory / f"frame_{self.count:05d}.png",
            dpi=self.dpi,
            facecolor=self.fig.get_facecolor(),
        )
        self.count += 1


@contextmanager
def frame_sink(fig: Figure, path: PathLike, spec: AnimationSpec) -> Iterator[Callable[[], None]]:
    """Open a writer for ``spec`` at ``path`` and yield ``grab()``, which
    appends the figure's current state as a frame. The file is finalised
    when the block exits.

    ``format="auto"`` without ffmpeg writes a PNG sequence and warns
    (``FFmpegFallbackWarning``); ``"mp4"`` without ffmpeg raises
    ``FFmpegNotFoundError``.
    """
    fmt = resolve_format(spec)
    out = output_path(path, fmt)
    if fmt == "png":
        if spec.format == "auto":
            warnings.warn(
                f"ffmpeg not found, so writing a PNG sequence to {out} instead of an MP4; "
                f"for MP4, {_INSTALL_HINT}",
                FFmpegFallbackWarning,
                stacklevel=3,
            )
        png = _PngSequenceWriter(fig, out, spec.dpi)
        yield png.grab_frame
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    fps = int(spec.fps or 30)
    if fmt == "gif":
        gif = _SharedPaletteGifWriter(fps=fps)
        with gif.saving(fig, str(out), dpi=spec.dpi):
            yield gif.grab_frame
        return
    exe = require_ffmpeg()
    extra = ["-pix_fmt", spec.pix_fmt, "-crf", str(spec.crf), "-movflags", "+faststart"]
    rc: Dict[Any, Any] = {_MPL_FFMPEG_KEY: exe, "savefig.bbox": None}
    with mpl.rc_context(rc):
        mp4 = FFMpegWriter(fps=fps, codec=spec.codec, extra_args=extra)
        with mp4.saving(fig, str(out), dpi=spec.dpi):
            yield mp4.grab_frame


def encode(
    fig: Figure,
    draw: Callable[[Figure, float], None],
    path: PathLike,
    spec: AnimationSpec,
    progress: Optional[ProgressCallback] = None,
) -> RenderResult:
    """Draw and write every frame of ``spec``: ``draw(fig, t)`` for each
    frame time ``t`` in ``spec.frame_times()``, then grab. ``progress(i,
    n)`` is called after each frame."""
    t0 = time.perf_counter()
    times = spec.frame_times()
    fmt = resolve_format(spec)
    with frame_sink(fig, path, spec) as grab:
        for i, t in enumerate(times):
            draw(fig, float(t))
            grab()
            if progress is not None:
                progress(i + 1, len(times))
    width, height = spec.size
    assert spec.fps is not None
    return RenderResult(
        path=output_path(path, fmt),
        format=fmt,
        width=width,
        height=height,
        fps=int(spec.fps),
        n_frames=len(times),
        elapsed_s=time.perf_counter() - t0,
    )
