"""AnimationSpec: the output settings of an animation (size, rate, length, codec)."""

from __future__ import annotations

import math
import numbers
from dataclasses import dataclass, field
from typing import Literal, Optional, Tuple

import numpy as np

AnimationFormat = Literal["mp4", "gif", "png", "auto"]
_FORMATS: Tuple[str, ...] = ("mp4", "gif", "png", "auto")

# Every frame is laid out on a 19.2 x 10.8 in canvas (1920 x 1080 px at
# 100 dpi); another 16:9 size renders the same layout at another dpi.
DESIGN_WIDTH_PX = 1920
DESIGN_HEIGHT_PX = 1080
DESIGN_DPI = 100.0
FIGSIZE_IN = (DESIGN_WIDTH_PX / DESIGN_DPI, DESIGN_HEIGHT_PX / DESIGN_DPI)

# Defaults by format: full HD video and frames; a reduced-size GIF preview.
_DEFAULT_SIZE = {"gif": (640, 360)}
_DEFAULT_FPS = {"gif": 12}
_MAX_FRAMES = 36_000  # 20 min at 30 fps


def _is_int(v: object) -> bool:
    return isinstance(v, numbers.Integral) and not isinstance(v, bool)


@dataclass(frozen=True)
class AnimationSpec:
    """Output settings of an animation.

    Parameters
    ----------
    width, height : int, optional
        Frame size in pixels, 16:9 and even (H.264 with ``yuv420p`` needs
        even sizes). Default 1920 x 1080, or 640 x 360 for ``format="gif"``
        (a preview). Only 1080p is tested; other 16:9 sizes render the same
        layout scaled.
    fps : int, optional
        Frames per second, 1 to 120. Default 30, or 12 for a GIF.
    duration_s : float
        Clip length in seconds, ``> 0`` (default 10). The clip has
        ``round(duration_s * fps)`` frames, at least 2, the first at
        ``t = 0`` and the last at ``t = 1``.
    codec : str
        ffmpeg video encoder for MP4 (default ``"libx264"``, H.264).
    crf : int
        Constant rate factor, 0 to 51 (default 18, visually lossless).
    pix_fmt : str
        ffmpeg output pixel format (default ``"yuv420p"``, which every
        player decodes).
    format : {"mp4", "gif", "png", "auto"}
        ``"mp4"`` (default) needs ffmpeg; ``"gif"`` uses Pillow; ``"png"``
        writes one PNG per frame into a directory; ``"auto"`` is MP4 when
        ffmpeg is found and a PNG sequence otherwise.

    Raises ``ValueError`` for a value out of range or of the wrong type.

    >>> AnimationSpec().n_frames, AnimationSpec(format="gif").width
    (300, 640)
    """

    width: Optional[int] = None
    height: Optional[int] = None
    fps: Optional[int] = None
    duration_s: float = 10.0
    codec: str = "libx264"
    crf: int = 18
    pix_fmt: str = "yuv420p"
    format: AnimationFormat = "mp4"
    # Filled in __post_init__; not an argument.
    dpi: float = field(init=False, repr=False, compare=False, default=DESIGN_DPI)

    def __post_init__(self) -> None:
        if self.format not in _FORMATS:
            raise ValueError(
                f"AnimationSpec.format must be one of {', '.join(_FORMATS)} (got {self.format!r})"
            )
        w0, h0 = _DEFAULT_SIZE.get(self.format, (DESIGN_WIDTH_PX, DESIGN_HEIGHT_PX))
        if self.width is None and self.height is None:
            object.__setattr__(self, "width", w0)
            object.__setattr__(self, "height", h0)
        elif self.width is None or self.height is None:
            # One side given: the other follows from 16:9.
            if self.width is None:
                self._check_px("height", self.height)
                assert self.height is not None
                object.__setattr__(self, "width", self.height * 16 // 9)
            else:
                self._check_px("width", self.width)
                object.__setattr__(self, "height", self.width * 9 // 16)
        if self.fps is None:
            object.__setattr__(self, "fps", _DEFAULT_FPS.get(self.format, 30))
        w, h = self.width, self.height
        self._check_px("width", w)
        self._check_px("height", h)
        assert w is not None and h is not None
        if w * 9 != h * 16:
            raise ValueError(
                f"AnimationSpec frames must be 16:9 (got {w} x {h}); e.g. 1920 x 1080, "
                "1280 x 720 or 640 x 360"
            )
        if w % 2 or h % 2:
            raise ValueError(f"AnimationSpec width and height must be even (got {w} x {h})")
        if not _is_int(self.fps) or not 1 <= int(self.fps) <= 120:  # type: ignore[arg-type]
            raise ValueError(
                f"AnimationSpec.fps must be an integer from 1 to 120 (got {self.fps!r})"
            )
        if (
            isinstance(self.duration_s, bool)
            or not isinstance(self.duration_s, numbers.Real)
            or not math.isfinite(self.duration_s)
            or self.duration_s <= 0
        ):
            raise ValueError(
                f"AnimationSpec.duration_s must be a finite number > 0 (got {self.duration_s!r})"
            )
        if self.n_frames > _MAX_FRAMES:
            raise ValueError(
                f"AnimationSpec gives {self.n_frames} frames, more than {_MAX_FRAMES}; "
                "shorten duration_s or lower fps"
            )
        if not _is_int(self.crf) or not 0 <= self.crf <= 51:
            raise ValueError(
                f"AnimationSpec.crf must be an integer from 0 to 51 (got {self.crf!r})"
            )
        for name in ("codec", "pix_fmt"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or value.startswith("-"):
                raise ValueError(f"AnimationSpec.{name} must be an ffmpeg name (got {value!r})")
        # A dpi a hair above w / 19.2 in, so the truncation to whole pixels
        # in the Agg canvas gives exactly w x h.
        object.__setattr__(self, "dpi", w / FIGSIZE_IN[0] * (1.0 + 1e-12))

    @staticmethod
    def _check_px(name: str, v: object) -> None:
        if not _is_int(v) or not 16 <= int(v) <= 7680:  # type: ignore[call-overload]
            raise ValueError(f"AnimationSpec.{name} must be an integer from 16 to 7680 (got {v!r})")

    @property
    def n_frames(self) -> int:
        """Number of frames: ``round(duration_s * fps)``, at least 2."""
        assert self.fps is not None
        return max(2, int(round(self.duration_s * self.fps)))

    def frame_times(self) -> np.ndarray:
        """The normalised time ``t`` of each frame, ``0`` to ``1`` inclusive."""
        return np.linspace(0.0, 1.0, self.n_frames)

    @property
    def size(self) -> Tuple[int, int]:
        """``(width, height)`` in pixels."""
        assert self.width is not None and self.height is not None
        return int(self.width), int(self.height)
