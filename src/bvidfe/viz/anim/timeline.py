"""Timelines: model outputs at keyframes, plus the schedule that maps the
normalised clip time ``t`` (0 to 1) onto them.

A scene's ``build_timeline`` runs the model once and stores only its
outputs, at the keyframe values of the swept parameter (impact energies for
the energy sweep). ``draw_frame(fig, timeline, t)`` then draws any ``t``
without calling the model: between two keyframes it interpolates the
*drawing* (an ellipse's size, a marker's position), while every number it
prints is a keyframe's model output.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, ClassVar, Dict, Tuple

import numpy as np
from numpy.typing import ArrayLike

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.figure import Figure

Easing = Callable[[float], float]


def linear(u: float) -> float:
    """No easing."""
    return u


def ease_in_out(u: float) -> float:
    """Cubic ease-in-out: starts and stops with zero speed."""
    return 4.0 * u**3 if u < 0.5 else 1.0 - (-2.0 * u + 2.0) ** 3 / 2.0


def ease_out(u: float) -> float:
    """Cubic ease-out: starts fast, settles with zero speed."""
    return 1.0 - (1.0 - u) ** 3


def smoothstep(u: float) -> float:
    """Hermite smoothstep, ``3u^2 - 2u^3``."""
    return u * u * (3.0 - 2.0 * u)


EASINGS: Dict[str, Easing] = {
    "linear": linear,
    "ease_in_out": ease_in_out,
    "ease_out": ease_out,
    "smoothstep": smoothstep,
}


def clamp01(u: float) -> float:
    """``u`` clipped to ``[0, 1]`` (NaN raises)."""
    if math.isnan(u):
        raise ValueError("t is NaN")
    return 0.0 if u < 0.0 else 1.0 if u > 1.0 else float(u)


@dataclass(frozen=True)
class Phase:
    """One segment of a clip: its name, relative length and easing."""

    name: str
    weight: float
    easing: str = "linear"

    def __post_init__(self) -> None:
        if not (math.isfinite(self.weight) and self.weight > 0):
            raise ValueError(f"Phase.weight must be > 0 (got {self.weight!r})")
        if self.easing not in EASINGS:
            raise ValueError(f"unknown easing {self.easing!r}; choose from {', '.join(EASINGS)}")


@dataclass(frozen=True)
class Schedule:
    """Consecutive phases that share the clip time ``t`` in proportion to
    their weights.

    >>> s = Schedule((Phase("hold", 1.0), Phase("sweep", 3.0)))
    >>> s.at(0.0), s.at(0.25), s.at(1.0)
    (('hold', 0.0), ('sweep', 0.0), ('sweep', 1.0))
    """

    phases: Tuple[Phase, ...]

    def __post_init__(self) -> None:
        if not self.phases:
            raise ValueError("a Schedule needs at least one phase")
        names = [p.name for p in self.phases]
        if len(set(names)) != len(names):
            raise ValueError(f"phase names must be unique (got {names})")

    def bounds(self, name: str) -> Tuple[float, float]:
        """``(t_start, t_end)`` of phase ``name``."""
        total = sum(p.weight for p in self.phases)
        start = 0.0
        for p in self.phases:
            end = start + p.weight / total
            if p.name == name:
                return start, end
            start = end
        raise KeyError(name)

    def at(self, t: float) -> Tuple[str, float]:
        """``(phase name, eased progress 0..1 within it)`` at clip time ``t``.
        A boundary belongs to the later phase, except ``t = 1``."""
        t = clamp01(t)
        total = sum(p.weight for p in self.phases)
        start = 0.0
        for k, p in enumerate(self.phases):
            end = start + p.weight / total
            last = k == len(self.phases) - 1
            if t < end or last:
                u = 0.0 if end <= start else clamp01((t - start) / (end - start))
                return p.name, EASINGS[p.easing](u)
            start = end
        raise AssertionError("unreachable")  # pragma: no cover

    def progress(self, t: float, name: str) -> float:
        """Eased progress of phase ``name`` at ``t``: 0 before it, 1 after."""
        t0, t1 = self.bounds(name)
        t = clamp01(t)
        if t <= t0:
            return 0.0
        if t >= t1:
            return 1.0
        easing = next(p.easing for p in self.phases if p.name == name)
        return EASINGS[easing]((t - t0) / (t1 - t0))


def bracket(keys: ArrayLike, x: float) -> Tuple[int, int, float]:
    """``(i, j, w)`` with ``keys[i] <= x <= keys[j]``, ``j = i + 1`` (or
    ``i == j`` at or past an end) and ``w`` the linear weight of ``keys[j]``.
    ``keys`` must be sorted ascending.

    >>> bracket([0.0, 1.0, 3.0], 2.0)
    (1, 2, 0.5)
    """
    k = np.asarray(keys, dtype=float)
    if k.size == 0:
        raise ValueError("no keyframes")
    if x <= k[0]:
        return 0, 0, 0.0
    if x >= k[-1]:
        n = int(k.size) - 1
        return n, n, 0.0
    j = int(np.searchsorted(k, x, side="right"))
    i = j - 1
    span = k[j] - k[i]
    w = 0.0 if span <= 0 else float((x - k[i]) / span)
    return i, j, w


def nearest(keys: ArrayLike, x: float) -> int:
    """Index of the keyframe nearest ``x`` (the lower one on a tie)."""
    i, j, w = bracket(keys, x)
    return j if w > 0.5 else i


def lerp(a: float, b: float, w: float) -> float:
    """``a + (b - a) w``."""
    return a + (b - a) * w


class Timeline:
    """Base of the scene timelines. A subclass is a frozen dataclass of
    model outputs with a ``schedule``; ``draw(fig, t)`` draws the frame at
    clip time ``t`` with the scene's ``draw_frame``."""

    scene: ClassVar[str] = ""
    schedule: Schedule

    def draw(self, fig: "Figure", t: float) -> None:  # pragma: no cover - abstract
        raise NotImplementedError
