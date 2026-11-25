"""1080p animations of what BVID-FE computes.

Three scenes, each in two steps: ``build_timeline`` runs the model once and
keeps its outputs; ``draw_frame(fig, timeline, t)`` draws the frame at clip
time ``t`` (0 to 1) without calling the model.

* ``energy-sweep`` (:func:`render_energy_sweep`): the delaminations, damage
  area, dent depth and knockdown as the impact energy rises.
* ``damage-stack`` (:func:`render_damage_stack`): the delaminations through
  the thickness as an exploded 3D stack, with the camera orbiting.
* ``cai-ramp`` (:func:`render_cai_ramp`): the fe3d load ramp of the damaged
  panel to its residual strength, each element coloured by its load ratio
  to local failure (two fe3d solves, then linear scaling).

Every frame carries the citation of BVID-FE and its version in the footer
(:mod:`bvidfe.citation`). MP4 output needs ffmpeg: ``pip install
"bvidfe[anim]"`` bundles one (see :func:`find_ffmpeg`).

>>> from bvidfe.viz.anim import AnimationSpec
>>> AnimationSpec().size, AnimationSpec().fps
((1920, 1080), 30)
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

from bvidfe.analysis.config import AnalysisConfig
from bvidfe.core.geometry import PanelGeometry
from bvidfe.damage.state import DamageState
from bvidfe.viz.anim.frame import new_figure, render_rgba, save_png
from bvidfe.viz.anim.scenes.cai_ramp import CaiRampTimeline, RatioKind
from bvidfe.viz.anim.scenes.cai_ramp import build_timeline as build_cai_ramp
from bvidfe.viz.anim.scenes.damage_stack import DamageStackTimeline
from bvidfe.viz.anim.scenes.damage_stack import build_timeline as build_damage_stack
from bvidfe.viz.anim.scenes.energy_sweep import EnergySweepTimeline
from bvidfe.viz.anim.scenes.energy_sweep import build_timeline as build_energy_sweep
from bvidfe.viz.anim.spec import AnimationSpec
from bvidfe.viz.anim.timeline import Timeline
from bvidfe.viz.anim.writer import (
    FFmpegFallbackWarning,
    FFmpegNotFoundError,
    ProgressCallback,
    RenderResult,
    encode,
    find_ffmpeg,
)

PathLike = Union[str, Path]

SCENES: Dict[str, str] = {
    "energy-sweep": "Damage growing with impact energy, with area, dent and knockdown curves",
    "damage-stack": "Exploded 3D stack of the delaminations, camera orbiting",
    "cai-ramp": "fe3d load ramp to failure, elements coloured by their load ratio",
}

# Encoding one 1080p frame on top of drawing it (s): the RGBA buffer piped
# to ffmpeg, x264 working in parallel.
_ENCODE_S_PER_FRAME = 0.04

__all__ = [
    "AnimationSpec",
    "CaiRampTimeline",
    "DamageStackTimeline",
    "EnergySweepTimeline",
    "FFmpegFallbackWarning",
    "FFmpegNotFoundError",
    "RenderResult",
    "SCENES",
    "Timeline",
    "build_cai_ramp",
    "build_damage_stack",
    "build_energy_sweep",
    "encode_timeline",
    "estimate_render_seconds",
    "find_ffmpeg",
    "list_animations",
    "preview_frames",
    "render_cai_ramp",
    "render_damage_stack",
    "render_energy_sweep",
]


def list_animations() -> Dict[str, str]:
    """The scene names and what each shows."""
    return dict(SCENES)


def encode_timeline(
    timeline: Timeline,
    path: PathLike,
    spec: Optional[AnimationSpec] = None,
    progress: Optional[ProgressCallback] = None,
) -> RenderResult:
    """Draw every frame of ``timeline`` and write it to ``path`` in
    ``spec``'s format (default 1080p MP4). Returns what was written."""
    spec = spec or AnimationSpec()
    fig = new_figure(spec)
    return encode(fig, timeline.draw, path, spec, progress)


def estimate_render_seconds(
    timeline: Timeline,
    spec: Optional[AnimationSpec] = None,
    frame_seconds: Optional[float] = None,
) -> float:
    """A rough wall time (s) to encode ``timeline`` with ``spec`` (default
    1080p MP4, 10 s at 30 fps): frames x (drawing time + encoding).

    ``frame_seconds`` is the drawing time of one frame, e.g. measured on
    this machine while drawing :func:`preview_frames`; by default the
    scene's own estimate (``timeline.frame_seconds_estimate()``, measured
    on one core of a Linux container), scaled with the frame area.
    """
    spec = spec or AnimationSpec()
    w, h = spec.size
    area = max(0.35, (w * h) / (1920 * 1080))
    per = timeline.frame_seconds_estimate() * area if frame_seconds is None else frame_seconds
    encode = _ENCODE_S_PER_FRAME * area if spec.format in ("mp4", "auto") else 0.0
    return float(spec.n_frames * (per + encode))


def preview_frames(
    timeline: Timeline,
    t: Sequence[float] = (0.0, 0.5, 1.0),
    spec: Optional[AnimationSpec] = None,
    out_dir: Optional[PathLike] = None,
    stem: str = "frame",
) -> List[np.ndarray]:
    """Stills of ``timeline`` at clip times ``t``, as ``(height, width, 4)``
    RGBA arrays at ``spec``'s size (default 1920 x 1080). With ``out_dir``
    each is also written there as ``<stem>_t<t>.png`` (e.g.
    ``frame_t0.50.png``)."""
    spec = spec or AnimationSpec()
    fig = new_figure(spec)
    frames = []
    for ti in t:
        timeline.draw(fig, float(ti))
        rgba = render_rgba(fig)
        frames.append(rgba)
        if out_dir is not None:
            save_png(rgba, Path(out_dir) / f"{stem}_t{float(ti):.2f}.png")
    return frames


def render_energy_sweep(
    cfg: AnalysisConfig,
    energies: Optional[Sequence[float]] = None,
    path: PathLike = "energy_sweep.mp4",
    spec: Optional[AnimationSpec] = None,
    *,
    energy_min: Optional[float] = None,
    energy_max: Optional[float] = None,
    n_energies: Optional[int] = None,
    progress: Optional[ProgressCallback] = None,
) -> EnergySweepTimeline:
    """Render the energy-sweep animation of ``cfg`` to ``path`` and return
    its timeline.

    ``cfg`` is impact-driven, tier ``empirical`` or ``semi_analytical``
    (``fe3d`` is refused: it would solve two FE models per energy).
    ``energies`` (J) are the keyframes, each one ``BvidAnalysis`` run;
    default 60 (empirical) or 40 (semi_analytical) from 0.8 x the onset
    energy (``energy_min``) to the config's energy (``energy_max``), plus
    the onset energy. See :class:`AnimationSpec` for ``spec``.
    """
    tl = build_energy_sweep(
        cfg, energies, energy_min=energy_min, energy_max=energy_max, n_energies=n_energies
    )
    encode_timeline(tl, path, spec, progress)
    return tl


def render_damage_stack(
    source: Union[AnalysisConfig, DamageState],
    path: PathLike = "damage_stack.mp4",
    spec: Optional[AnimationSpec] = None,
    *,
    layup_deg: Optional[Sequence[float]] = None,
    ply_thickness_mm: Optional[Union[float, Sequence[float]]] = None,
    panel: Optional[PanelGeometry] = None,
    intro: bool = True,
    progress: Optional[ProgressCallback] = None,
) -> DamageStackTimeline:
    """Render the exploded-stack animation to ``path`` and return its
    timeline.

    ``source`` is an :class:`AnalysisConfig` (the impact mapping runs once)
    or a :class:`DamageState` such as a C-scan, with ``layup_deg``,
    ``ply_thickness_mm`` and ``panel`` describing its laminate (optional).
    ``intro=False`` skips the assembly, for a loopable orbit.
    """
    tl = build_damage_stack(
        source,
        layup_deg=layup_deg,
        ply_thickness_mm=ply_thickness_mm,
        panel=panel,
        intro=intro,
    )
    encode_timeline(tl, path, spec, progress)
    return tl


def render_cai_ramp(
    cfg: AnalysisConfig,
    path: PathLike = "cai_ramp.mp4",
    spec: Optional[AnimationSpec] = None,
    *,
    damage: Optional[DamageState] = None,
    ratio: RatioKind = "averaged",
    progress: Optional[ProgressCallback] = None,
) -> CaiRampTimeline:
    """Render the fe3d load-ramp animation of ``cfg`` to ``path`` and return
    its timeline.

    The fe3d meshes of ``cfg`` (``cfg.mesh``, default 5 mm; ``cfg.tier`` is
    not read) are solved once, two solves (four on a biaxial config), as
    ``BvidAnalysis(cfg).run()`` with ``tier="fe3d"`` solves them; every
    frame scales those states. ``ratio="pointwise"`` colours by each
    element's lowest Gauss point instead of the disk-averaged check that
    sets first-ply failure. Raises ``ValueError`` when the 3D view would
    draw more than ``cai_ramp.MAX_FACES`` faces, and ``FE3DSizeError``
    above the fe3d size cap.
    """
    tl = build_cai_ramp(cfg, damage, ratio=ratio)
    encode_timeline(tl, path, spec, progress)
    return tl
