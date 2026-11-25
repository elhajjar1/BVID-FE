"""A1: impact damage growing with the impact energy.

Left, the plan view of the panel from the impact face, with each
interface's delamination ellipse (the geometry of
:func:`bvidfe.viz.plots_2d.plot_damage_map`) and the impact window when
``impact_panel`` is set. Right, the projected damage area, the dent depth
(with the BVID detectability threshold) and the knockdown against energy,
each with a marker that moves with the frame's energy.

``build_timeline`` runs :class:`~bvidfe.analysis.bvid.BvidAnalysis` once per
keyframe energy; ``draw_frame`` interpolates the ellipses between two
keyframes (an interface that first delaminates fades in), while every
number it prints is the nearest keyframe's model output.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import ClassVar, List, Optional, Sequence, Tuple

import matplotlib.patches as mpatches
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from bvidfe.analysis.bvid import BvidAnalysis
from bvidfe.analysis.config import AnalysisConfig
from bvidfe.damage.state import DelaminationEllipse
from bvidfe.report.ncr import BVID_DENT_THRESHOLD_MM
from bvidfe.viz.anim import frame as fr
from bvidfe.viz.anim.scenes._common import damage_region, impact_point, onset_energy_J
from bvidfe.viz.anim.timeline import (
    Phase,
    Schedule,
    Timeline,
    bracket,
    lerp,
    nearest,
    smoothstep,
)
from bvidfe.viz.style import (
    BAND,
    BASELINE,
    CATEGORICAL,
    MUTED,
    STATUS,
    SURFACE,
    TEXT,
    TEXT_2,
    series_color,
    series_linestyle,
)

SCENE = "energy-sweep"
SUPPORTED_TIERS = ("empirical", "semi_analytical")
DEFAULT_N_ENERGIES = {"empirical": 60, "semi_analytical": 40}
ONSET_FRACTION = 0.8  # the default sweep starts at 0.8 x the onset energy

DEFAULT_SCHEDULE = Schedule(
    (
        Phase("hold_start", 0.06),
        Phase("sweep", 0.86, "ease_in_out"),
        Phase("hold_end", 0.08),
    )
)


@dataclass(frozen=True, eq=False)
class EnergyKeyframe:
    """The model outputs at one keyframe energy (a ``BvidAnalysis`` run)."""

    energy_J: float
    dpa_mm2: float
    dent_mm: float
    knockdown: float
    residual_MPa: float
    pristine_MPa: float
    fiber_break_radius_mm: float
    delaminations: Tuple[DelaminationEllipse, ...]
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True, eq=False)
class EnergySweepTimeline(Timeline):
    """Keyframes of an energy sweep and the clip schedule.

    ``energies`` are the keyframe energies (ascending); the onset energy is
    one of them when it lies inside the sweep. The clip holds at the first
    energy, sweeps to the last (eased) and holds there.
    """

    scene: ClassVar[str] = SCENE
    frame_seconds: ClassVar[float] = 0.22  # measured at 1080p, one core

    config: AnalysisConfig
    keyframes: Tuple[EnergyKeyframe, ...]
    onset_J: float
    schedule: Schedule = DEFAULT_SCHEDULE
    dent_threshold_mm: float = BVID_DENT_THRESHOLD_MM

    @property
    def energies(self) -> np.ndarray:
        return np.array([k.energy_J for k in self.keyframes])

    @property
    def tier(self) -> str:
        return self.config.tier

    def energy_at(self, t: float) -> float:
        """The drawing energy at clip time ``t``."""
        e = self.energies
        return lerp(float(e[0]), float(e[-1]), self.schedule.progress(t, "sweep"))

    def series(self, name: str) -> np.ndarray:
        """A keyframe quantity (``"dpa_mm2"``, ``"dent_mm"``, ``"knockdown"``...)
        as an array over the keyframes."""
        return np.array([getattr(k, name) for k in self.keyframes], dtype=float)

    def draw(self, fig: Figure, t: float) -> None:
        draw_frame(fig, self, t)


# ---- Model -----------------------------------------------------------------


def check_config(cfg: AnalysisConfig) -> None:
    """``ValueError`` unless ``cfg`` is an impact-driven config of a
    supported tier."""
    if not isinstance(cfg, AnalysisConfig):
        raise ValueError(f"expected an AnalysisConfig (got {type(cfg).__name__})")
    if cfg.impact is None:
        raise ValueError(
            "the energy sweep needs an impact-driven config (AnalysisConfig.impact); "
            "a C-scan damage state has no energy to sweep"
        )
    if cfg.tier not in SUPPORTED_TIERS:
        raise ValueError(
            f"the energy sweep runs one analysis per keyframe energy, which tier "
            f"{cfg.tier!r} is too slow for (seconds to minutes per energy); use "
            f"{' or '.join(SUPPORTED_TIERS)}"
        )


def energy_grid(
    cfg: AnalysisConfig,
    energy_min: Optional[float] = None,
    energy_max: Optional[float] = None,
    n_energies: Optional[int] = None,
) -> np.ndarray:
    """The default keyframe energies: ``n_energies`` evenly spaced from
    ``energy_min`` (default 0.8 x the onset energy) to ``energy_max``
    (default the config's impact energy, or 4 x the onset energy when that
    is at or below the onset), plus the onset energy when inside. Default
    ``n_energies``: 60 (empirical), 40 (semi_analytical)."""
    check_config(cfg)
    assert cfg.impact is not None
    onset = onset_energy_J(cfg)
    lo = ONSET_FRACTION * onset if energy_min is None else float(energy_min)
    if energy_max is None:
        hi = cfg.impact.energy_J if cfg.impact.energy_J > 1.05 * onset else 4.0 * onset
    else:
        hi = float(energy_max)
    n = DEFAULT_N_ENERGIES[cfg.tier] if n_energies is None else n_energies
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError(f"n_energies must be an integer >= 2 (got {n!r})")
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo < 0 or hi <= lo:
        raise ValueError(
            f"the energy range must satisfy 0 <= energy_min < energy_max (got {lo:g} to {hi:g} J)"
        )
    return _with_onset(np.linspace(lo, hi, int(n)), onset)


def _with_onset(energies: np.ndarray, onset: float) -> np.ndarray:
    e = np.unique(np.asarray(energies, dtype=float))
    if e[0] < onset < e[-1] and not np.any(np.isclose(e, onset, rtol=0, atol=1e-9)):
        e = np.sort(np.append(e, onset))
    return e


def _keyframe(cfg: AnalysisConfig, energy: float) -> EnergyKeyframe:
    assert cfg.impact is not None
    r = BvidAnalysis(replace(cfg, impact=replace(cfg.impact, energy_J=float(energy)))).run()
    return EnergyKeyframe(
        energy_J=float(energy),
        dpa_mm2=float(r.dpa_mm2),
        dent_mm=float(r.damage.dent_depth_mm),
        knockdown=float(r.knockdown),
        residual_MPa=float(r.residual_strength_MPa),
        pristine_MPa=float(r.pristine_strength_MPa),
        fiber_break_radius_mm=float(r.damage.fiber_break_radius_mm),
        delaminations=tuple(r.damage.delaminations),
        warnings=tuple(r.warnings),
    )


def build_timeline(
    cfg: AnalysisConfig,
    energies: Optional[Sequence[float]] = None,
    *,
    energy_min: Optional[float] = None,
    energy_max: Optional[float] = None,
    n_energies: Optional[int] = None,
    schedule: Schedule = DEFAULT_SCHEDULE,
) -> EnergySweepTimeline:
    """Run ``cfg`` at each keyframe energy and return the timeline.

    ``energies`` (J, at least two, ``>= 0``) replaces the default grid of
    :func:`energy_grid`; the onset energy is added when it lies inside.
    ``cfg`` must be impact-driven with tier ``empirical`` or
    ``semi_analytical`` (``ValueError`` otherwise; fe3d would solve two FE
    models per energy). Impact-model warnings are shown once per sweep.
    """
    from bvidfe.sweep.parametric_sweep import _dedupe_impact_warnings

    check_config(cfg)
    onset = onset_energy_J(cfg)
    if energies is None:
        grid = energy_grid(cfg, energy_min, energy_max, n_energies)
    else:
        arr = np.asarray(energies, dtype=float).ravel()
        if arr.size < 2 or not np.all(np.isfinite(arr)) or np.any(arr < 0):
            raise ValueError("energies must be at least two finite values >= 0 (J)")
        grid = _with_onset(arr, onset)
        if grid.size < 2:
            raise ValueError("energies must hold at least two distinct values")
    with _dedupe_impact_warnings():
        keyframes = tuple(_keyframe(cfg, float(e)) for e in grid)
    return EnergySweepTimeline(config=cfg, keyframes=keyframes, onset_J=onset, schedule=schedule)


# ---- Drawing -----------------------------------------------------------------

# Layout, px on the 1920 x 1080 canvas.
_PLAN = (110, 214, 1030, 730)  # left, top, width, height
_LEGEND_X = 1188
_CURVES_X, _CURVES_W = 1478, 394
_BLOCK_TOPS = (176, 452, 728)
_BLOCK_AXES_DY, _BLOCK_AXES_H = 62, 160


@dataclass(frozen=True)
class _DrawnEllipse:
    interface_index: int
    centroid_mm: Tuple[float, float]
    major_mm: float
    minor_mm: float
    orientation_deg: float
    fade: float  # 0 (not yet delaminated) to 1


def interpolate_ellipses(
    a: Sequence[DelaminationEllipse], b: Sequence[DelaminationEllipse], w: float
) -> List[_DrawnEllipse]:
    """The ellipses drawn a fraction ``w`` of the way from keyframe ``a`` to
    ``b``: semi-axes lerped per interface (from zero for an interface that
    is new in ``b``, which also fades in)."""
    by_a = {d.interface_index: d for d in a}
    by_b = {d.interface_index: d for d in b}
    out = []
    for idx in sorted(set(by_a) | set(by_b)):
        da, db = by_a.get(idx), by_b.get(idx)
        ref = db if db is not None else da
        assert ref is not None
        if da is not None and db is not None:
            cx = lerp(da.centroid_mm[0], db.centroid_mm[0], w)
            cy = lerp(da.centroid_mm[1], db.centroid_mm[1], w)
            centroid = (cx, cy)
        else:
            centroid = ref.centroid_mm
        major = lerp(da.major_mm if da else 0.0, db.major_mm if db else 0.0, w)
        minor = lerp(da.minor_mm if da else 0.0, db.minor_mm if db else 0.0, w)
        fade = smoothstep(lerp(1.0 if da else 0.0, 1.0 if db else 0.0, w))
        if major <= 0 or minor <= 0 or fade <= 0:
            continue
        out.append(_DrawnEllipse(idx, centroid, major, minor, ref.orientation_deg, fade))
    return out


def _interfaces(tl: EnergySweepTimeline) -> List[int]:
    found = {d.interface_index for k in tl.keyframes for d in k.delaminations}
    return sorted(found)


def _rgba(colour: str, alpha: float) -> Tuple[float, float, float, float]:
    from matplotlib.colors import to_rgba

    return to_rgba(colour, alpha)


def plan_region(tl: EnergySweepTimeline) -> Tuple[Tuple[float, float, float, float], bool]:
    """The plan view's region: the whole panel, or a zoom around the impact
    point when even the last keyframe's damage is small on the panel."""
    _, _, width, height = _PLAN
    return damage_region(
        tl.config.panel,
        tl.keyframes[-1].delaminations,
        impact_point(tl.config),
        width / height,
    )


def _plan_view(fig: Figure, tl: EnergySweepTimeline, energy: float) -> None:
    cfg = tl.config
    panel = cfg.panel
    left, top, width, height = _PLAN
    title = "Delaminations, plan view from the impact face"
    if plan_region(tl)[1]:
        title += ", zoomed to the damage"
    fr.text_px(fig, left, 186, title, fr.PANEL_PX, fr.TEXT, weight="bold", ha="left")
    ax = fr.px_axes(fig, left, top, width, height)
    (rx0, ry0, rx1, ry1), zoomed = plan_region(tl)
    pad = 0.0 if zoomed else 0.03 * max(panel.Lx_mm, panel.Ly_mm)
    ax.set_xlim(rx0 - pad, rx1 + pad)
    ax.set_ylim(ry0 - pad, ry1 + pad)
    ax.set_aspect("equal", adjustable="box", anchor="NW")
    fr.style_axes(ax, grid=False)
    ax.set_xlabel("x (mm), load direction", fontsize=fr.pt(fr.AXIS_PX), color=TEXT_2)
    ax.set_ylabel("y (mm)", fontsize=fr.pt(fr.AXIS_PX), color=TEXT_2)
    ax.add_patch(
        mpatches.Rectangle(
            (0, 0),
            panel.Lx_mm,
            panel.Ly_mm,
            facecolor=BAND,
            edgecolor=TEXT_2,
            linewidth=fr.lw(2.0),
            zorder=1,
        )
    )
    if cfg.impact_panel is not None:
        q = cfg.impact_panel
        x0 = (panel.Lx_mm - q.Lx_mm) / 2
        y0 = (panel.Ly_mm - q.Ly_mm) / 2
        ax.add_patch(
            mpatches.Rectangle(
                (x0, y0),
                q.Lx_mm,
                q.Ly_mm,
                fill=False,
                edgecolor=TEXT_2,
                linewidth=fr.lw(1.8),
                linestyle=(0, (6, 4)),
                zorder=2,
            )
        )
        ax.text(
            x0 + 1.0,
            y0 + q.Ly_mm - 1.0,
            f"impact window {q.Lx_mm:g} × {q.Ly_mm:g} mm",
            fontsize=fr.pt(fr.CAPTION_PX),
            color=TEXT_2,
            ha="left",
            va="top",
            zorder=6,
        )

    i, j, w = bracket(tl.energies, energy)
    drawn = interpolate_ellipses(tl.keyframes[i].delaminations, tl.keyframes[j].delaminations, w)
    # Back face first, so the interfaces nearer the impact face draw on top.
    for d in sorted(drawn, key=lambda e: -e.interface_index):
        colour = series_color(d.interface_index)
        ax.add_patch(
            mpatches.Ellipse(
                xy=d.centroid_mm,
                width=2 * d.major_mm,
                height=2 * d.minor_mm,
                angle=d.orientation_deg,
                facecolor=_rgba(colour, 0.10 * d.fade),
                edgecolor=_rgba(colour, d.fade),
                linewidth=fr.lw(2.6),
                linestyle=series_linestyle(d.interface_index),
                zorder=3,
            )
        )
    r_fb = lerp(tl.keyframes[i].fiber_break_radius_mm, tl.keyframes[j].fiber_break_radius_mm, w)
    loc = impact_point(cfg)
    if r_fb > 0:
        ax.add_patch(
            mpatches.Circle(
                loc,
                r_fb,
                facecolor=_rgba(TEXT, 0.25),
                edgecolor=TEXT,
                hatch="////",
                linewidth=fr.lw(1.5),
                zorder=4,
            )
        )
    ax.plot(
        *loc,
        marker="+",
        markersize=fr.pt(26),
        markeredgewidth=fr.lw(3),
        color=TEXT,
        zorder=5,
    )
    if energy < tl.onset_J:
        ax.text(
            loc[0],
            loc[1] - 0.08 * (ry1 - ry0),
            f"below damage onset ({tl.onset_J:.1f} J)",
            fontsize=fr.pt(28),
            color=TEXT_2,
            ha="center",
            va="top",
            zorder=6,
            bbox={"facecolor": BAND, "edgecolor": "none", "pad": 6},
        )


def _legend(fig: Figure, tl: EnergySweepTimeline, energy: float) -> None:
    """Interface key: colour and line style, the two ply angles; muted until
    the interface has delaminated."""
    layup = tl.config.layup_deg
    x = _LEGEND_X
    fr.text_px(fig, x, 232, "Interface", fr.AXIS_PX, TEXT, weight="bold", ha="left")
    fr.text_px(fig, x, 258, "plies above / below", fr.CAPTION_PX, TEXT_2, ha="left")
    i, j, w = bracket(tl.energies, energy)
    on = {
        d.interface_index
        for d in interpolate_ellipses(
            tl.keyframes[i].delaminations, tl.keyframes[j].delaminations, w
        )
    }
    ifaces = _interfaces(tl)
    step = min(32.0, (930 - 296) / max(1, len(ifaces)))
    ax = fr.px_axes(fig, x, 280, 260, 930 - 280)
    ax.set_xlim(0, 260)
    ax.set_ylim(930 - 280, 0)
    ax.set_axis_off()
    for k, idx in enumerate(ifaces):
        y = 16 + k * step
        active = idx in on
        ax.plot(
            [0, 44],
            [y, y],
            color=series_color(idx) if active else BASELINE,
            linewidth=fr.lw(3.5),
            linestyle=series_linestyle(idx),
            solid_capstyle="butt",
            dash_capstyle="butt",
        )
        a, b = layup[idx], layup[idx + 1]
        ax.text(
            58,
            y,
            f"{idx:>2}   {a:g}° / {b:g}°",
            fontsize=fr.pt(fr.CAPTION_PX),
            color=TEXT_2 if active else MUTED,
            va="center",
            ha="left",
        )
    if not ifaces:
        ax.text(0, 16, "no delamination\nin this sweep", fontsize=fr.pt(16), color=MUTED, va="top")


def _block(
    fig: Figure,
    top: float,
    title: str,
    readout: str,
    caption: str,
    tl: EnergySweepTimeline,
    values: np.ndarray,
    energy: float,
    ylim: Tuple[float, float],
    xlabel: bool,
    caption_icon: str = "",
) -> Axes:
    """One mini-plot against energy: a title row with the keyframe's value,
    a caption row, then the curve (past solid, future faint) and a marker
    at the frame's energy."""
    x0 = _CURVES_X - 8
    fr.text_px(fig, x0, top + 14, title, fr.PANEL_PX, TEXT, weight="bold", ha="left")
    fr.text_px(fig, _CURVES_X + _CURVES_W, top + 14, readout, fr.PANEL_PX, TEXT, ha="right")
    if caption_icon:
        fr.text_px(fig, x0, top + 44, "\u25b2", 16, caption_icon, ha="left")
        fr.text_px(fig, x0 + 22, top + 44, caption, fr.CAPTION_PX, TEXT, weight="bold", ha="left")
    else:
        fr.text_px(fig, x0, top + 44, caption, fr.CAPTION_PX, TEXT_2, ha="left")
    ax = fr.px_axes(fig, _CURVES_X + 62, top + _BLOCK_AXES_DY, _CURVES_W - 62, _BLOCK_AXES_H)
    fr.style_axes(ax)
    e = tl.energies
    colour = CATEGORICAL[0]
    ax.plot(e, values, color=colour, alpha=0.25, linewidth=fr.lw(2.5), zorder=2)
    past = e <= energy
    v_now = float(np.interp(energy, e, values))
    ax.plot(
        np.append(e[past], energy),
        np.append(values[past], v_now),
        color=colour,
        linewidth=fr.lw(3.5),
        zorder=3,
    )
    ax.plot(
        [energy],
        [v_now],
        marker="o",
        markersize=fr.pt(16),
        markerfacecolor=colour,
        markeredgecolor=SURFACE,
        markeredgewidth=fr.lw(2.5),
        zorder=5,
        clip_on=False,
    )
    if e[0] < tl.onset_J < e[-1]:
        ax.axvline(tl.onset_J, color=BASELINE, linewidth=fr.lw(1.8), zorder=1)
    ax.set_xlim(float(e[0]), float(e[-1]))
    ax.set_ylim(*ylim)
    ax.locator_params(axis="y", nbins=4)
    ax.locator_params(axis="x", nbins=6)
    if xlabel:
        ax.set_xlabel("Impact energy (J)", fontsize=fr.pt(fr.AXIS_PX), color=TEXT_2)
    return ax


def _curves(fig: Figure, tl: EnergySweepTimeline, energy: float, kf: EnergyKeyframe) -> None:
    dpa = tl.series("dpa_mm2")
    dent = tl.series("dent_mm")
    kd = tl.series("knockdown")
    thr = tl.dent_threshold_mm
    e = tl.energies
    ax = _block(
        fig,
        _BLOCK_TOPS[0],
        "Damage area",
        f"{kf.dpa_mm2:,.0f} mm\u00b2",
        "projected area of the delaminations (mm\u00b2)",
        tl,
        dpa,
        energy,
        (0.0, max(1.0, float(dpa.max())) * 1.1),
        xlabel=False,
    )
    if e[0] < tl.onset_J < e[-1]:
        ax.text(
            tl.onset_J,
            ax.get_ylim()[1],
            f" onset {tl.onset_J:.1f} J",
            fontsize=fr.pt(fr.CAPTION_PX),
            color=TEXT_2,
            ha="left",
            va="top",
        )
    visible = kf.dent_mm >= thr
    ax = _block(
        fig,
        _BLOCK_TOPS[1],
        "Dent depth",
        f"{kf.dent_mm:.2f} mm",
        (
            f"now visible: dent \u2265 {thr:.2f} mm"
            if visible
            else f"below the {thr:.2f} mm detectability threshold"
        ),
        tl,
        dent,
        energy,
        (0.0, max(float(dent.max()) * 1.1, thr * 1.35)),
        xlabel=False,
        caption_icon=STATUS["serious"] if visible else "",
    )
    ax.axhline(thr, color=TEXT_2, linewidth=fr.lw(1.8), linestyle=(0, (5, 3)), zorder=1)
    ax.text(
        float(e[0]) + 0.02 * float(e[-1] - e[0]),
        thr,
        f"threshold {thr:.2f} mm",
        fontsize=fr.pt(fr.CAPTION_PX),
        color=TEXT_2,
        ha="left",
        va="bottom",
    )
    load = "CAI" if tl.config.loading == "compression" else "TAI"
    _block(
        fig,
        _BLOCK_TOPS[2],
        f"Knockdown, {load}",
        f"{kf.knockdown:.3f}",
        "residual / pristine strength",
        tl,
        kd,
        energy,
        (0.0, 1.05),
        xlabel=True,
    )


def draw_frame(fig: Figure, timeline: EnergySweepTimeline, t: float) -> None:
    """Draw the frame at clip time ``t`` (0 to 1) into ``fig`` (cleared
    first). Pure: no model call."""
    tl = timeline
    cfg = tl.config
    with fr.style_context():
        fig.clear()
        fig.set_facecolor(SURFACE)
        energy = tl.energy_at(t)
        kf = tl.keyframes[nearest(tl.energies, energy)]
        assert cfg.impact is not None
        subtitle = (
            f"{fr.material_name(cfg)}  {fr.layup_label(cfg.layup_deg)}  ·  "
            f"{fr.impactor_label(cfg.impact)}  ·  "
            f"{cfg.panel.Lx_mm:g} × {cfg.panel.Ly_mm:g} mm panel"
        )
        load = "CAI" if cfg.loading == "compression" else "TAI"
        fr.draw_title_band(
            fig,
            "Impact damage as the impact energy rises",
            subtitle,
            readout=f"{kf.energy_J:.1f} J",
            readout_sub=f"impact energy  \u00b7  {tl.tier} tier, {load}",
        )
        fr.draw_caveats(fig, fr.honesty_notes(cfg))
        _plan_view(fig, tl, energy)
        _legend(fig, tl, energy)
        _curves(fig, tl, energy, kf)
        fr.draw_footer(
            fig,
            fr.config_summary(cfg, tl.tier)
            + f"  ·  {len(tl.keyframes)} model runs; frames between them interpolated",
        )
