"""A2: the delaminations through the thickness, as an exploded 3D stack.

Every ply is a sheet at its own height, the impact face (ply 0) on top and
spread apart so the interfaces between them show; each interface's
delamination ellipse is filled in its interface colour (the colours of the
energy sweep and :func:`bvidfe.viz.plots_2d.plot_damage_map`'s geometry),
each ply carries short lines in its fibre direction, the impactor sits
above the dent on the impact face, and the camera orbits the stack. An
optional intro drops the plies in one by one.

The damage comes from a config (the impact mapping, or the config's own
damage state) or directly from a :class:`~bvidfe.damage.state.DamageState`
such as a C-scan. ``draw_frame`` is pure: the timeline holds the damage.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar, Dict, List, Optional, Sequence, Tuple, Union, cast

import numpy as np
from matplotlib.colors import LightSource, to_rgba
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from bvidfe.analysis.bvid import BvidAnalysis, _resolve_material
from bvidfe.analysis.config import AnalysisConfig
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate, _normalise_ply_thicknesses
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.viz.anim import frame as fr
from bvidfe.viz.anim.scenes._common import damage_region, impact_point, onset_energy_J
from bvidfe.viz.anim.timeline import Phase, Schedule, Timeline, clamp01, ease_out, smoothstep
from bvidfe.viz.style import (
    BAND,
    SURFACE,
    TEXT,
    TEXT_2,
    series_color,
    series_linestyle,
)

SCENE = "damage-stack"

# Fraction of the clip the assembly intro takes, then the 360 deg orbit.
INTRO_SCHEDULE = Schedule((Phase("assemble", 0.22), Phase("orbit", 0.78, "ease_in_out")))
ORBIT_SCHEDULE = Schedule((Phase("orbit", 1.0, "linear"),))

ELEV_DEG = 28.0
AZIM0_DEG = -58.0
STACK_HEIGHT_FRACTION = 0.55  # drawn stack height / longest panel side
DEFAULT_PANEL = PanelGeometry(150.0, 100.0)
IMPACTOR_GREY = "#9a9892"  # a neutral steel, between MUTED and BASELINE
DENT_GREY = "#b9b6ad"


@dataclass(frozen=True, eq=False)
class DamageStackTimeline(Timeline):
    """The damage state of one impact (or inspection) and the stack it
    sits in; the clip schedule assembles the stack (optional) and orbits.

    ``layup_deg`` is ``None`` when only a damage state was given: the plies
    are then drawn without fibre directions.
    """

    scene: ClassVar[str] = SCENE
    frame_seconds: ClassVar[float] = 0.20  # measured at 1080p, one core

    damage: DamageState
    n_plies: int
    panel: PanelGeometry
    layup_deg: Optional[Tuple[float, ...]] = None
    ply_thicknesses_mm: Optional[Tuple[float, ...]] = None
    impact: Optional[ImpactEvent] = None
    impact_point_mm: Optional[Tuple[float, float]] = None
    onset_J: Optional[float] = None
    config: Optional[AnalysisConfig] = None
    source: str = "impact model"
    schedule: Schedule = INTRO_SCHEDULE

    @property
    def intro(self) -> bool:
        return any(p.name == "assemble" for p in self.schedule.phases)

    def draw(self, fig: Figure, t: float) -> None:
        draw_frame(fig, self, t)


# ---- Model -------------------------------------------------------------------


def _infer_panel(damage: DamageState) -> PanelGeometry:
    """The 150 x 100 mm coupon, enlarged when the damage reaches past it."""
    lx, ly = DEFAULT_PANEL.Lx_mm, DEFAULT_PANEL.Ly_mm
    for d in damage.delaminations:
        ex, ey = d.half_extents_mm()
        lx = max(lx, 1.1 * (d.centroid_mm[0] + ex))
        ly = max(ly, 1.1 * (d.centroid_mm[1] + ey))
    return PanelGeometry(lx, ly)


def build_timeline(
    source: Union[AnalysisConfig, DamageState],
    *,
    layup_deg: Optional[Sequence[float]] = None,
    ply_thickness_mm: Optional[Union[float, Sequence[float]]] = None,
    panel: Optional[PanelGeometry] = None,
    intro: bool = True,
) -> DamageStackTimeline:
    """The damage to draw and its stack.

    ``source`` is an :class:`AnalysisConfig` (impact-driven: the impact
    mapping runs once, the tier does not matter; inspection-driven: its
    ``damage``) or a :class:`DamageState`, e.g. from
    :func:`~bvidfe.damage.io.load_cscan_json`. For a bare damage state,
    ``layup_deg`` and ``ply_thickness_mm`` give the plies (without them the
    plies are drawn without fibre directions) and ``panel`` the outline
    (default the 150 x 100 mm coupon, enlarged to hold the damage).
    ``intro=False`` leaves out the assembly, giving a loopable orbit.
    """
    schedule = INTRO_SCHEDULE if intro else ORBIT_SCHEDULE
    if isinstance(source, AnalysisConfig):
        if layup_deg is not None or ply_thickness_mm is not None or panel is not None:
            raise ValueError("layup_deg, ply_thickness_mm and panel come from the config")
        cfg = source
        lam = Laminate(_resolve_material(cfg.material), cfg.layup_deg, cfg.ply_thickness_mm)
        damage = BvidAnalysis(cfg)._resolve_damage(lam)
        onset = None
        point = None
        if cfg.impact is not None:
            onset = onset_energy_J(cfg)
            point = impact_point(cfg)
        return DamageStackTimeline(
            damage=damage,
            n_plies=len(cfg.layup_deg),
            panel=cfg.panel,
            layup_deg=tuple(float(a) for a in cfg.layup_deg),
            ply_thicknesses_mm=tuple(lam.ply_thicknesses_mm),
            impact=cfg.impact,
            impact_point_mm=point,
            onset_J=onset,
            config=cfg,
            source="impact model" if cfg.impact is not None else "supplied damage state",
            schedule=schedule,
        )
    if not isinstance(source, DamageState):
        raise ValueError(
            f"expected an AnalysisConfig or a DamageState (got {type(source).__name__})"
        )
    damage = source
    last = max((d.interface_index for d in damage.delaminations), default=-1)
    if layup_deg is not None:
        layup = tuple(float(a) for a in layup_deg)
        if not layup or not all(math.isfinite(a) for a in layup):
            raise ValueError("layup_deg must be a non-empty list of finite angles")
        if last > len(layup) - 2:
            raise ValueError(
                f"interface_index {last} is past the last interface of the "
                f"{len(layup)}-ply layup ({len(layup) - 2})"
            )
        n = len(layup)
    else:
        layup = None
        n = max(2, last + 2)
    thicknesses = None
    if ply_thickness_mm is not None:
        thicknesses = tuple(_normalise_ply_thicknesses(ply_thickness_mm, n))
    centroids = {d.centroid_mm for d in damage.delaminations}
    point = next(iter(centroids)) if len(centroids) == 1 else None
    return DamageStackTimeline(
        damage=damage,
        n_plies=n,
        panel=panel if panel is not None else _infer_panel(damage),
        layup_deg=layup,
        ply_thicknesses_mm=thicknesses,
        impact_point_mm=point,
        source="supplied damage state",
        schedule=schedule,
    )


# ---- Geometry -------------------------------------------------------------------


@dataclass(frozen=True)
class _Stack:
    """Drawn heights: ply ``k`` at ``z[k]`` (ply 0, the impact face, on top)."""

    z: np.ndarray
    gap: float
    height: float
    exaggeration: Optional[float]


def stack_region(tl: DamageStackTimeline) -> Tuple[Tuple[float, float, float, float], bool]:
    """``((x0, y0, x1, y1), zoomed)``: the part of the panel the plies are
    drawn over, the whole panel or, for damage small on it, a cut-out
    around the damage."""
    # Wide enough that the impactor does not hide the stack.
    min_half = 1.1 * tl.impact.impactor.diameter_mm if tl.impact is not None else 0.0
    return damage_region(
        tl.panel,
        tl.damage.delaminations,
        tl.impact_point_mm,
        tl.panel.Lx_mm / tl.panel.Ly_mm,
        min_half,
    )


def _stack(tl: DamageStackTimeline) -> _Stack:
    n = tl.n_plies
    (x0, y0, x1, y1), _ = stack_region(tl)
    height = STACK_HEIGHT_FRACTION * max(x1 - x0, y1 - y0)
    gap = height / max(1, n - 1)
    z = (n - 1 - np.arange(n)) * gap
    exag = None
    if tl.ply_thicknesses_mm:
        exag = gap / float(np.mean(tl.ply_thicknesses_mm))
    return _Stack(z=z, gap=gap, height=height, exaggeration=exag)


def _interface_z(st: _Stack, i: int) -> float:
    return float(0.5 * (st.z[i] + st.z[i + 1]))


def _ellipse_xy(d: DelaminationEllipse, n_pts: int = 72) -> np.ndarray:
    th = np.linspace(0.0, 2.0 * np.pi, n_pts, endpoint=False)
    c, s = math.cos(math.radians(d.orientation_deg)), math.sin(math.radians(d.orientation_deg))
    x, y = d.major_mm * np.cos(th), d.minor_mm * np.sin(th)
    return np.column_stack((c * x - s * y + d.centroid_mm[0], s * x + c * y + d.centroid_mm[1]))


def fibre_segments(
    angle_deg: float,
    box: Tuple[float, float, float, float],
    spacing: float,
) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Segments of parallel lines at ``angle_deg`` from x, ``spacing``
    apart, clipped to the rectangle ``box = (x0, y0, x1, y1)``."""
    x0, y0, x1, y1 = box
    a = math.radians(angle_deg)
    d = np.array([math.cos(a), math.sin(a)])
    nrm = np.array([-d[1], d[0]])
    corners = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], dtype=float)
    proj = corners @ nrm
    out = []
    c0 = math.ceil(proj.min() / spacing) * spacing
    for c in np.arange(c0, proj.max(), spacing):
        p0 = c * nrm
        lo, hi = -np.inf, np.inf
        for k, (vmin, vmax) in enumerate(((x0, x1), (y0, y1))):
            if abs(d[k]) < 1e-12:
                if not vmin <= p0[k] <= vmax:
                    lo, hi = 1.0, 0.0
                continue
            t1, t2 = (vmin - p0[k]) / d[k], (vmax - p0[k]) / d[k]
            lo, hi = max(lo, min(t1, t2)), min(hi, max(t1, t2))
        if hi - lo > 1e-9:
            a0, a1 = p0 + lo * d, p0 + hi * d
            out.append(((float(a0[0]), float(a0[1])), (float(a1[0]), float(a1[1]))))
    return out


SWATCH_FRACTION = 0.2  # fibre swatch side / shorter panel side


def corner_boxes(lx: float, ly: float) -> List[Tuple[float, float, float, float]]:
    """The four corner squares ``(x0, y0, x1, y1)`` of a ply that carry its
    fibre-direction swatch: corners (0, 0), (lx, 0), (lx, ly), (0, ly)."""
    w = SWATCH_FRACTION * min(lx, ly)
    return [(0.0, 0.0, w, w), (lx - w, 0.0, lx, w), (lx - w, ly - w, lx, ly), (0.0, ly - w, w, ly)]


def corner_weights(lx: float, ly: float, azim_deg: float) -> np.ndarray:
    """How strongly each corner's swatch is drawn at camera azimuth
    ``azim_deg``: 1 for the corner farthest from the camera, fading to 0
    for the others. A far corner projects above its ply's centre, behind
    the plies over it, so its swatch never covers a delamination below."""
    a = math.radians(azim_deg)
    eye = np.array([math.cos(a), math.sin(a)])
    corners = np.array([[-lx, -ly], [lx, -ly], [lx, ly], [-lx, ly]], dtype=float) / 2.0
    far = -(corners @ eye) / np.linalg.norm(corners, axis=1)
    return np.array([smoothstep(clamp01(1.0 - (far.max() - f) / 0.3)) for f in far])


def _lambert(X: np.ndarray, Y: np.ndarray, Z: np.ndarray, colour: str, alpha: float) -> np.ndarray:
    """Per-face colours of a surface lit from a fixed direction (two-sided:
    a face is lit whichever way its normal points)."""
    p = np.stack((X, Y, Z), axis=-1)
    n = np.cross(p[1:, :-1] - p[:-1, :-1], p[:-1, 1:] - p[:-1, :-1])
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    n = np.divide(n, norm, out=np.zeros_like(n), where=norm > 0)
    light = np.array([-0.45, -0.35, 0.82])
    light /= np.linalg.norm(light)
    shade = 0.55 + 0.45 * np.abs(n @ light)
    rgba = np.empty(shade.shape + (4,))
    rgba[..., :3] = np.clip(np.array(to_rgba(colour)[:3]) * shade[..., None], 0, 1)
    rgba[..., 3] = alpha
    return rgba


def _revolve(profile_r: np.ndarray, profile_z: np.ndarray, cx: float, cy: float, n: int = 40):
    th = np.linspace(0.0, 2.0 * np.pi, n)
    r, th2 = np.meshgrid(profile_r, th, indexing="ij")
    z = np.repeat(profile_z[:, None], n, axis=1)
    return cx + r * np.cos(th2), cy + r * np.sin(th2), z


def _impactor_surface(impact: ImpactEvent, cx: float, cy: float, z_tip: float):
    """``(X, Y, Z)`` of the impactor, its tip at height ``z_tip``."""
    imp = impact.impactor
    if imp.shape == "conical":
        h = 18.0
        r_top = h * math.tan(math.radians(30.0))
        r = np.linspace(0.0, r_top, 12)
        z = z_tip + r / math.tan(math.radians(30.0))
        r = np.append(r, [r_top * 0.999, 0.0])
        z = np.append(z, [z_tip + h, z_tip + h])
        return _revolve(r, z, cx, cy)
    radius = imp.diameter_mm / 2.0
    if imp.shape == "flat":
        r = np.array([0.0, radius, radius, 0.0])
        z = np.array([z_tip, z_tip, z_tip + 1.6 * radius, z_tip + 1.6 * radius])
        return _revolve(r, z, cx, cy)
    phi = np.linspace(0.0, np.pi, 24)
    return _revolve(radius * np.sin(phi), z_tip + radius * (1 - np.cos(phi)), cx, cy)


# ---- Drawing --------------------------------------------------------------------

_VIEW = (40, 160, 1380, 810)  # left, top, width, height of the 3D view
_KEY_X = 1440


def _ply_progress(tl: DamageStackTimeline, t: float, k: int) -> float:
    """0 (not yet placed) to 1 (in place) for ply ``k`` during the intro;
    the back face (last ply) is in place from the first frame and the impact
    face lands at about 85% of the intro, before the impactor comes down."""
    if not tl.intro:
        return 1.0
    u = tl.schedule.progress(t, "assemble")
    n = tl.n_plies
    m = 3.0  # each ply falls over 3 slots, one slot after the ply below
    rank = n - 1 - k  # 0 for the back face
    slot = u * max(1, n - 1) / 0.85
    return ease_out(clamp01((slot + m - rank) / m))


def _azimuth(tl: DamageStackTimeline, t: float) -> float:
    return AZIM0_DEG + 360.0 * tl.schedule.progress(t, "orbit")


def _draw_3d(fig: Figure, tl: DamageStackTimeline, t: float) -> None:
    left, top, width, height = _VIEW
    ax = cast(Axes3D, fig.add_axes(fr.px_rect(left, top, width, height), projection="3d"))
    ax.computed_zorder = False
    ax.set_proj_type("ortho")
    ax.set_axis_off()
    ax.set_facecolor(SURFACE)
    (ox, oy, x1, y1), _ = stack_region(tl)
    lx, ly = x1 - ox, y1 - oy
    st = _stack(tl)
    drop = 0.7 * st.height
    radius = tl.impact.impactor.diameter_mm / 2 if tl.impact is not None else 8.0
    top_space = 2.6 * radius + st.gap
    ax.set_xlim(ox, x1)
    ax.set_ylim(oy, y1)
    ax.set_zlim(-0.04 * st.height, st.height + top_space)
    ax.set_box_aspect((lx, ly, st.height * 1.04 + top_space), zoom=1.1)
    ax.view_init(elev=ELEV_DEG, azim=_azimuth(tl, t), roll=0)

    order = 0

    def z_of(k: int) -> float:
        return float(st.z[k] + drop * (1.0 - _ply_progress(tl, t, k)))

    by_iface: Dict[int, List[DelaminationEllipse]] = {}
    for d in tl.damage.delaminations:  # a C-scan may hold several per interface
        by_iface.setdefault(d.interface_index, []).append(d)
    rect = np.array([[ox, oy], [x1, oy], [x1, y1], [ox, y1]], dtype=float)
    azim = _azimuth(tl, t)
    weights = corner_weights(lx, ly, azim)
    for k in range(tl.n_plies - 1, -1, -1):  # bottom (back face) up
        p = _ply_progress(tl, t, k)
        if p <= 0:
            continue
        z = z_of(k)
        sheet = Poly3DCollection(
            [np.column_stack((rect, np.full(4, z)))],
            facecolors=[to_rgba(BAND, 0.10 * p)],
            edgecolors=[to_rgba(TEXT_2, 0.55 * p)],
            linewidths=fr.lw(1.2),
        )
        sheet.set_zorder(order)
        order += 1
        ax.add_collection3d(sheet)
        if tl.layup_deg is not None:
            for box, wgt in zip(corner_boxes(lx, ly), weights):
                if wgt <= 0.01:
                    continue
                bx0, by0, bx1, by1 = box[0] + ox, box[1] + oy, box[2] + ox, box[3] + oy
                sq = np.array([[bx0, by0, z], [bx1, by0, z], [bx1, by1, z], [bx0, by1, z]])
                square = Poly3DCollection(
                    [sq],
                    facecolors=[to_rgba(SURFACE, 0.9 * p * wgt)],
                    edgecolors=[to_rgba(TEXT_2, 0.7 * p * wgt)],
                    linewidths=fr.lw(1.0),
                )
                square.set_zorder(order)
                order += 1
                ax.add_collection3d(square)
                segs = fibre_segments(tl.layup_deg[k], (bx0, by0, bx1, by1), (bx1 - bx0) / 4.0)
                lines = Line3DCollection(
                    [[(a[0], a[1], z), (b[0], b[1], z)] for a, b in segs],
                    colors=[to_rgba(TEXT, 0.8 * p * wgt)],
                    linewidths=fr.lw(1.8),
                )
                lines.set_zorder(order)
                order += 1
                ax.add_collection3d(lines)
        if k == tl.n_plies - 1 and tl.damage.fiber_break_radius_mm > 0 and tl.impact_point_mm:
            cx, cy = tl.impact_point_mm
            th = np.linspace(0, 2 * np.pi, 48, endpoint=False)
            r_fb = tl.damage.fiber_break_radius_mm
            disc = Poly3DCollection(
                [np.column_stack((cx + r_fb * np.cos(th), cy + r_fb * np.sin(th), np.full(48, z)))],
                facecolors=[to_rgba(TEXT, 0.55 * p)],
                edgecolors=[to_rgba(TEXT, p)],
                linewidths=fr.lw(1.5),
            )
            disc.set_zorder(order)
            order += 1
            ax.add_collection3d(disc)
        # The interface above ply k (between ply k - 1 and k): shown once
        # the ply above it is in place.
        i = k - 1
        if i >= 0:
            pa = _ply_progress(tl, t, i)
            for d in by_iface.get(i, []):
                # At its final height (half a gap above ply k), fading in
                # as the ply above lands.
                zi = z + 0.5 * st.gap
                colour = series_color(i)
                poly = Poly3DCollection(
                    [np.column_stack((_ellipse_xy(d), np.full(72, zi)))],
                    facecolors=[to_rgba(colour, 0.72 * pa)],
                    edgecolors=[to_rgba(colour, pa)],
                    linewidths=fr.lw(2.2),
                    linestyles=[series_linestyle(i)],
                )
                poly.set_zorder(order)
                order += 1
                ax.add_collection3d(poly)

    # The dent on the impact face and the impactor above it.
    p_top = _ply_progress(tl, t, 0)
    point = tl.impact_point_mm
    light = LightSource(azdeg=315, altdeg=50)
    if point is not None and tl.damage.dent_depth_mm > 0 and p_top > 0:
        cx, cy = point
        r_d = max(radius, 4.0)
        exag = st.exaggeration if st.exaggeration is not None else 25.0
        # A shallow dimple: the real depth spans several plies.
        depth = min(tl.damage.dent_depth_mm * exag, 0.35 * st.gap, 0.25 * r_d)
        rr = np.linspace(0.0, r_d, 10)
        X, Y, Z = _revolve(rr, z_of(0) - depth * (1 - (rr / r_d) ** 2) ** 2, cx, cy, n=36)
        surf = ax.plot_surface(
            X,
            Y,
            Z,
            facecolors=_lambert(X, Y, Z, DENT_GREY, 0.9 * p_top),
            shade=False,
            linewidth=0,
            antialiased=False,
        )
        surf.set_zorder(order)
        order += 1
    if tl.impact is not None and point is not None:
        u_imp = 1.0
        if tl.intro:
            u_imp = ease_out(clamp01((tl.schedule.progress(t, "assemble") - 0.75) / 0.25))
        if u_imp > 0:
            z_tip = float(st.z[0]) + 0.35 * st.gap + drop * (1.0 - u_imp)
            X, Y, Z = _impactor_surface(tl.impact, point[0], point[1], z_tip)
            body = ax.plot_surface(
                X,
                Y,
                Z,
                color=IMPACTOR_GREY,
                alpha=u_imp,
                shade=True,
                lightsource=light,
                linewidth=0,
                antialiased=False,
            )
            body.set_zorder(order)
            order += 1
    # Face labels at the stack's left on screen: a point beyond the panel
    # opposite the camera's right vector projects to a fixed screen side.
    a = math.radians(azim)
    right = np.array([-math.sin(a), math.cos(a)])
    reach = 0.5 * math.hypot(lx, ly) * 1.02
    lab_x, lab_y = np.array([ox + lx / 2, oy + ly / 2]) - reach * right
    labels = ((0, "impact face, ply 0"), (tl.n_plies - 1, f"back face, ply {tl.n_plies - 1}"))
    for k, label in labels:
        if _ply_progress(tl, t, k) >= 1.0:
            ax.text(
                float(lab_x),
                float(lab_y),
                z_of(k),
                label + "  ",
                fontsize=fr.pt(fr.CAPTION_PX),
                color=TEXT_2,
                ha="right",
                va="center",
                zorder=order + 1,
            )
    if not tl.damage.delaminations:
        ax.text2D(
            0.5,
            0.12,
            _no_damage_text(tl),
            transform=ax.transAxes,
            fontsize=fr.pt(28),
            color=TEXT_2,
            ha="center",
            va="center",
            bbox={"facecolor": BAND, "edgecolor": "none", "pad": 6},
        )


def _no_damage_text(tl: DamageStackTimeline) -> str:
    if tl.onset_J is not None:
        return f"below damage onset ({tl.onset_J:.1f} J): no delamination"
    return "no delamination in this damage state"


def _key(fig: Figure, tl: DamageStackTimeline, t: float) -> None:
    x = _KEY_X
    d = tl.damage
    fr.text_px(fig, x, 186, "Damage", fr.PANEL_PX, TEXT, weight="bold", ha="left")
    rows = [
        ("damage area", f"{d.projected_damage_area_mm2:,.0f} mm²"),
        ("dent depth", f"{d.dent_depth_mm:.2f} mm"),
        ("delaminated interfaces", f"{len({e.interface_index for e in d.delaminations})}"),
    ]
    if tl.onset_J is not None:
        rows.append(("damage onset", f"{tl.onset_J:.1f} J"))
    if d.fiber_break_radius_mm > 0:
        rows.append(("fibre-break core radius", f"{d.fiber_break_radius_mm:.1f} mm"))
    y = 222.0
    for label, value in rows:
        fr.text_px(fig, x, y, label, fr.CAPTION_PX + 2, TEXT_2, ha="left")
        fr.text_px(fig, fr.W - fr.MARGIN, y, value, fr.CAPTION_PX + 2, TEXT, ha="right")
        y += 28
    y += 22
    fr.text_px(fig, x, y, "Interface", fr.AXIS_PX, TEXT, weight="bold", ha="left")
    fr.text_px(fig, fr.W - fr.MARGIN, y, "2a × 2b", fr.CAPTION_PX, TEXT_2, ha="right")
    y += 24
    fr.text_px(
        fig,
        x,
        y,
        "plies above / below" if tl.layup_deg is not None else "ply angles not given",
        fr.CAPTION_PX,
        TEXT_2,
        ha="left",
    )
    y += 18
    dels = sorted(d.delaminations, key=lambda e: e.interface_index)
    bottom = 950.0
    step = min(30.0, (bottom - y) / max(1, len(dels)))
    ax = fr.px_axes(fig, x, y, fr.W - fr.MARGIN - x, bottom - y)
    ax.set_xlim(0, fr.W - fr.MARGIN - x)
    ax.set_ylim(bottom - y, 0)
    ax.set_axis_off()
    for k, e in enumerate(dels):
        yy = 14 + k * step
        idx = e.interface_index
        p = _ply_progress(tl, t, idx)
        ax.plot(
            [0, 44],
            [yy, yy],
            color=to_rgba(series_color(idx), max(0.25, p)),
            linewidth=fr.lw(3.5),
            linestyle=series_linestyle(idx),
            solid_capstyle="butt",
            dash_capstyle="butt",
        )
        if tl.layup_deg is not None:
            label = f"{idx:>2}   {tl.layup_deg[idx]:g}° / {tl.layup_deg[idx + 1]:g}°"
        else:
            label = f"{idx:>2}"
        ax.text(58, yy, label, fontsize=fr.pt(fr.CAPTION_PX), color=TEXT_2, va="center")
        ax.text(
            fr.W - fr.MARGIN - x,
            yy,
            f"{2 * e.major_mm:.0f} × {2 * e.minor_mm:.0f} mm",
            fontsize=fr.pt(fr.CAPTION_PX),
            color=TEXT_2,
            va="center",
            ha="right",
        )


def _scale_note(tl: DamageStackTimeline) -> str:
    st = _stack(tl)
    if st.exaggeration is not None:
        note = (
            f"Not to scale through the thickness: plies spread ×{st.exaggeration:.0f}; "
            "the dent is drawn on the impact face, its depth shortened"
        )
    else:
        note = "Not to scale through the thickness: plies spread apart"
    (x0, y0, x1, y1), zoomed = stack_region(tl)
    if zoomed:
        note += (
            f".  Plies cut to {x1 - x0:.0f} × {y1 - y0:.0f} mm around the damage "
            f"(panel {tl.panel.Lx_mm:g} × {tl.panel.Ly_mm:g} mm)"
        )
    return note


def draw_frame(fig: Figure, timeline: DamageStackTimeline, t: float) -> None:
    """Draw the frame at clip time ``t`` (0 to 1) into ``fig`` (cleared
    first). Pure: no model call."""
    tl = timeline
    cfg = tl.config
    with fr.style_context():
        fig.clear()
        fig.set_facecolor(SURFACE)
        if cfg is not None:
            parts = [f"{fr.material_name(cfg)}  {fr.layup_label(cfg.layup_deg)}"]
            if cfg.impact is not None:
                parts.append(fr.impactor_label(cfg.impact))
            parts.append(f"{cfg.panel.Lx_mm:g} × {cfg.panel.Ly_mm:g} mm panel")
            subtitle = "  ·  ".join(parts)
            caveats = fr.honesty_notes(cfg)
            summary = fr.config_summary(cfg)
        else:
            lay = fr.layup_label(tl.layup_deg) if tl.layup_deg else f"{tl.n_plies} plies"
            subtitle = (
                f"Supplied damage state  ·  {lay}  ·  "
                f"{tl.panel.Lx_mm:g} × {tl.panel.Ly_mm:g} mm panel"
            )
            caveats = []
            summary = subtitle
        if tl.impact is not None:
            readout = f"{tl.impact.energy_J:.1f} J"
            readout_sub = "impact energy  ·  damage from the impact model (every tier)"
        else:
            readout = "C-scan" if cfg is None else "inspection"
            readout_sub = "damage as supplied (no impact model)"
        fr.draw_title_band(
            fig,
            "Delaminations through the thickness, exploded view",
            subtitle,
            readout=readout,
            readout_sub=readout_sub,
        )
        fr.draw_caveats(fig, caveats)
        _draw_3d(fig, tl, t)
        fr.text_px(fig, fr.MARGIN, 986, _scale_note(tl), fr.CAPTION_PX, TEXT_2, ha="left")
        _key(fig, tl, t)
        fr.draw_footer(fig, summary)
