"""A3: the fe3d load ramp of a damaged panel, to its residual strength.

fe3d is linear elastic up to first-ply failure, so one solve per mesh (the
damaged panel and the undamaged one it is normalised by) gives every load
step: at an applied stress ``sigma`` the displacement is ``sigma /
sigma_ref`` times the reference one, and each element's load ratio to its
own failure is ``(sigma / sigma_ref) / c_e`` (:mod:`bvidfe.analysis.fe_fields`).
The ratio uses the disk-averaged multiplier by default, so the first
element to reach 1 is the one at which fe3d sets first-ply failure.

Left, the damaged panel in 3D, impact face up, cut open along the load
through the first element to fail: only the faces the camera sees are
drawn (the impact face, the cut face and one end face of the remaining
block), so mplot3d's depth sorting cannot misorder them. Each element is
coloured by its load ratio, the deformation and the thickness are
exaggerated (both factors printed), and the delamination footprint (the
fe3d damage zone) is outlined on the impact face. Below it, top views of
the highest ratio through the thickness, damaged and undamaged at the same
applied stress. Right, stress against end shortening, the highest ratio of
each panel against 1, and the results ``BvidAnalysis(cfg).run()`` reports
for the config with ``tier="fe3d"``.

When delamination growth governs (its closed-form stress below first-ply
failure, as in ``run()``), the ramp ends at the growth stress, with no FE
field of the growth itself. Stresses are printed on the scale ``run()``
reports (fe3d's applied stress times the pristine strength over the
undamaged panel's first-ply failure).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar, List, Literal, Optional, Tuple, cast

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import Normalize, to_rgba
from matplotlib.figure import Figure
from matplotlib.patheffects import AbstractPathEffect, withStroke
from mpl_toolkits.mplot3d import Axes3D, proj3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from bvidfe.analysis.config import AnalysisConfig, MeshParams
from bvidfe.analysis.fe_fields import Fe3dFields, fe3d_field_history
from bvidfe.analysis.fe_mesh import estimate_fe_mesh_size
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.viz.anim import frame as fr
from bvidfe.viz.anim.timeline import Phase, Schedule, Timeline
from bvidfe.viz.style import (
    BASELINE,
    CATEGORICAL,
    GRID,
    STATUS,
    SURFACE,
    TEXT,
    TEXT_2,
    sequential_cmap,
)

SCENE = "cai-ramp"
RatioKind = Literal["averaged", "pointwise"]

DEFAULT_SCHEDULE = Schedule(
    (
        Phase("hold_start", 0.05),
        Phase("ramp", 0.80, "ease_in_out"),
        Phase("hold_end", 0.15),  # 1.5 s of a 10 s clip at the residual strength
    )
)

#: Most faces the 3D view draws (impact face, cut face and end face of the
#: block it shows), about 0.4 s a frame. A 150 x 100 mm, 16-ply panel draws
#: about 1,100 faces at 5 mm and 2,800 at 2.5 mm; the fe3d size cap
#: (``FE3D_MAX_DOF``) keeps most meshes well below this.
MAX_FACES = 12_000

ELEV_DEG = 24.0
AZIM_DEG = -58.0  # camera on the -y side; +180 deg when the block kept is y below the cut
DEFORMATION_TARGET = 0.08  # largest drawn displacement / longest panel side, at the ramp's end
THICKNESS_TARGET = 0.12  # drawn thickness / longest panel side
FAILED = 1.0 - 1e-9  # a load ratio at or above this is drawn as failed

# Drawing time of one 1080p frame on one core (s), measured on a Linux
# container: 0.26 s at 1,100 and 2,800 faces, 0.30 s at 5,200.
_FRAME_S_BASE = 0.24
_FRAME_S_PER_FACE = 1.2e-5


# ---- Timeline -------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class CaiRampView:
    """What the 3D view draws: the visible faces of the block of the
    damaged mesh on one side of the cut, and the damage-zone outline.

    ``faces`` are ``(n, 4)`` node indices and ``parent`` their elements:
    the impact-face quads of the kept rows, the cut face (row ``cut_row``,
    every layer) and the end face the camera sees. ``zone_xy`` holds the
    footprint outline's polylines on the kept side (mm, undeformed) and
    ``zone_u_ref`` the reference displacement at their points, interpolated
    from the impact-face nodes.
    """

    faces: np.ndarray
    parent: np.ndarray
    n_top: int
    n_cut: int
    cut_row: int
    cut_y_mm: float
    keep_high: bool
    azim_deg: float
    deformation_scale: float
    thickness_scale: float
    zone_xy: Tuple[np.ndarray, ...]
    zone_u_ref: Tuple[np.ndarray, ...]
    zone_full_xy: Tuple[np.ndarray, ...]
    x_edges: np.ndarray
    y_edges: np.ndarray
    c_column_damaged: np.ndarray
    c_column_undamaged: np.ndarray

    @property
    def n_faces(self) -> int:
        return int(self.faces.shape[0])


@dataclass(frozen=True, eq=False)
class CaiRampTimeline(Timeline):
    """The fe3d fields of a config (two solves, four on a biaxial config)
    and the clip schedule: hold unloaded, ramp the applied stress to where
    the damaged panel reaches its residual strength, hold there.

    ``ratio`` is the multiplier the colours use: ``"averaged"`` (the
    disk-averaged check that sets first-ply failure) or ``"pointwise"``
    (the lowest Gauss point of each element, which can exceed 1 before
    first-ply failure)."""

    scene: ClassVar[str] = SCENE

    config: AnalysisConfig
    fields: Fe3dFields
    view: CaiRampView
    ratio: RatioKind = "averaged"
    schedule: Schedule = DEFAULT_SCHEDULE

    @property
    def end_stress_MPa(self) -> float:
        """Applied (fe3d) stress at the end of the ramp."""
        return float(self.fields.residual_applied_stress_MPa)

    def stress_at(self, t: float) -> float:
        """The applied (fe3d) stress at clip time ``t``: 0 to
        ``end_stress_MPa``."""
        return self.end_stress_MPa * self.schedule.progress(t, "ramp")

    def reported_stress_at(self, t: float) -> float:
        """The stress at clip time ``t`` on the scale ``run()`` reports; the
        residual strength itself at the end of the ramp."""
        if self.schedule.progress(t, "ramp") >= 1.0:
            return float(self.fields.residual_strength_MPa)
        return float(self.fields.reported_stress(self.stress_at(t)))

    def load_ratio(self, stress_MPa: float, damaged: bool = True) -> np.ndarray:
        """Per element, the load ratio to local failure at ``stress_MPa``
        (applied scale), with this timeline's multiplier."""
        panel = self.fields.damaged if damaged else self.fields.undamaged
        return panel.load_ratio(stress_MPa, pointwise=self.ratio == "pointwise")

    def max_ratio(self, stress_MPa: float, damaged: bool = True) -> float:
        return float(self.load_ratio(stress_MPa, damaged).max())

    @property
    def tier(self) -> str:
        return "fe3d"

    def frame_seconds_estimate(self) -> float:
        """Drawing time of one frame on one core (s), from the face count."""
        return _FRAME_S_BASE + _FRAME_S_PER_FACE * self.view.n_faces

    def draw(self, fig: Figure, t: float) -> None:
        draw_frame(fig, self, t)


# ---- Model -----------------------------------------------------------------


def face_count_bound(cfg: AnalysisConfig) -> int:
    """An upper bound on the faces the 3D view of ``cfg``'s fe3d mesh draws
    (without solving): the whole impact face, one long side and one end,
    each element face through the thickness."""
    stats = estimate_fe_mesh_size(cfg)
    nx, ny, nz = stats["nx"], stats["ny"], stats["nz"]
    return int(nx * ny + nx * nz + ny * nz)


def check_config(cfg: AnalysisConfig, max_faces: int = MAX_FACES) -> None:
    """``ValueError`` unless ``cfg`` is an :class:`AnalysisConfig` whose
    3D view stays within ``max_faces``."""
    if not isinstance(cfg, AnalysisConfig):
        raise ValueError(f"expected an AnalysisConfig (got {type(cfg).__name__})")
    n = face_count_bound(cfg)
    if n > max_faces:
        size = (cfg.mesh or MeshParams()).in_plane_size_mm
        raise ValueError(
            f"the fe3d load ramp would draw up to {n:,} faces a frame at the "
            f"{size:g} mm mesh, more than the {max_faces:,} the animation allows "
            f"(about {_FRAME_S_BASE + _FRAME_S_PER_FACE * n:.1f} s a frame); use a "
            "coarser in-plane mesh size (MeshParams.in_plane_size_mm) or fewer "
            "elements per ply"
        )


def _nice_floor(v: float) -> float:
    """The largest of 1, 2, 5 x 10^n not above ``v`` (at least 1)."""
    if not math.isfinite(v) or v <= 1.0:
        return 1.0
    p = 10.0 ** math.floor(math.log10(v))
    for m in (5.0, 2.0, 1.0):
        if m * p <= v:
            return m * p
    return p  # pragma: no cover


def _ellipse_polygon(d: DelaminationEllipse):
    from shapely import affinity
    from shapely.geometry import Point

    unit = Point(d.centroid_mm).buffer(1.0, 32)
    return affinity.rotate(
        affinity.scale(unit, d.major_mm, d.minor_mm, origin=d.centroid_mm),
        d.orientation_deg,
        origin=d.centroid_mm,
    )


def _lines(geom) -> List[np.ndarray]:
    """The polylines of a (multi)line geometry, each ``(n, 2)``."""
    if geom.is_empty:
        return []
    parts = getattr(geom, "geoms", [geom])
    out = []
    for g in parts:
        if g.geom_type in ("LineString", "LinearRing"):
            xy = np.asarray(g.coords, dtype=float)
            if len(xy) >= 2:
                out.append(xy[:, :2])
        else:
            out += _lines(g)
    return out


def damage_zone_outline(
    damage: DamageState, box: Tuple[float, float, float, float]
) -> List[np.ndarray]:
    """The outline of the delamination footprint (the union of every
    interface's ellipse, fe3d's damage zone) inside ``box = (x0, y0, x1,
    y1)``, as polylines in mm."""
    if not damage.delaminations:
        return []
    from shapely.geometry import box as shp_box
    from shapely.ops import unary_union

    zone = unary_union([_ellipse_polygon(d) for d in damage.delaminations])
    return _lines(zone.boundary.intersection(shp_box(*box)))


def _impact_face_interpolator(fields: Fe3dFields, nx: int, ny: int):
    """``f(xy) -> (n, 3)``: the damaged panel's reference displacement on
    the impact face (z = 0), bilinear in the structured node grid."""
    from scipy.interpolate import RegularGridInterpolator

    mesh = fields.damaged.mesh
    n_layer = (nx + 1) * (ny + 1)
    xyz = mesh.node_coords[:n_layer]
    xs = xyz[: nx + 1, 0]
    ys = xyz[:: nx + 1, 1][: ny + 1]
    u = fields.damaged.u_ref.reshape(-1, 3)[:n_layer].reshape(ny + 1, nx + 1, 3)
    interp = RegularGridInterpolator((ys, xs), u, bounds_error=False, fill_value=None)

    def at(xy: np.ndarray) -> np.ndarray:
        return np.asarray(interp(np.column_stack((xy[:, 1], xy[:, 0]))), dtype=float)

    return at


def _view(fields: Fe3dFields, ratio: RatioKind) -> CaiRampView:
    d, u = fields.damaged, fields.undamaged
    mesh = d.mesh
    nz, ny, nx = d.grid_shape
    conn = mesh.element_connectivity
    e = np.arange(mesh.n_elements)
    k_idx, rem = np.divmod(e, ny * nx)
    j_idx, i_idx = np.divmod(rem, nx)
    _, j_f, _ = _unravel(d.fpf_element, d.grid_shape)
    cfg = fields.config
    lx, ly = cfg.panel.Lx_mm, cfg.panel.Ly_mm
    dels = fields.damage.delaminations
    cy = float(np.mean([el.centroid_mm[1] for el in dels])) if dels else ly / 2.0
    row_lo = float(mesh.node_coords[conn[d.fpf_element], 1].min())
    row_hi = float(mesh.node_coords[conn[d.fpf_element], 1].max())
    keep_high = cy >= 0.5 * (row_lo + row_hi)
    azim = AZIM_DEG if keep_high else AZIM_DEG + 180.0
    kept = (j_idx >= j_f) if keep_high else (j_idx <= j_f)
    top = kept & (k_idx == 0)
    cut = j_idx == j_f
    # The end the camera sees: x high from the -y side, x low from +y.
    end = kept & (i_idx == (nx - 1 if keep_high else 0))
    cut_face = [0, 1, 5, 4] if keep_high else [3, 2, 6, 7]
    end_face = [1, 2, 6, 5] if keep_high else [0, 3, 7, 4]
    faces = np.concatenate(
        (conn[top][:, [0, 1, 2, 3]], conn[cut][:, cut_face], conn[end][:, end_face])
    ).astype(np.int64)
    parent = np.concatenate((e[top], e[cut], e[end])).astype(np.int64)
    cut_y = row_lo if keep_high else row_hi

    # Deformation and thickness exaggeration, from the ramp's end state.
    u_end = d.displacement(fields.residual_applied_stress_MPa)
    u_max = float(np.abs(u_end[:, :2]).max()) if u_end.size else 0.0
    side = max(lx, ly)
    deform = _nice_floor(DEFORMATION_TARGET * side / u_max) if u_max > 0 else 1.0
    h = float(mesh.node_coords[:, 2].max())
    thick = max(1.0, float(round(THICKNESS_TARGET * side / h))) if h > 0 else 1.0

    box = (0.0, cut_y, lx, ly) if keep_high else (0.0, 0.0, lx, cut_y)
    zone = damage_zone_outline(fields.damage, box)
    full = damage_zone_outline(fields.damage, (0.0, 0.0, lx, ly))
    at = _impact_face_interpolator(fields, nx, ny)
    zone_u = tuple(at(z) for z in zone)

    n_layer = (nx + 1) * (ny + 1)
    xyz = mesh.node_coords[:n_layer]
    pointwise = ratio == "pointwise"
    c_d = d.c_pointwise if pointwise else d.c_averaged
    c_u = u.c_pointwise if pointwise else u.c_averaged
    return CaiRampView(
        faces=_frozen(faces),
        parent=_frozen(parent),
        n_top=int(top.sum()),
        n_cut=int(cut.sum()),
        cut_row=int(j_f),
        cut_y_mm=float(cut_y),
        keep_high=bool(keep_high),
        azim_deg=float(azim),
        deformation_scale=float(deform),
        thickness_scale=float(thick),
        zone_xy=tuple(_frozen(z) for z in zone),
        zone_u_ref=tuple(_frozen(z) for z in zone_u),
        zone_full_xy=tuple(_frozen(z) for z in full),
        x_edges=_frozen(np.array(xyz[: nx + 1, 0], dtype=float)),
        y_edges=_frozen(np.array(xyz[:: nx + 1, 1][: ny + 1], dtype=float)),
        c_column_damaged=_frozen(np.asarray(c_d).reshape(nz, ny, nx).min(axis=0)),
        c_column_undamaged=_frozen(np.asarray(c_u).reshape(nz, ny, nx).min(axis=0)),
    )


def _frozen(a: np.ndarray) -> np.ndarray:
    a = np.array(a, copy=True)
    a.setflags(write=False)
    return a


def _unravel(e: int, grid: Tuple[int, int, int]) -> Tuple[int, int, int]:
    """``(k, j, i)`` of element ``e`` in the ``(nz, ny, nx)`` grid."""
    _, ny, nx = grid
    k, rem = divmod(int(e), ny * nx)
    j, i = divmod(rem, nx)
    return k, j, i


def build_timeline(
    cfg: AnalysisConfig,
    damage: Optional[DamageState] = None,
    *,
    ratio: RatioKind = "averaged",
    schedule: Schedule = DEFAULT_SCHEDULE,
    max_faces: int = MAX_FACES,
) -> CaiRampTimeline:
    """Solve ``cfg``'s fe3d meshes once (:func:`fe3d_field_history`: two
    solves, four on a biaxial config) and return the timeline.

    ``cfg.tier`` is not read: the ramp is fe3d's whatever the tier, on
    ``cfg.mesh`` (default :class:`MeshParams`). ``damage`` defaults to the
    damage ``run()`` resolves (``cfg.damage`` or the impact mapping).
    ``ratio`` picks the multiplier the colours use (see
    :class:`CaiRampTimeline`). Raises ``ValueError`` when the 3D view would
    draw more than ``max_faces`` faces (checked before solving), and
    ``FE3DSizeError`` above the fe3d size cap, as ``run()`` does.
    """
    if ratio not in ("averaged", "pointwise"):
        raise ValueError(f"ratio must be 'averaged' or 'pointwise' (got {ratio!r})")
    check_config(cfg, max_faces)
    fields = fe3d_field_history(cfg, damage)
    view = _view(fields, ratio)
    if view.n_faces > max_faces:  # pragma: no cover - the bound above holds
        raise ValueError(f"{view.n_faces:,} faces is more than {max_faces:,}")
    return CaiRampTimeline(config=cfg, fields=fields, view=view, ratio=ratio, schedule=schedule)


# ---- Drawing -----------------------------------------------------------------

# Layout, px on the 1920 x 1080 canvas.
_VIEW_BOX = (40, 258, 1150, 392)  # left, top, width, height of the 3D view
_MAPS_TOP = 784
_MAP_H = 196
_MAP_W_MAX = 330
_MAP_LEFTS = (96, 486)
_CBAR = (872, 784, 22, 196)
_KEY_X = 994
_RIGHT_X = 1268
_CURVE = (1340, 252, 520, 268)
_GAUGE_TOP = 600
_RESULT_TOP = 776

_CMAP = sequential_cmap()
_NORM = Normalize(0.0, 1.0)
_CRITICAL = STATUS["critical"]
_HALO: List[AbstractPathEffect] = [withStroke(linewidth=fr.lw(5.0), foreground=SURFACE)]


def _load_word(cfg: AnalysisConfig) -> str:
    return "CAI" if cfg.loading == "compression" else "TAI"


def _ply_label(tl: CaiRampTimeline, e: int) -> str:
    mesh = tl.fields.damaged.mesh
    return f"ply {int(mesh.ply_indices[e])} ({float(mesh.ply_angles_deg[e]):g}°)"


def _element_centre(tl: CaiRampTimeline, e: int) -> Tuple[float, float]:
    xyz = tl.fields.damaged.mesh.node_coords[tl.fields.damaged.mesh.element_connectivity[e]]
    return float(xyz[:, 0].mean()), float(xyz[:, 1].mean())


def _drawn_points(tl: CaiRampTimeline, xyz: np.ndarray, u: np.ndarray, m: float) -> np.ndarray:
    """Drawn coordinates of points at ``xyz`` with reference displacement
    ``u`` at multiplier ``m``: displaced by the exaggerated displacement,
    the impact face (z = 0) on top and the thickness stretched."""
    v = tl.view
    h = float(tl.fields.damaged.mesh.node_coords[:, 2].max())
    p = xyz + v.deformation_scale * m * u
    p[:, 2] = (h - p[:, 2]) * v.thickness_scale
    return p


def _project(ax: Axes3D, pts: np.ndarray) -> np.ndarray:
    """Projected 2D data coordinates of 3D points in ``ax``."""
    x, y, _ = proj3d.proj_transform(pts[:, 0], pts[:, 1], pts[:, 2], ax.get_proj())
    return np.column_stack((x, y))


def _fit_zoom(
    ax: Axes3D, corners: np.ndarray, aspect: Tuple[float, ...], box_px: Tuple[float, float]
) -> None:
    """Zoom the 3D axes so the projected ``corners`` fill 94% of a box of
    ``box_px`` (width, height in design pixels) at the axes' centre."""
    ax.set_box_aspect(aspect, zoom=1.0)
    xy = _project(ax, corners)
    # The square axes spans its 2D view limits; projected units -> px.
    lim = ax.viewLim
    side = float(_VIEW_BOX[2])
    w = float(np.ptp(xy[:, 0])) / lim.width * side or 1.0
    h = float(np.ptp(xy[:, 1])) / lim.height * side or 1.0
    ax.set_box_aspect(aspect, zoom=0.94 * min(box_px[0] / w, box_px[1] / h))


def _draw_3d(fig: Figure, tl: CaiRampTimeline, stress: float, failed: bool) -> None:
    f, v = tl.fields, tl.view
    d = f.damaged
    mesh = d.mesh
    cfg = tl.config
    m = d.multiplier(stress)
    left, top, width, height = _VIEW_BOX
    # mplot3d draws in a square: a square of the view box's width, centred
    # on it, transparent, the drawing zoomed to fill the box.
    cy = top + height / 2.0
    ax = cast(
        Axes3D, fig.add_axes(fr.px_rect(left, cy - width / 2.0, width, width), projection="3d")
    )
    ax.patch.set_visible(False)
    ax.computed_zorder = False
    ax.set_proj_type("ortho")
    ax.set_axis_off()
    ax.set_facecolor(SURFACE)
    lx, ly = cfg.panel.Lx_mm, cfg.panel.Ly_mm
    h = float(mesh.node_coords[:, 2].max())
    y0, y1 = (v.cut_y_mm, ly) if v.keep_high else (0.0, v.cut_y_mm)
    zt = h * v.thickness_scale
    pad = 0.04 * max(lx, ly)
    ax.set_xlim(-pad, lx + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_zlim(0.0, zt)
    ax.view_init(elev=ELEV_DEG, azim=v.azim_deg, roll=0)
    aspect = (lx + 2 * pad, (y1 - y0) + 2 * pad, zt)
    corners = np.array(
        [[x, y, z] for x in (-pad, lx + pad) for y in (y0, y1) for z in (0.0, zt)], dtype=float
    )
    _fit_zoom(ax, corners, aspect, (width, height))

    nodes = _drawn_points(tl, mesh.node_coords, d.u_ref.reshape(-1, 3), m)
    ratio = tl.load_ratio(stress)[v.parent]
    colours = _CMAP(_NORM(np.clip(ratio, 0.0, 1.0)))
    edge_alpha = 0.45 if v.n_faces < 2500 else 0.25
    polys = Poly3DCollection(
        nodes[v.faces],
        facecolors=colours,
        edgecolors=[to_rgba(SURFACE, edge_alpha)],
        linewidths=fr.lw(0.8),
    )
    polys.set_zorder(1)
    ax.add_collection3d(polys)
    # The block's silhouette edges, in ink.
    outline = _block_outline(tl, nodes)
    rim = Line3DCollection(outline, colors=[to_rgba(TEXT_2, 0.8)], linewidths=fr.lw(1.4))
    rim.set_zorder(2)
    ax.add_collection3d(rim)

    # The damage zone on the impact face (z = 0), following the deformation.
    for xy, uz in zip(v.zone_xy, v.zone_u_ref):
        pts = _drawn_points(tl, np.column_stack((xy, np.zeros(len(xy)))), uz, m)
        (ln,) = ax.plot(
            pts[:, 0],
            pts[:, 1],
            pts[:, 2],
            color=TEXT,
            linewidth=fr.lw(2.2),
            linestyle=(0, (5, 3)),
            zorder=3,
        )
        ln.set_path_effects(_HALO)

    # Failed elements (ratio >= 1): their visible faces outlined in the
    # critical status colour; before that, the first element to fail in ink.
    hot = ratio >= FAILED
    if hot.any():
        ring = Poly3DCollection(
            nodes[v.faces[hot]],
            facecolors="none",
            edgecolors=[_CRITICAL],
            linewidths=fr.lw(3.5),
        )
        ring.set_zorder(4)
        ax.add_collection3d(ring)
    first = (v.parent == d.fpf_element) & ~hot
    if first.any():
        ring = Poly3DCollection(
            nodes[v.faces[first]],
            facecolors="none",
            edgecolors=[TEXT],
            linewidths=fr.lw(2.0),
        )
        ring.set_zorder(4)
        ax.add_collection3d(ring)

    _annotate_3d(ax, tl, nodes, failed)


def _block_outline(tl: CaiRampTimeline, nodes: np.ndarray) -> List[np.ndarray]:
    """The edges of the drawn block: the impact face's border and the cut
    face's and end face's outer edges, as polylines of drawn points."""
    v = tl.view
    nz, ny, nx = tl.fields.damaged.grid_shape
    n_layer = (nx + 1) * (ny + 1)

    def node(kz: int, jy: int, ix: int) -> int:
        return kz * n_layer + jy * (nx + 1) + ix

    j0 = v.cut_row if v.keep_high else 0
    j1 = ny if v.keep_high else v.cut_row + 1
    jc = v.cut_row if v.keep_high else v.cut_row + 1  # the node row of the cut
    ie = nx if v.keep_high else 0  # the node column of the end face
    lines = []
    # Impact face border (z = 0 nodes, layer 0).
    ring = (
        [node(0, j0, i) for i in range(nx + 1)]
        + [node(0, j, nx) for j in range(j0, j1 + 1)]
        + [node(0, j1, i) for i in range(nx, -1, -1)]
        + [node(0, j, 0) for j in range(j1, j0 - 1, -1)]
    )
    lines.append(nodes[ring])
    # The cut face's bottom edge (back face) and the vertical edges.
    lines.append(nodes[[node(nz, jc, i) for i in range(nx + 1)]])
    for ix in (0, nx):
        lines.append(nodes[[node(kz, jc, ix) for kz in range(nz + 1)]])
    far = j1 if v.keep_high else j0
    lines.append(nodes[[node(kz, far, ie) for kz in range(nz + 1)]])
    lines.append(
        nodes[
            [node(nz, j, ie) for j in (range(jc, far + 1) if v.keep_high else range(far, jc + 1))]
        ]
    )
    return lines


def _annotate_3d(ax: Axes3D, tl: CaiRampTimeline, nodes: np.ndarray, failed: bool) -> None:
    """Labels on the 3D view: the impact face, the cut, the first element
    to fail and the load arrows."""
    f, v = tl.fields, tl.view
    d = f.damaged
    conn = d.mesh.element_connectivity
    e_f = d.fpf_element
    # The first element to fail: the centre of its face on the cut.
    on_cut = np.flatnonzero(v.parent[v.n_top : v.n_top + v.n_cut] == e_f)
    face = v.faces[v.n_top + int(on_cut[0])] if on_cut.size else conn[e_f][:4]
    centre = _project(ax, nodes[face].mean(axis=0, keepdims=True))[0]
    pts2 = _project(ax, nodes[v.faces.ravel()])
    # 30 px below the lowest point of the block (projected units per px:
    # the square axes spans its view limits).
    below = float(pts2[:, 1].min()) - 30.0 * ax.viewLim.height / float(_VIEW_BOX[2])
    colour = _CRITICAL if failed else TEXT
    label = (
        f"first-ply failure: {_ply_label(tl, e_f)}"
        if failed
        else f"first element to fail: {_ply_label(tl, e_f)}"
    )
    if tl.fields.governing == "delamination_growth" and not failed:
        label = f"highest load ratio: {_ply_label(tl, e_f)}"
    ax.annotate(
        label,
        xy=(float(centre[0]), float(centre[1])),
        xytext=(float(centre[0]), below),
        textcoords="data",
        fontsize=fr.pt(fr.CAPTION_PX + 2),
        color=colour,
        weight="bold" if failed else "normal",
        ha="center",
        va="top",
        arrowprops={
            "arrowstyle": "-",
            "color": colour,
            "linewidth": fr.lw(1.6),
            "shrinkA": 2,
            "shrinkB": 3,
        },
        zorder=6,
        annotation_clip=False,
    )
    # The impact face, labelled at its far left corner on screen.
    cfg = tl.config
    lx, ly = cfg.panel.Lx_mm, cfg.panel.Ly_mm
    far_y = ly if v.keep_high else 0.0
    _, ny, nx = d.grid_shape
    corner = (ny * (nx + 1)) if v.keep_high else nx  # node (x 0, y max) or (x max, y 0)
    cxy = _project(ax, nodes[[corner]])[0]
    ax.annotate(
        "impact face, ply 0",
        xy=(float(cxy[0]), float(cxy[1])),
        xytext=(-10, 0),
        textcoords="offset points",
        fontsize=fr.pt(fr.CAPTION_PX),
        color=TEXT_2,
        ha="right",
        va="center",
        annotation_clip=False,
    )
    # Load arrows at mid-thickness on both ends: inward in compression,
    # outward in tension.
    y_mid = 0.5 * (v.cut_y_mm + far_y)
    zt = float(d.mesh.node_coords[:, 2].max()) * v.thickness_scale
    tension = cfg.loading != "compression"
    reach = 0.10 * lx
    for x_end, sgn in ((0.0, -1.0), (lx, 1.0)):
        near = _project(ax, np.array([[x_end + sgn * 0.015 * lx, y_mid, 0.5 * zt]]))[0]
        away = _project(ax, np.array([[x_end + sgn * (0.015 * lx + reach), y_mid, 0.5 * zt]]))[0]
        start, end = (near, away) if tension else (away, near)
        ax.annotate(
            "",
            xy=(float(end[0]), float(end[1])),
            xytext=(float(start[0]), float(start[1])),
            arrowprops={
                "arrowstyle": "-|>,head_length=0.7,head_width=0.4",
                "color": TEXT_2,
                "linewidth": fr.lw(3.0),
            },
            zorder=5,
            annotation_clip=False,
        )
        beyond = away + 0.3 * (away - near)
        ax.text2D(
            float(beyond[0]),
            float(beyond[1]),
            "σ",
            transform=ax.transData,
            fontsize=fr.pt(fr.AXIS_PX),
            color=TEXT_2,
            ha="center",
            va="center",
        )


def _map(
    fig: Figure,
    tl: CaiRampTimeline,
    left: float,
    values: np.ndarray,
    title: str,
    damaged: bool,
    failed: bool,
) -> Axes:
    """A top view of a per-column value (the highest ratio through the
    thickness) on the panel."""
    v = tl.view
    cfg = tl.config
    lx, ly = cfg.panel.Lx_mm, cfg.panel.Ly_mm
    width = min(_MAP_H * lx / ly, float(_MAP_W_MAX))
    height = width * ly / lx
    fr.text_px(fig, left, _MAPS_TOP - 16, title, fr.CAPTION_PX + 2, TEXT_2, ha="left")
    ax = fr.px_axes(fig, left, _MAPS_TOP, width, height)
    ax.pcolormesh(
        v.x_edges,
        v.y_edges,
        np.clip(values, 0.0, 1.0),
        cmap=_CMAP,
        norm=_NORM,
        shading="flat",
        rasterized=False,
    )
    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    ax.set_aspect("equal", adjustable="box", anchor="NW")
    fr.style_axes(ax, grid=False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(BASELINE)
    ax.set_xticks([0, lx])
    ax.set_yticks([0, ly])
    ax.set_xticklabels(["0", f"{lx:g} mm"])
    ax.set_yticklabels(["0", f"{ly:g}"])
    ax.tick_params(labelsize=fr.pt(14), pad=2)
    if damaged:
        for xy in v.zone_full_xy:
            (ln,) = ax.plot(
                xy[:, 0],
                xy[:, 1],
                color=TEXT,
                linewidth=fr.lw(1.8),
                linestyle=(0, (4, 2.5)),
                zorder=3,
            )
            ln.set_path_effects(_HALO)
        cut = ax.axhline(
            v.cut_y_mm,
            color=TEXT_2,
            linewidth=fr.lw(1.6),
            linestyle=(0, (1.5, 2.5)),
            zorder=3,
        )
        cut.set_path_effects(_HALO)
        e_f = tl.fields.damaged.fpf_element
        _, j, i = _unravel(e_f, tl.fields.damaged.grid_shape)
        x0, x1 = v.x_edges[i], v.x_edges[i + 1]
        y0, y1 = v.y_edges[j], v.y_edges[j + 1]
        ax.plot(
            [x0, x1, x1, x0, x0],
            [y0, y0, y1, y1, y0],
            color=_CRITICAL if failed else TEXT,
            linewidth=fr.lw(3.0 if failed else 1.8),
            zorder=4,
        )
    return ax


def _colorbar(fig: Figure) -> None:
    left, top, width, height = _CBAR
    ax = fr.px_axes(fig, left, top, width, height)
    grad = np.linspace(0.0, 1.0, 256)[:, None]
    ax.imshow(grad, aspect="auto", cmap=_CMAP, norm=_NORM, origin="lower", extent=(0, 1, 0, 1))
    ax.set_xticks([])
    ax.yaxis.tick_right()
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1  fails"])
    ax.tick_params(axis="y", labelsize=fr.pt(14), colors=TEXT_2, length=3, width=fr.lw(1.2))
    for s in ax.spines.values():
        s.set_color(BASELINE)
    ax.axhline(1.0, color=_CRITICAL, linewidth=fr.lw(4.0), clip_on=False)
    fr.text_px(fig, left - 4, top - 24, "λ / c_e", fr.CAPTION_PX + 2, TEXT_2, ha="left")


def _maps(fig: Figure, tl: CaiRampTimeline, stress: float, failed: bool) -> None:
    v = tl.view
    m = tl.fields.damaged.multiplier(stress)
    mu = tl.fields.undamaged.multiplier(stress)
    with np.errstate(divide="ignore"):
        dmg = m / v.c_column_damaged
        und = mu / v.c_column_undamaged
    fr.text_px(
        fig,
        fr.MARGIN,
        744,
        "Top views: the highest load ratio through the thickness, at the same applied stress",
        fr.AXIS_PX,
        TEXT,
        weight="bold",
        ha="left",
    )
    _map(fig, tl, _MAP_LEFTS[0], dmg, "damaged panel", True, failed)
    _map(fig, tl, _MAP_LEFTS[1], und, "undamaged panel", False, False)
    _colorbar(fig)
    # Key of the marks on the damaged map.
    rows = [
        ((0, (4, 2.5)), TEXT, "delamination footprint"),
        ((0, (1.5, 2.5)), TEXT_2, "cut of the 3D view"),
        ("solid", _CRITICAL if failed else TEXT, "first element to fail"),
    ]
    key = fr.px_axes(fig, _KEY_X, _MAPS_TOP, 240, 110)
    key.set_xlim(0, 240)
    key.set_ylim(110, 0)
    key.set_axis_off()
    for k, (style, colour, text) in enumerate(rows):
        yy = 10 + 30 * k
        key.plot([0, 34], [yy, yy], color=colour, linewidth=fr.lw(2.2), linestyle=style)
        key.text(44, yy, text, fontsize=fr.pt(fr.CAPTION_PX), color=TEXT_2, va="center")


def _curve(fig: Figure, tl: CaiRampTimeline, stress: float, failed: bool) -> None:
    f = tl.fields
    d, u = f.damaged, f.undamaged
    cfg = tl.config
    tension = cfg.loading != "compression"
    left, top, width, height = _CURVE
    fr.text_px(
        fig,
        _RIGHT_X,
        186,
        "Stress against end shortening" if not tension else "Stress against end extension",
        fr.PANEL_PX,
        TEXT,
        weight="bold",
        ha="left",
    )
    scale = f.pristine_strength_MPa / u.fpf_stress_MPa
    cap = fr.text_px(
        fig,
        _RIGHT_X,
        216,
        f"stress as BVID-FE reports it: fe3d's × {scale:.3f}, so the undamaged panel "
        "fails at the pristine strength",
        fr.CAPTION_PX,
        TEXT_2,
        ha="left",
    )
    fr.fit_width(fig, cap, fr.W - fr.MARGIN - _RIGHT_X)
    ax = fr.px_axes(fig, left, top, width, height)
    fr.style_axes(ax)
    end = tl.end_stress_MPa
    fpf_rep = f.fpf_residual_MPa
    y_top = 1.25 * max(f.residual_strength_MPa, min(fpf_rep, 1.6 * f.residual_strength_MPa))
    sig_top = f.applied_stress(y_top)
    x_top = d.end_shortening_mm(sig_top) * 1.25
    blue, orange = CATEGORICAL[0], CATEGORICAL[1]
    # Undamaged, faint and clipped by the axes.
    su = np.array([0.0, sig_top * 4.0])
    ax.plot(
        [u.end_shortening_mm(s) for s in su],
        [f.reported_stress(s) for s in su],
        color=orange,
        linewidth=fr.lw(2.2),
        zorder=2,
    )
    # Damaged: the ramp, solid to the frame's stress, faint to its end.
    ax.plot(
        [0.0, d.end_shortening_mm(end)],
        [0.0, f.reported_stress(end)],
        color=blue,
        alpha=0.25,
        linewidth=fr.lw(2.5),
        zorder=2,
    )
    ax.plot(
        [0.0, d.end_shortening_mm(stress)],
        [0.0, f.reported_stress(stress)],
        color=blue,
        linewidth=fr.lw(3.5),
        zorder=3,
    )
    for panel, colour, size in ((u, orange, 12), (d, blue, 16)):
        ax.plot(
            [panel.end_shortening_mm(stress)],
            [f.reported_stress(stress)],
            marker="o",
            markersize=fr.pt(size),
            markerfacecolor=colour,
            markeredgecolor=SURFACE,
            markeredgewidth=fr.lw(2.5),
            zorder=5,
            clip_on=False,
        )
    # First-ply failure of the damaged panel, and the growth stress.
    if math.isfinite(fpf_rep) and fpf_rep <= y_top:
        reached = f.governing == "first_ply_failure"
        ax.plot(
            [d.end_shortening_mm(d.fpf_stress_MPa)],
            [fpf_rep],
            marker="X",
            markersize=fr.pt(20),
            markerfacecolor=_CRITICAL if reached else SURFACE,
            markeredgecolor=_CRITICAL,
            markeredgewidth=fr.lw(2.0),
            zorder=6,
            clip_on=False,
        )
        ax.annotate(
            "first-ply failure" if reached else "first-ply failure,\nnot reached",
            xy=(d.end_shortening_mm(d.fpf_stress_MPa), fpf_rep),
            xytext=(14, 0 if reached else -6),
            textcoords="offset points",
            fontsize=fr.pt(fr.CAPTION_PX),
            color=TEXT,
            ha="left",
            va="center" if reached else "top",
            path_effects=_HALO,
            zorder=7,
        )
    if f.governing == "delamination_growth":
        g = f.growth_stress_MPa
        # From the axis to the end of the ramp, where the growth stops it.
        ax.plot(
            [0.0, d.end_shortening_mm(end)],
            [g, g],
            color=TEXT_2,
            linewidth=fr.lw(1.8),
            linestyle=(0, (5, 3)),
            zorder=1,
        )
        ax.text(
            0.02 * x_top,
            g,
            f"delamination growth {g:.1f} MPa",
            fontsize=fr.pt(fr.CAPTION_PX),
            color=TEXT_2,
            ha="left",
            va="bottom",
            path_effects=_HALO,
        )
    ax.set_xlim(0.0, x_top)
    ax.set_ylim(0.0, y_top)
    ax.locator_params(axis="x", nbins=5)
    ax.locator_params(axis="y", nbins=5)
    word = "extension" if tension else "shortening"
    ax.set_xlabel(f"End {word} (mm)", fontsize=fr.pt(fr.AXIS_PX), color=TEXT_2)
    ax.set_ylabel("|σx| (MPa)", fontsize=fr.pt(fr.AXIS_PX), color=TEXT_2)
    # Direct labels: below the damaged line, above the undamaged one.
    s_mid = 0.55 * end
    ax.annotate(
        "damaged",
        xy=(d.end_shortening_mm(s_mid), f.reported_stress(s_mid)),
        xytext=(10, -8),
        textcoords="offset points",
        fontsize=fr.pt(fr.CAPTION_PX),
        color=TEXT_2,
        ha="left",
        va="top",
        path_effects=_HALO,
    )
    s_lab = f.applied_stress(0.93 * y_top)
    ax.annotate(
        "undamaged",
        xy=(u.end_shortening_mm(s_lab), f.reported_stress(s_lab)),
        xytext=(-10, 0),
        textcoords="offset points",
        fontsize=fr.pt(fr.CAPTION_PX),
        color=TEXT_2,
        ha="right",
        va="center",
        path_effects=_HALO,
    )


def ratio_label(r: float) -> str:
    """A load ratio to two decimals, rounded down below failure so that
    "1.00" means an element has failed."""
    if r >= FAILED:
        return f"{r:.2f}"
    return f"{math.floor(r * 100.0) / 100.0:.2f}"


def _gauges(fig: Figure, tl: CaiRampTimeline, stress: float) -> None:
    top = _GAUGE_TOP
    fr.text_px(
        fig,
        _RIGHT_X,
        top,
        "Highest load ratio λ / c_e",
        fr.PANEL_PX,
        TEXT,
        weight="bold",
        ha="left",
    )
    what = "disk-averaged" if tl.ratio == "averaged" else "pointwise (Gauss point)"
    fr.text_px(
        fig,
        _RIGHT_X,
        top + 30,
        f"applied load over the element's own failure load ({what}); 1 = it fails",
        fr.CAPTION_PX,
        TEXT_2,
        ha="left",
    )
    bar_x, bar_w = _RIGHT_X + 130, 420.0
    for k, (label, damaged) in enumerate((("damaged", True), ("undamaged", False))):
        y = top + 70 + 44 * k
        r = tl.max_ratio(stress, damaged)
        fr.text_px(fig, _RIGHT_X, y, label, fr.CAPTION_PX + 2, TEXT_2, ha="left")
        ax = fr.px_axes(fig, bar_x, y - 13, bar_w, 26)
        ax.set_xlim(0.0, 1.1)
        ax.set_ylim(0, 1)
        ax.set_axis_off()
        ax.barh([0.5], [1.1], height=1.0, color=GRID)
        ax.barh([0.5], [min(r, 1.1)], height=1.0, color=_CMAP(_NORM(min(r, 1.0))))
        ax.axvline(1.0, color=_CRITICAL, linewidth=fr.lw(3.0))
        fr.text_px(
            fig,
            fr.W - fr.MARGIN,
            y,
            ratio_label(r),
            fr.CAPTION_PX + 4,
            _CRITICAL if r >= FAILED else TEXT,
            ha="right",
            weight="bold" if r >= FAILED else "normal",
        )
    fr.text_px(
        fig,
        bar_x + bar_w * 1.0 / 1.1,
        top + 70 + 44 + 30,
        "1",
        fr.CAPTION_PX,
        _CRITICAL,
        ha="center",
    )


def _results(fig: Figure, tl: CaiRampTimeline) -> None:
    f = tl.fields
    top = _RESULT_TOP
    load = _load_word(tl.config)
    fr.text_px(
        fig, _RIGHT_X, top, f"fe3d result ({load})", fr.PANEL_PX, TEXT, weight="bold", ha="left"
    )
    if f.governing == "delamination_growth":
        gov = "delamination growth"
        if f.critical_sublaminate is not None:
            gov += f", interface {f.critical_sublaminate}"
    else:
        gov = "first-ply failure"
    rows = [
        ("residual strength", f"{f.residual_strength_MPa:.1f} MPa"),
        ("knockdown", f"{f.knockdown:.3f}"),
        ("pristine strength", f"{f.pristine_strength_MPa:.1f} MPa"),
        ("governed by", gov),
    ]
    if f.governing == "delamination_growth":
        rows.append(("first-ply failure (not reached)", f"{f.fpf_residual_MPa:.1f} MPa"))
    y = top + 38.0
    for label, value in rows:
        fr.text_px(fig, _RIGHT_X, y, label, fr.CAPTION_PX + 2, TEXT_2, ha="left")
        fr.text_px(fig, fr.W - fr.MARGIN, y, value, fr.CAPTION_PX + 2, TEXT, ha="right")
        y += 30


def _event(tl: CaiRampTimeline) -> str:
    f = tl.fields
    e_f = f.damaged.fpf_element
    if f.governing == "delamination_growth":
        iface = f.critical_sublaminate
        where = f" at interface {iface}" if iface is not None else ""
        return (
            f"Delamination growth governs at {f.residual_strength_MPa:.1f} MPa{where} "
            "(closed form, no FE field): the ramp ends here"
        )
    x, y = _element_centre(tl, e_f)
    if not f.damage.delaminations:
        return (
            f"Below damage onset: the panel reaches the pristine strength, "
            f"{f.residual_strength_MPa:.1f} MPa"
        )
    if f.fpf_knockdown >= 1.0:
        return (
            f"The damage does not lower first-ply failure: the panel reaches the pristine "
            f"strength, {f.residual_strength_MPa:.1f} MPa"
        )
    return (
        f"First-ply failure at {f.residual_strength_MPa:.1f} MPa: {_ply_label(tl, e_f)} "
        f"at x = {x:.1f}, y = {y:.1f} mm"
    )


def _scale_note(tl: CaiRampTimeline) -> str:
    v = tl.view
    side = "y above" if v.keep_high else "y below"
    return (
        f"Displacements ×{v.deformation_scale:g}, thickness ×{v.thickness_scale:g}.  "
        f"Cut along the load at y = {v.cut_y_mm:g} mm through the first element to fail; "
        f"the part at {side} the cut is shown, impact face up"
    )


def _mesh_label(tl: CaiRampTimeline) -> str:
    mesh = tl.config.mesh or MeshParams()
    n = tl.fields.damaged.mesh.n_elements
    per = "element" if mesh.elements_per_ply == 1 else "elements"
    return (
        f"fe3d mesh {mesh.in_plane_size_mm:g} mm, {mesh.elements_per_ply} {per} per ply, "
        f"{n:,} elements"
    )


def draw_frame(fig: Figure, timeline: CaiRampTimeline, t: float) -> None:
    """Draw the frame at clip time ``t`` (0 to 1) into ``fig`` (cleared
    first). Pure: no model call."""
    tl = timeline
    cfg = tl.config
    f = tl.fields
    with fr.style_context():
        fig.clear()
        fig.set_facecolor(SURFACE)
        stress = tl.stress_at(t)
        done = tl.schedule.progress(t, "ramp") >= 1.0
        failed = done and f.governing == "first_ply_failure" and tl.max_ratio(stress) >= FAILED
        load = _load_word(cfg)
        parts = [f"{fr.material_name(cfg)}  {fr.layup_label(cfg.layup_deg)}"]
        if cfg.impact is not None:
            parts.append(f"{cfg.impact.energy_J:g} J, {fr.impactor_label(cfg.impact)}")
        else:
            parts.append("damage as supplied")
        parts.append(f"{cfg.panel.Lx_mm:g} × {cfg.panel.Ly_mm:g} mm panel")
        what = "Compression" if load == "CAI" else "Tension"
        fr.draw_title_band(
            fig,
            f"{what} after impact: the fe3d load ramp to failure",
            "  ·  ".join(parts),
            readout=f"{tl.reported_stress_at(t):.1f} MPa",
            readout_sub=f"applied stress |σx|  ·  fe3d tier, {load}",
        )
        caveats = fr.honesty_notes(cfg)
        caveats.append(
            "fe3d: qualitative; linear elastic to first-ply failure, no progressive damage"
        )
        fr.draw_caveats(fig, caveats)
        title = "Damaged panel, cut open through the first element to fail"
        if not f.damage.delaminations:
            title = (
                "Panel below damage onset (no delamination), cut through the first "
                "element to fail"
            )
        fr.text_px(fig, fr.MARGIN, 186, title, fr.PANEL_PX, TEXT, weight="bold", ha="left")
        if done:
            colour = _CRITICAL if failed else STATUS["serious"]
            fr.text_px(fig, fr.MARGIN, 222, "▲", 16, colour, ha="left")
            ev = fr.text_px(
                fig,
                fr.MARGIN + 24,
                222,
                _event(tl),
                fr.CAPTION_PX + 2,
                TEXT,
                weight="bold",
                ha="left",
            )
            fr.fit_width(fig, ev, 1150 - 24)
        _draw_3d(fig, tl, stress, failed)
        note = fr.text_px(fig, fr.MARGIN, 712, _scale_note(tl), fr.CAPTION_PX, TEXT_2, ha="left")
        fr.fit_width(fig, note, 1180 - fr.MARGIN)
        _maps(fig, tl, stress, failed)
        _curve(fig, tl, stress, failed)
        _gauges(fig, tl, stress)
        _results(fig, tl)
        fr.draw_footer(
            fig,
            fr.config_summary(cfg, "fe3d")
            + f"  ·  {_mesh_label(tl)}  ·  one solve per panel, scaled linearly",
        )
