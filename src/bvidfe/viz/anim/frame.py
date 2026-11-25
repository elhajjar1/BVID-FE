"""The 1080p frame template shared by every scene.

A frame is laid out in pixels on a 1920 x 1080 canvas (19.2 x 10.8 in at
100 dpi); another 16:9 size renders the same layout at another dpi. From
the top:

* the title band (y 0-118): title, subtitle, and a readout on the right;
* the caveat row (y 124-156): honest labels, e.g. "experimental, not
  validated" for a biaxial load;
* the content area (y 166-1000), laid out by the scene;
* the footer (y 1006-1080): the citation with the version, then a summary
  of the configuration.

Type scale at 1080p: title 34 px, panel titles 24 px, axis labels 20 px,
ticks and captions 16 px. Fonts are matplotlib's bundled DejaVu Sans, so
frames do not depend on the host's fonts.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterator, Optional, Sequence, Tuple, Union, cast

import matplotlib as mpl
import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from bvidfe.citation import frame_credit, short_credit
from bvidfe.viz.anim.spec import DESIGN_HEIGHT_PX, DESIGN_WIDTH_PX, FIGSIZE_IN, AnimationSpec
from bvidfe.viz.style import BASELINE, GRID, MUTED, STATUS, SURFACE, TEXT, TEXT_2

if TYPE_CHECKING:  # pragma: no cover
    from bvidfe.analysis.config import AnalysisConfig
    from bvidfe.impact.mapping import ImpactEvent

W, H = DESIGN_WIDTH_PX, DESIGN_HEIGHT_PX
MARGIN = 48  # left and right page margin, px

# Type scale, px at 1080p.
TITLE_PX = 34
SUBTITLE_PX = 20
PANEL_PX = 24
AXIS_PX = 20
TICK_PX = 16
CAPTION_PX = 16
READOUT_PX = 34

# Vertical bands, px from the top.
TITLE_BAND_BOTTOM = 118
CAVEAT_Y = 140
CONTENT_TOP = 166
CONTENT_BOTTOM = 1000
FOOTER_TOP = 1006

CITATION_GID = "bvidfe-citation"
COMPACT_BELOW_PX = 1279  # narrower frames get the compact footer
SUMMARY_GID = "bvidfe-summary"
TITLE_GID = "bvidfe-title"
CAVEAT_GID = "bvidfe-caveat"


def pt(px: float) -> float:
    """Font size in points for ``px`` pixels at the 100 dpi design canvas."""
    return px * 72.0 / 100.0


def lw(px: float) -> float:
    """Line width in points for ``px`` pixels."""
    return px * 72.0 / 100.0


_RC = {
    "font.family": "DejaVu Sans",
    "font.size": pt(TICK_PX),
    "mathtext.fontset": "dejavusans",
    "text.color": TEXT,
    "axes.edgecolor": BASELINE,
    "axes.labelcolor": TEXT_2,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "xtick.color": TEXT_2,
    "ytick.color": TEXT_2,
    "xtick.labelsize": pt(TICK_PX),
    "ytick.labelsize": pt(TICK_PX),
    "axes.labelsize": pt(AXIS_PX),
    "axes.linewidth": lw(1.5),
    "xtick.major.width": lw(1.5),
    "ytick.major.width": lw(1.5),
    "xtick.major.size": 4.0,
    "ytick.major.size": 4.0,
    "axes.unicode_minus": True,
    "lines.solid_capstyle": "round",
    "hatch.linewidth": lw(1.5),
    "savefig.bbox": None,
    "path.simplify": False,
}


@contextmanager
def style_context() -> Iterator[None]:
    """matplotlib rcParams for frames, applied only inside the block."""
    with mpl.rc_context(cast(Dict[Any, Any], _RC)):
        yield


def new_figure(spec: Optional[AnimationSpec] = None) -> Figure:
    """A frame figure (no pyplot): 19.2 x 10.8 in at the dpi of ``spec``
    (1920 x 1080 px by default), with an Agg canvas."""
    spec = spec or AnimationSpec()
    fig = Figure(figsize=FIGSIZE_IN, dpi=spec.dpi, facecolor=SURFACE)
    FigureCanvasAgg(fig)
    return fig


def px_rect(
    left: float, top: float, width: float, height: float
) -> Tuple[float, float, float, float]:
    """A pixel box (from the top left of the canvas) as figure fractions."""
    return (left / W, 1.0 - (top + height) / H, width / W, height / H)


def px_axes(fig: Figure, left: float, top: float, width: float, height: float, **kw: Any) -> Axes:
    """Axes occupying the pixel box ``left, top, width, height``."""
    return fig.add_axes(px_rect(left, top, width, height), **kw)


def text_px(
    fig: Figure,
    x: float,
    y: float,
    s: str,
    size_px: float = CAPTION_PX,
    color: str = TEXT,
    **kw: Any,
) -> Text:
    """Figure text at pixel ``(x, y)`` from the top left."""
    kw.setdefault("va", "center")
    return fig.text(x / W, 1.0 - y / H, s, fontsize=pt(size_px), color=color, **kw)


def hline_px(fig: Figure, y: float, x0: float = MARGIN, x1: float = W - MARGIN, **kw: Any) -> None:
    """A horizontal hairline across the canvas at pixel row ``y``."""
    kw.setdefault("color", GRID)
    kw.setdefault("linewidth", lw(1.5))
    fig.add_artist(
        mpl.lines.Line2D(
            [x0 / W, x1 / W], [1.0 - y / H, 1.0 - y / H], transform=fig.transFigure, **kw
        )
    )


def text_width_px(fig: Figure, t: Text) -> float:
    """Drawn width of ``t`` in design pixels (1920-wide canvas); 0 when the
    figure's canvas cannot measure text (not an Agg canvas)."""
    get_renderer = getattr(fig.canvas, "get_renderer", None)
    if get_renderer is None:
        return 0.0
    return float(t.get_window_extent(get_renderer()).width) * 100.0 / fig.dpi


def fit_width(fig: Figure, t: Text, max_px: float, min_px: float = 12.0) -> None:
    """Shrink ``t`` to fit ``max_px`` (design pixels), down to ``min_px``
    type; past that, shorten it with an ellipsis."""
    width = text_width_px(fig, t)
    if width <= max_px:
        return
    size = float(t.get_fontsize())
    new = max(pt(min_px), size * max_px / width * 0.99)
    t.set_fontsize(new)
    width = text_width_px(fig, t)
    s = t.get_text()
    while width > max_px and len(s) > 4:
        s = s[: max(4, int(len(s) * max_px / width) - 2)].rstrip()
        t.set_text(s + "\u2026")
        width = text_width_px(fig, t)


def draw_title_band(
    fig: Figure,
    title: str,
    subtitle: str,
    readout: str = "",
    readout_sub: str = "",
) -> None:
    """Title and subtitle on the left; a large readout (e.g. the current
    energy) and its caption (e.g. the model tier) on the right."""
    text_px(fig, MARGIN, 44, title, TITLE_PX, TEXT, weight="bold", ha="left", gid=TITLE_GID)
    sub = text_px(fig, MARGIN, 90, subtitle, SUBTITLE_PX, TEXT_2, ha="left")
    if readout:
        text_px(fig, W - MARGIN, 44, readout, READOUT_PX, TEXT, weight="bold", ha="right")
    room = float(W - 2 * MARGIN)
    if readout_sub:
        right = text_px(fig, W - MARGIN, 90, readout_sub, SUBTITLE_PX, TEXT_2, ha="right")
        room -= text_width_px(fig, right) + 48
    fit_width(fig, sub, room)
    hline_px(fig, TITLE_BAND_BOTTOM)


def draw_caveats(fig: Figure, caveats: Sequence[str]) -> None:
    """Honest labels on the caveat row, each after a status-coloured mark
    (the text stays in ink: the mark is the icon). When they do not fit
    the width of the frame, every label is set smaller (down to 13 px)."""
    size = 18.0
    room = float(W - 2 * MARGIN)
    for attempt in range(3):
        x = float(MARGIN)
        texts = []
        for k, note in enumerate(caveats):
            colour = STATUS["critical"] if k == 0 else STATUS["serious"]
            mark = text_px(fig, x, CAVEAT_Y, "▲", size * 16 / 18, colour, ha="left")
            t = text_px(
                fig, x + 24, CAVEAT_Y, note, size, TEXT, weight="bold", ha="left", gid=CAVEAT_GID
            )
            texts += [mark, t]
            x += 24 + (text_width_px(fig, t) or 0.6 * size * len(note)) + 48 * size / 18
        used = x - 48 * size / 18 - MARGIN
        if used <= room or size <= 13.0 or attempt == 2:
            return
        for t in texts:
            t.remove()
        size = max(13.0, size * room / used * 0.99)


def draw_footer(fig: Figure, summary: str) -> None:
    """The citation footer, on every frame, then the configuration summary.

    A frame narrower than 1280 px (a GIF preview) gets the short credit in
    larger type instead, without the summary, so it stays readable."""
    hline_px(fig, FOOTER_TOP)
    if fig.dpi * FIGSIZE_IN[0] < COMPACT_BELOW_PX:
        t = text_px(fig, MARGIN, 1043, short_credit(), 30, TEXT_2, ha="left", gid=CITATION_GID)
        fit_width(fig, t, W - 2 * MARGIN)
        return
    text_px(fig, MARGIN, 1030, frame_credit(), CAPTION_PX, TEXT_2, ha="left", gid=CITATION_GID)
    t = text_px(fig, MARGIN, 1058, summary, CAPTION_PX, MUTED, ha="left", gid=SUMMARY_GID)
    fit_width(fig, t, W - 2 * MARGIN)


def style_axes(ax: Axes, grid: bool = True) -> None:
    """Recessive chrome: left and bottom baselines, hairline grid."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(colors=TEXT_2, labelsize=pt(TICK_PX), length=4, width=lw(1.5))
    if grid:
        ax.grid(True, color=GRID, linewidth=lw(1.2))
        ax.set_axisbelow(True)


def render_rgba(fig: Figure) -> np.ndarray:
    """Draw ``fig`` and return its pixels, an ``(H, W, 4)`` uint8 array."""
    with style_context():
        canvas = fig.canvas
        assert isinstance(canvas, FigureCanvasAgg)
        canvas.draw()
        return np.asarray(canvas.buffer_rgba()).copy()


def save_png(rgba: np.ndarray, path: Union[str, Path]) -> Path:
    """Write an RGBA frame to ``path`` as PNG (exact pixels, no resampling)."""
    from PIL import Image

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgba, mode="RGBA").save(p)
    return p


def find_text(fig: Figure, gid: str) -> list[str]:
    """The strings of the figure texts with ``gid`` (for tests and checks)."""
    return [t.get_text() for t in fig.texts if t.get_gid() == gid]


# ---- Labels from a configuration ----------------------------------------------


def _angle(a: float) -> str:
    return f"{a:g}"


def layup_label(layup: Sequence[float]) -> str:
    """Compact stacking notation: ``[0/45/-45/90]4``, ``[0/45/-45/90]s``,
    ``[0/45/-45/90]2s``, or the plies listed out (shortened past 12).

    >>> layup_label([0, 45, -45, 90] * 4)
    '[0/45/-45/90]4'
    >>> layup_label([0, 45, -45, 90, 90, -45, 45, 0])
    '[0/45/-45/90]s'
    """
    plies = [float(a) for a in layup]
    n = len(plies)

    def repeat(seq: list[float]) -> str:
        m = len(seq)
        for size in range(1, m // 2 + 1):
            if m % size == 0 and seq == seq[:size] * (m // size):
                unit = "/".join(_angle(a) for a in seq[:size])
                return f"[{unit}]{m // size}"
        return "[" + "/".join(_angle(a) for a in seq) + "]"

    if n >= 2 and n % 2 == 0 and plies[: n // 2] == plies[n // 2 :][::-1]:
        return repeat(plies[: n // 2]) + "s"
    label = repeat(plies)
    if label.endswith("]") and n > 12:
        head = "/".join(_angle(a) for a in plies[:4])
        tail = "/".join(_angle(a) for a in plies[-2:])
        return f"[{head}/.../{tail}] ({n} plies)"
    return label


def material_name(cfg: "AnalysisConfig") -> str:
    m = cfg.material
    return m if isinstance(m, str) else str(getattr(m, "name", "custom material"))


def impactor_label(impact: "ImpactEvent") -> str:
    """``16 mm hemispherical impactor, 5.5 kg``."""
    imp = impact.impactor
    if imp.shape == "conical":
        tip = "conical impactor"
    else:
        tip = f"{imp.diameter_mm:g} mm {imp.shape} impactor"
    return f"{tip}, {impact.mass_kg:g} kg"


def panel_label(cfg: "AnalysisConfig") -> str:
    p = cfg.panel
    label = f"{p.Lx_mm:g} × {p.Ly_mm:g} mm panel ({p.boundary.replace('_', ' ')})"
    if cfg.impact_panel is not None:
        q = cfg.impact_panel
        label += f", impacted on a {q.Lx_mm:g} × {q.Ly_mm:g} mm window"
    return label


def thickness_label(cfg: "AnalysisConfig") -> str:
    from bvidfe.core.laminate import _normalise_ply_thicknesses

    ts = _normalise_ply_thicknesses(cfg.ply_thickness_mm, len(cfg.layup_deg))
    if all(abs(t - ts[0]) < 1e-12 for t in ts):
        return f"{len(ts)} plies × {ts[0]:g} mm"
    return f"{len(ts)} plies, {sum(ts):.3g} mm"


def config_summary(cfg: "AnalysisConfig", tier: Optional[str] = None) -> str:
    """One line describing the configuration, for the footer."""
    parts = [
        material_name(cfg),
        layup_label(cfg.layup_deg),
        thickness_label(cfg),
        panel_label(cfg),
    ]
    if cfg.impact is not None:
        parts.append(impactor_label(cfg.impact))
    parts.append("CAI" if cfg.loading == "compression" else "TAI")
    if cfg.transverse_stress_ratio != 0:
        parts.append(f"σy/σx = {cfg.transverse_stress_ratio:g}")
    if tier:
        parts.append(f"{tier} tier")
    return "  ·  ".join(parts)


def honesty_notes(cfg: "AnalysisConfig") -> list[str]:
    """Caveats every frame of ``cfg`` must carry."""
    notes = []
    if cfg.transverse_stress_ratio != 0:
        notes.append(
            f"Biaxial load σy/σx = {cfg.transverse_stress_ratio:g}: " "experimental, not validated"
        )
    if cfg.loading == "tension":
        notes.append("TAI: not validated (one measured condition)")
    return notes
