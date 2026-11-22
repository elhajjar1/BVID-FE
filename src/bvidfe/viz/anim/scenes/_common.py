"""Helpers shared by the scenes: the impact point, the onset energy and the
region of the panel a view frames."""

from __future__ import annotations

from typing import Iterable, Optional, Tuple

from bvidfe.analysis.bvid import _resolve_material
from bvidfe.analysis.config import AnalysisConfig
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.damage.state import DelaminationEllipse
from bvidfe.impact.olsson import onset_energy

Box = Tuple[float, float, float, float]  # x0, y0, x1, y1 in mm

# A view zooms onto the damage when the damage (with a margin) spans less
# than this fraction of the panel in both directions.
ZOOM_FRACTION = 0.45


def laminate(cfg: AnalysisConfig) -> Laminate:
    return Laminate(_resolve_material(cfg.material), cfg.layup_deg, cfg.ply_thickness_mm)


def onset_energy_J(cfg: AnalysisConfig) -> float:
    """The damage onset energy of ``cfg``'s impact, on the impacted panel
    (``impact_panel`` or ``panel``), as the impact mapping computes it."""
    if cfg.impact is None:
        raise ValueError("an impact-driven config (AnalysisConfig.impact) is needed")
    pan = cfg.impact_panel or cfg.panel
    loc = cfg.impact.location_xy_mm or (pan.Lx_mm / 2, pan.Ly_mm / 2)
    return float(onset_energy(laminate(cfg), pan, cfg.impact.impactor, location_xy_mm=loc))


def impact_point(cfg: AnalysisConfig) -> Tuple[float, float]:
    """The impact point in the frame of ``cfg.panel`` (the impacted panel is
    concentric with it)."""
    assert cfg.impact is not None
    pan = cfg.impact_panel or cfg.panel
    x, y = cfg.impact.location_xy_mm or (pan.Lx_mm / 2, pan.Ly_mm / 2)
    return x + (cfg.panel.Lx_mm - pan.Lx_mm) / 2, y + (cfg.panel.Ly_mm - pan.Ly_mm) / 2


def damage_region(
    panel: PanelGeometry,
    ellipses: Iterable[DelaminationEllipse],
    centre: Optional[Tuple[float, float]],
    aspect: float,
    min_half_mm: float = 0.0,
) -> Tuple[Box, bool]:
    """``(box, zoomed)``: the region a view of ``panel`` frames. The whole
    panel, unless the ellipses with a margin span less than
    ``ZOOM_FRACTION`` of it both ways: then a box of width / height
    ``aspect`` around ``centre`` (default the damage's centre) holding them,
    at least ``min_half_mm`` from the centre each way.
    """
    whole: Box = (0.0, 0.0, panel.Lx_mm, panel.Ly_mm)
    ells = list(ellipses)
    if not ells:
        return whole, False
    xs, ys = [], []
    for e in ells:
        ex, ey = e.half_extents_mm()
        xs += [e.centroid_mm[0] - ex, e.centroid_mm[0] + ex]
        ys += [e.centroid_mm[1] - ey, e.centroid_mm[1] + ey]
    if centre is None:
        centre = (0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)))
    cx, cy = centre
    hx = max(abs(v - cx) for v in xs) * 1.3 + 2.0
    hy = max(abs(v - cy) for v in ys) * 1.3 + 2.0
    hx, hy = max(hx, min_half_mm), max(hy, min_half_mm)
    if hx / hy < aspect:
        hx = hy * aspect
    else:
        hy = hx / aspect
    if 2 * hx >= ZOOM_FRACTION * panel.Lx_mm or 2 * hy >= ZOOM_FRACTION * panel.Ly_mm:
        return whole, False
    return (cx - hx, cy - hy, cx + hx, cy + hy), True
