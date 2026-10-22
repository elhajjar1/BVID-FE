"""BVID damage state: per-interface elliptical delaminations + dent + fiber-break core."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union


def _ellipse_polygon(e: "DelaminationEllipse", n_pts: int = 72) -> Polygon:
    """Approximate the ellipse footprint as a polygon in the panel frame."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_pts, endpoint=False)
    c = math.cos(math.radians(e.orientation_deg))
    s = math.sin(math.radians(e.orientation_deg))
    x = e.major_mm * np.cos(theta)
    y = e.minor_mm * np.sin(theta)
    xr = c * x - s * y + e.centroid_mm[0]
    yr = s * x + c * y + e.centroid_mm[1]
    return Polygon(zip(xr.tolist(), yr.tolist()))


@dataclass
class DelaminationEllipse:
    """A single delamination at a specific ply interface, in the panel XY frame.

    Attributes
    ----------
    interface_index : int
        Interface between plies ``interface_index`` and ``interface_index + 1``
        (0-based, in ``layup_deg`` order); ``>= 0``. ``AnalysisConfig``
        rejects an index past the last interface of its layup.
    centroid_mm : tuple of float
        Ellipse centre ``(x, y)`` in mm, in the frame of the panel analysed
        (``0 <= x <= Lx_mm``, x the load direction).
    major_mm, minor_mm : float
        Semi-axes in mm (``> 0``); ``area_mm2 = pi * major_mm * minor_mm``.
    orientation_deg : float
        Angle of the ``major_mm`` axis from the panel x axis, counterclockwise.

    Raises ``ValueError`` for a non-positive semi-axis or a negative
    ``interface_index``.
    """

    interface_index: int
    centroid_mm: Tuple[float, float]
    major_mm: float
    minor_mm: float
    orientation_deg: float

    def __post_init__(self) -> None:
        if self.major_mm <= 0 or self.minor_mm <= 0:
            raise ValueError("ellipse semi-axes must be positive")
        if self.interface_index < 0:
            raise ValueError("interface_index must be >= 0")

    @property
    def area_mm2(self) -> float:
        return math.pi * self.major_mm * self.minor_mm

    def __repr__(self) -> str:
        return (
            f"DelaminationEllipse(iface={self.interface_index}, "
            f"a={self.major_mm:.2g}mm, b={self.minor_mm:.2g}mm, "
            f"theta={self.orientation_deg:.0f}deg)"
        )


@dataclass
class DamageState:
    """Full BVID damage description handed between impact mapping and solvers.

    Attributes
    ----------
    delaminations : list of DelaminationEllipse
        The delaminations; empty for an undamaged panel.
    dent_depth_mm : float
        Residual dent depth on the impact face (mm). Stored and reported;
        no tier's strength uses it.
    fiber_break_radius_mm : float
        Radius (mm) of the fibre-break core; fe3d softens the elements
        within it of any delamination centroid. ``0`` disables it.

    ``projected_damage_area_mm2`` (the DPA) is the area of the union of the
    footprints, each drawn as a 72-sided polygon (about 0.1% below the
    exact ellipse area).

    Examples
    --------
    >>> d = DamageState([DelaminationEllipse(3, (75.0, 50.0), 20.0, 10.0, 45.0)])
    >>> round(d.projected_damage_area_mm2)
    628
    """

    delaminations: List[DelaminationEllipse] = field(default_factory=list)
    dent_depth_mm: float = 0.0
    fiber_break_radius_mm: float = 0.0

    @property
    def projected_damage_area_mm2(self) -> float:
        """Union of all delamination footprints (plan view, mm^2)."""
        if not self.delaminations:
            return 0.0
        polys = [_ellipse_polygon(e) for e in self.delaminations]
        return float(unary_union(polys).area)

    @property
    def per_interface_area(self) -> Dict[int, float]:
        """Sum of ellipse areas at each interface (no per-interface union)."""
        out: Dict[int, float] = {}
        for e in self.delaminations:
            out[e.interface_index] = out.get(e.interface_index, 0.0) + e.area_mm2
        return out

    def __repr__(self) -> str:
        return (
            f"DamageState(n_delam={len(self.delaminations)}, "
            f"dent={self.dent_depth_mm:.3g}mm, "
            f"fbr={self.fiber_break_radius_mm:.3g}mm)"
        )
