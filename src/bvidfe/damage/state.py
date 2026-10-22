"""BVID damage state: per-interface elliptical delaminations + dent + fiber-break core."""

from __future__ import annotations

import math
import numbers
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

from bvidfe._validation import require_finite, require_non_negative, require_positive, require_xy


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
    """A single delamination at a specific ply interface, in the panel XY frame."""

    interface_index: int
    centroid_mm: Tuple[float, float]
    major_mm: float
    minor_mm: float
    orientation_deg: float

    def __post_init__(self) -> None:
        if isinstance(self.interface_index, bool) or not isinstance(
            self.interface_index, numbers.Integral
        ):
            raise ValueError(
                f"interface_index must be an int (got {type(self.interface_index).__name__}: "
                f"{self.interface_index!r})"
            )
        if self.interface_index < 0:
            raise ValueError("interface_index must be >= 0")
        self.centroid_mm = require_xy(self.centroid_mm, "centroid_mm")
        try:
            require_positive(self.major_mm, "major_mm")
            require_positive(self.minor_mm, "minor_mm")
        except ValueError as exc:
            raise ValueError(f"ellipse semi-axes must be positive and finite: {exc}") from exc
        require_finite(self.orientation_deg, "orientation_deg")

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
    """Full BVID damage description handed between impact mapping and solvers."""

    delaminations: List[DelaminationEllipse] = field(default_factory=list)
    dent_depth_mm: float = 0.0
    fiber_break_radius_mm: float = 0.0

    def __post_init__(self) -> None:
        require_non_negative(self.dent_depth_mm, "dent_depth_mm")
        require_non_negative(self.fiber_break_radius_mm, "fiber_break_radius_mm")
        for i, e in enumerate(self.delaminations):
            if not isinstance(e, DelaminationEllipse):
                raise ValueError(
                    f"delaminations[{i}] must be a DelaminationEllipse " f"(got {type(e).__name__})"
                )

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
