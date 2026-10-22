"""Panel and impactor geometry dataclasses for BVID-FE."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from bvidfe._validation import require_positive

BoundaryKind = Literal["clamped", "simply_supported", "free"]
ImpactorShape = Literal["hemispherical", "flat", "conical"]


@dataclass
class PanelGeometry:
    """In-plane rectangular panel with boundary condition.

    Lx_mm, Ly_mm: in-plane dimensions in millimeters.
    boundary: edge condition applied by downstream solvers.
    """

    Lx_mm: float
    Ly_mm: float
    boundary: BoundaryKind = "simply_supported"

    def __post_init__(self) -> None:
        require_positive(self.Lx_mm, "Lx_mm")
        require_positive(self.Ly_mm, "Ly_mm")
        if self.boundary not in ("clamped", "simply_supported", "free"):
            raise ValueError(f"boundary {self.boundary!r} not recognized")


@dataclass
class ImpactorGeometry:
    """Impactor used in a simulated drop/impact event.

    diameter_mm: impactor diameter in millimeters (16.0 is ASTM D7136 default).
    shape: 'hemispherical' (standard), 'flat' (cylindrical punch), or 'conical'.
    """

    diameter_mm: float = 16.0
    shape: ImpactorShape = "hemispherical"

    def __post_init__(self) -> None:
        require_positive(self.diameter_mm, "diameter_mm")
        if self.shape not in ("hemispherical", "flat", "conical"):
            raise ValueError(f"shape {self.shape!r} not recognized")
