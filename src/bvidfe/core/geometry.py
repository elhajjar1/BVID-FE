"""Panel and impactor geometry dataclasses for BVID-FE."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

BoundaryKind = Literal["clamped", "simply_supported", "free"]
ImpactorShape = Literal["hemispherical", "flat", "conical"]


@dataclass
class PanelGeometry:
    """In-plane rectangular panel with boundary condition.

    Lx_mm, Ly_mm: in-plane dimensions in millimeters (``> 0``). x is the
        load direction; the panel spans ``0 <= x <= Lx_mm``,
        ``0 <= y <= Ly_mm``.
    boundary: edge condition, ``"simply_supported"`` (default),
        ``"clamped"`` or ``"free"``. It sets the plate bending stiffness in
        the impact-onset model (x2.5 clamped, x0.4 free, with a
        ``UserWarning`` for free), the sublaminate and panel buckling
        factors (x1.9 clamped, x0.5 free), and the fe3d out-of-plane edge
        restraint. It does not enter the empirical tier's strength
        formulas, only its predicted damage.

    Raises ``ValueError`` for a non-positive dimension or an unknown
    ``boundary``.
    """

    Lx_mm: float
    Ly_mm: float
    boundary: BoundaryKind = "simply_supported"

    def __post_init__(self) -> None:
        if self.Lx_mm <= 0:
            raise ValueError(f"Lx_mm must be > 0 (got {self.Lx_mm})")
        if self.Ly_mm <= 0:
            raise ValueError(f"Ly_mm must be > 0 (got {self.Ly_mm})")
        if self.boundary not in ("clamped", "simply_supported", "free"):
            raise ValueError(f"boundary {self.boundary!r} not recognized")


@dataclass
class ImpactorGeometry:
    """Impactor used in a simulated drop/impact event.

    diameter_mm: impactor diameter in millimeters (``> 0``; 16.0 is the ASTM
        D7136 default). It changes the onset energy and scales the damage
        area by ``(16 / diameter_mm) ** 0.2``; not used for ``conical``.
    shape: 'hemispherical' (standard), 'flat' (cylindrical punch), or
        'conical' (60 deg included angle). Relative to a hemisphere the
        damage area is x0.2 for flat and x1.3 for conical.

    Raises ``ValueError`` for a non-positive diameter or an unknown shape.
    """

    diameter_mm: float = 16.0
    shape: ImpactorShape = "hemispherical"

    def __post_init__(self) -> None:
        if self.diameter_mm <= 0:
            raise ValueError(f"diameter_mm must be > 0 (got {self.diameter_mm})")
        if self.shape not in ("hemispherical", "flat", "conical"):
            raise ValueError(f"shape {self.shape!r} not recognized")
