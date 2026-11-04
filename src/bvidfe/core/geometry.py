"""Panel and impactor geometry dataclasses for BVID-FE."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from bvidfe._types import _BOUNDARY_KINDS, _IMPACTOR_SHAPES, BoundaryKind, ImpactorShape
from bvidfe._validation import dataclass_kwargs, require_choice, require_positive

# The Literal aliases live in bvidfe._types; re-exported here, where they
# have always been importable from.
__all__ = ["BoundaryKind", "ImpactorShape", "PanelGeometry", "ImpactorGeometry"]


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
        factors (x1.9 clamped, x0.5 free), and through them the
        delamination growth stress. It does not enter the empirical tier's
        strength formulas, only its predicted damage, nor fe3d's
        first-ply-failure solve, a membrane state whose out-of-plane edge
        supports carry no load (``bvidfe.solver.boundary.uniaxial_x_bcs``).

    Raises ``ValueError`` for a non-positive dimension or an unknown
    ``boundary`` (the message lists the valid values).
    """

    Lx_mm: float
    Ly_mm: float
    boundary: BoundaryKind = "simply_supported"

    def __post_init__(self) -> None:
        require_positive(self.Lx_mm, "Lx_mm")
        require_positive(self.Ly_mm, "Ly_mm")
        require_choice(self.boundary, _BOUNDARY_KINDS, "PanelGeometry.boundary")

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "PanelGeometry":
        """Build from a dict of the fields (``dataclasses.asdict`` output);
        ``boundary`` may be left out. ``ValueError`` for an unknown or
        missing key."""
        return cls(**dataclass_kwargs(cls, d, "panel"))


@dataclass
class ImpactorGeometry:
    """Impactor used in a simulated drop/impact event.

    diameter_mm: impactor diameter in millimeters (``> 0``; 16.0 is the ASTM
        D7136 default). It changes the onset energy and scales the damage
        area by ``(DIAMETER_REF_MM / diameter_mm) ** DIAMETER_DPA_EXPONENT``
        (``bvidfe.impact.mapping``); not used for ``conical``.
    shape: 'hemispherical' (standard), 'flat' (cylindrical punch), or
        'conical' (60 deg included angle). The shape changes the onset
        energy (contact stiffness) and, relative to a hemisphere, scales
        the damage area by ``bvidfe.impact.mapping.SHAPE_DPA_FACTOR[shape]``.

    Raises ``ValueError`` for a non-positive diameter or an unknown shape
    (the message lists the valid values).
    """

    diameter_mm: float = 16.0
    shape: ImpactorShape = "hemispherical"

    def __post_init__(self) -> None:
        require_positive(self.diameter_mm, "diameter_mm")
        require_choice(self.shape, _IMPACTOR_SHAPES, "ImpactorGeometry.shape")

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ImpactorGeometry":
        """Build from a dict of the fields, each optional. ``ValueError``
        for an unknown key."""
        return cls(**dataclass_kwargs(cls, d, "impactor"))
