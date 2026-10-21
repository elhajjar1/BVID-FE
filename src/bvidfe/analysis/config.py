"""Analysis configuration dataclasses."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Union

from bvidfe._types import _LOADING_MODES, _TIER_NAMES, LoadingMode, TierName
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.material import OrthotropicMaterial
from bvidfe.damage.state import DamageState
from bvidfe.impact.mapping import ImpactEvent


def _validate_choice(value: object, allowed: Iterable[str], field: str) -> None:
    """Raise ``ValueError`` if *value* is not a member of *allowed*.

    The error message lists the allowed values (sorted for deterministic
    output) and the offending value, prefixed with the *field* name, e.g.
    ``'tier must be one of [...]; got "fe3D"'``.
    """
    allowed_sorted = sorted(allowed)
    if value not in allowed_sorted:
        raise ValueError(f"{field} must be one of {allowed_sorted}; got {value!r}")


@dataclass
class MeshParams:
    """Mesh resolution parameters for the ``tier="fe3d"`` solve.

    Parameters
    ----------
    elements_per_ply : int
        Number of hex elements stacked through the thickness of each ply
        (dimensionless count, ``>= 1``). Total through-thickness elements
        ``nz = n_plies * elements_per_ply``. Default ``1`` (one element
        per ply) is adequate for first-ply-failure screening; increase
        for smoother through-thickness stress gradients at higher cost.
    in_plane_size_mm : float
        Target in-plane element edge length in millimetres (``> 0``). The
        mesh uses ``nx = ceil(Lx_mm / in_plane_size_mm)`` and
        ``ny = ceil(Ly_mm / in_plane_size_mm)`` elements in x and y, so
        smaller values resolve the stress concentration around the damage
        footprint more finely at quadratically higher element count and
        memory. Default ``5.0`` mm. On Streamlit Cloud (1 GB RAM) keep
        this ``>= 5`` mm.
    cohesive_zone_factor : float
        Reserved (``> 0``, default ``1.0``). Intended as a multiplier on
        the cohesive-zone characteristic length once true cohesive
        surfaces land (see README "Limitations"); **currently validated
        but not consumed** by ``build_fe_mesh`` — the present fe3d tier
        uses a component-wise stiffness-reduction model, not cohesive
        elements. Leave at the default.
    """

    elements_per_ply: int = 1
    in_plane_size_mm: float = 5.0
    cohesive_zone_factor: float = 1.0

    def __post_init__(self) -> None:
        if not isinstance(self.elements_per_ply, int) or self.elements_per_ply <= 0:
            raise ValueError(
                f"MeshParams.elements_per_ply must be a positive int "
                f"(got {self.elements_per_ply!r})"
            )
        if not (self.in_plane_size_mm > 0):
            raise ValueError(
                f"MeshParams.in_plane_size_mm must be > 0 " f"(got {self.in_plane_size_mm!r})"
            )
        if not (self.cohesive_zone_factor > 0):
            raise ValueError(
                f"MeshParams.cohesive_zone_factor must be > 0 "
                f"(got {self.cohesive_zone_factor!r})"
            )


@dataclass
class AnalysisConfig:
    """BVID analysis configuration. Provide exactly ONE of `impact` or `damage`.

    ``ply_thickness_mm`` may be either a single positive ``float`` (uniform
    laminate) or a list/tuple of positive floats with one entry per ply
    (matching ``len(layup_deg)``). Per-ply thicknesses let users model
    laminates that mix plies of different fabric weights or prepreg gauges.

    ``impact_panel`` is the panel the impact is applied on, when it is not
    the panel being tested: e.g. a plate impacted on a 60 x 60 mm clamped
    window and then trimmed to an 80 x 50 mm compression coupon. The impact
    is mapped to damage on ``impact_panel`` and the damage is placed on
    ``panel``, the two panels taken as concentric. Only valid with
    ``impact``; defaults to ``panel``.
    """

    material: Union[str, OrthotropicMaterial]
    layup_deg: List[float]
    ply_thickness_mm: Union[float, Sequence[float]]
    panel: PanelGeometry
    loading: LoadingMode = "compression"
    tier: TierName = "empirical"
    impact: Optional[ImpactEvent] = None
    damage: Optional[DamageState] = None
    mesh: Optional[MeshParams] = None
    impact_panel: Optional[PanelGeometry] = None

    def __post_init__(self) -> None:
        _validate_choice(self.tier, _TIER_NAMES, "tier")
        _validate_choice(self.loading, _LOADING_MODES, "loading")
        if (self.impact is None) == (self.damage is None):
            raise ValueError(
                "Provide exactly one of AnalysisConfig.impact or AnalysisConfig.damage"
            )
        if self.impact_panel is not None and self.impact is None:
            raise ValueError(
                "AnalysisConfig.impact_panel only applies to an impact; it cannot be "
                "combined with a given damage state"
            )
        # Validate ply_thickness_mm shape against layup_deg here so callers
        # get an early, descriptive error before Laminate construction.
        if isinstance(self.ply_thickness_mm, (list, tuple)):
            if len(self.ply_thickness_mm) != len(self.layup_deg):
                raise ValueError(
                    f"ply_thickness_mm sequence length "
                    f"({len(self.ply_thickness_mm)}) must equal the number "
                    f"of plies ({len(self.layup_deg)})."
                )
            for i, t in enumerate(self.ply_thickness_mm):
                if not (float(t) > 0):
                    raise ValueError(f"ply_thickness_mm[{i}] must be > 0 (got {t}).")
        else:
            if not (float(self.ply_thickness_mm) > 0):
                raise ValueError(f"ply_thickness_mm must be > 0 (got {self.ply_thickness_mm}).")
        if self.damage is not None:
            # An n-ply laminate has interfaces 0 .. n-2. A larger index used
            # to be silently ignored by the closed-form tiers and to raise
            # IndexError when fe3d placed it on the mesh.
            last = len(self.layup_deg) - 2
            for d in self.damage.delaminations:
                if d.interface_index > last:
                    raise ValueError(
                        f"delamination interface_index {d.interface_index} is outside "
                        f"the {len(self.layup_deg)}-ply layup (interfaces 0-{last})."
                    )
        if self.tier == "fe3d" and self.mesh is None:
            self.mesh = MeshParams()
