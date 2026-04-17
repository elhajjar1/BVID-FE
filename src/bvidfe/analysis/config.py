"""Analysis configuration dataclasses."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Literal, Optional, Union

from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.material import OrthotropicMaterial
from bvidfe.damage.state import DamageState
from bvidfe.impact.mapping import ImpactEvent


@dataclass
class MeshParams:
    """Mesh resolution parameters for the 3D FE tier."""

    elements_per_ply: int = 1
    in_plane_size_mm: float = 5.0
    cohesive_zone_factor: float = 1.0


@dataclass
class AnalysisConfig:
    """BVID analysis configuration. Provide exactly ONE of `impact` or `damage`."""

    material: Union[str, OrthotropicMaterial]
    layup_deg: List[float]
    ply_thickness_mm: float
    panel: PanelGeometry
    loading: Literal["compression", "tension"] = "compression"
    tier: Literal["empirical", "semi_analytical", "fe3d"] = "empirical"
    impact: Optional[ImpactEvent] = None
    damage: Optional[DamageState] = None
    mesh: Optional[MeshParams] = None

    def __post_init__(self) -> None:
        if (self.impact is None) == (self.damage is None):
            raise ValueError(
                "Provide exactly one of AnalysisConfig.impact or AnalysisConfig.damage"
            )
        if self.tier == "fe3d" and self.mesh is None:
            self.mesh = MeshParams()
