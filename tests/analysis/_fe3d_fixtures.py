"""Shared fe3d test inputs: the 8-ply quasi-isotropic IM7/8552 layup and an
impact on a 150 x 100 mm panel of it at a coarse 10 mm mesh."""

from __future__ import annotations

from bvidfe import AnalysisConfig, MeshParams
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent

QI8 = [0, 45, -45, 90, 90, -45, 45, 0]


def impact_cfg(loading: str, energy_J: float) -> AnalysisConfig:
    """fe3d config for a ``energy_J`` impact (16 mm hemisphere, 5.5 kg) on an
    8-ply QI IM7/8552 150 x 100 mm panel, one element per ply, 10 mm in plane."""
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        loading=loading,
        tier="fe3d",
        impact=ImpactEvent(energy_J, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0),
    )
