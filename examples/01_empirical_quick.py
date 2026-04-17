"""Minimal BVID-FE example: 30-J impact on IM7/8552 quasi-iso, empirical tier."""

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


def main() -> None:
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
        impact=ImpactEvent(
            energy_J=30.0,
            impactor=ImpactorGeometry(diameter_mm=16.0),
            mass_kg=5.5,
        ),
        loading="compression",
        tier="empirical",
    )
    result = BvidAnalysis(cfg).run()
    print(result.summary())


if __name__ == "__main__":
    main()
