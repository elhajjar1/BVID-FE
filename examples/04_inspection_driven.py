"""Load a C-scan JSON and run the damage-driven (inspection) path."""

from __future__ import annotations

from pathlib import Path

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry
from bvidfe.damage.io import load_cscan_json


def main() -> None:
    cscan = Path(__file__).parent / "sample_cscan.json"
    damage = load_cscan_json(cscan)
    print(f"Loaded {len(damage.delaminations)} delaminations from {cscan.name}")
    print(f"  projected damage area: {damage.projected_damage_area_mm2:.0f} mm^2")
    print(f"  dent depth: {damage.dent_depth_mm:.3f} mm")

    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
        damage=damage,  # <-- inspection-driven entry point
        loading="compression",
        tier="semi_analytical",
    )
    result = BvidAnalysis(cfg).run()
    print()
    print(result.summary())


if __name__ == "__main__":
    main()
