"""Compare all three modeling tiers (empirical / semi_analytical / fe3d) for
the same impact event. Prints a summary table and saves a bar chart."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from bvidfe.analysis import AnalysisConfig, BvidAnalysis  # noqa: E402
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry  # noqa: E402
from bvidfe.impact.mapping import ImpactEvent  # noqa: E402


OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    impact = ImpactEvent(
        energy_J=15.0,
        impactor=ImpactorGeometry(diameter_mm=16.0),
        mass_kg=5.5,
    )

    rows: list[tuple[str, float, float]] = []  # (tier, residual_MPa, knockdown)
    for tier in ("empirical", "semi_analytical", "fe3d"):
        cfg = AnalysisConfig(
            material="IM7/8552",
            layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
            ply_thickness_mm=0.152,
            panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
            impact=impact,
            loading="compression",
            tier=tier,
        )
        r = BvidAnalysis(cfg).run()
        rows.append((tier, r.residual_strength_MPa, r.knockdown))
        print(
            f"{tier:18s}  residual={r.residual_strength_MPa:7.1f} MPa  "
            f"knockdown={r.knockdown:.3f}"
        )

    tiers = [row[0] for row in rows]
    residuals = [row[1] for row in rows]
    knockdowns = [row[2] for row in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    ax1.bar(tiers, residuals, color=["#1f77b4", "#ff7f0e", "#2ca02c"])
    ax1.set_ylabel("Residual strength [MPa]")
    ax1.set_title("CAI residual strength by tier")
    ax1.grid(axis="y", linestyle="--", alpha=0.3)

    ax2.bar(tiers, knockdowns, color=["#1f77b4", "#ff7f0e", "#2ca02c"])
    ax2.set_ylabel("Knockdown (residual / pristine)")
    ax2.set_ylim(0, 1.0)
    ax2.set_title("Knockdown by tier")
    ax2.grid(axis="y", linestyle="--", alpha=0.3)

    fig.suptitle("BVID-FE tier comparison @ 15 J impact, IM7/8552 quasi-iso")
    fig.tight_layout()
    out = OUTPUT_DIR / "tier_comparison.png"
    fig.savefig(out, dpi=150)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
