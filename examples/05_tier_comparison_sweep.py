"""Tier-comparison knockdown sweep (the Python-API equivalent of the
File menu -> Compare Tiers GUI action).

Runs the empirical and semi_analytical tiers on the same impact event
across a sweep of energies, overlays their knockdown curves on one
figure, and also produces a CSV with both series.

empirical: Soutis knockdown scales explicitly with DPA (Olsson-predicted).
semi_analytical: Rayleigh-Ritz sublaminate buckling + Soutis envelope;
scales with the largest delamination ellipse.

fe3d is intentionally excluded — per the v0.2.0-dev release notes, its
knockdown is approximately flat vs. energy above the Olsson threshold
and a full fe3d sweep is multi-minute scale. Users who want fe3d spot
checks should run them on individual energies.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from bvidfe.analysis import AnalysisConfig, BvidAnalysis  # noqa: E402
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry  # noqa: E402
from bvidfe.impact.mapping import ImpactEvent  # noqa: E402


OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    base = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
        impact=ImpactEvent(
            energy_J=10.0,  # overridden per point
            impactor=ImpactorGeometry(diameter_mm=16.0),
            mass_kg=5.5,
        ),
        loading="compression",
        tier="empirical",
    )

    energies = list(np.linspace(3.0, 40.0, 12))
    rows: list[dict] = []
    for tier in ("empirical", "semi_analytical"):
        for E in energies:
            new_impact = replace(base.impact, energy_J=float(E))
            cfg = replace(base, impact=new_impact, tier=tier, mesh=None)
            result = BvidAnalysis(cfg).run()
            rows.append(
                {
                    "energy_J": float(E),
                    "tier": tier,
                    "knockdown": result.knockdown,
                    "residual_MPa": result.residual_strength_MPa,
                    "pristine_MPa": result.pristine_strength_MPa,
                    "dpa_mm2": result.dpa_mm2,
                }
            )
            print(
                f"  tier={tier:17s}  E={E:5.1f}J  KD={result.knockdown:.3f}  "
                f"residual={result.residual_strength_MPa:.1f} MPa"
            )

    df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / "tier_comparison_sweep.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nWrote {csv_path} ({len(df)} rows)")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    for tier, tier_df in df.groupby("tier"):
        ax.plot(
            tier_df["energy_J"],
            tier_df["knockdown"],
            marker="o",
            linewidth=2,
            label=tier,
        )
    ax.set_xlabel("Impact energy [J]")
    ax.set_ylabel("Knockdown (residual / pristine)")
    ax.set_title(
        "Knockdown-vs-energy: empirical vs. semi_analytical "
        "(IM7/8552 quasi-iso, 150x100)"
    )
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    png_path = OUTPUT_DIR / "tier_comparison_sweep.png"
    fig.savefig(png_path, dpi=150)
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
