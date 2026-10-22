"""Tier-comparison knockdown sweep.

Runs the empirical and semi_analytical tiers on the same impact event
across a sweep of energies, overlays their knockdown curves on one
figure, and also produces a CSV with both series.

empirical: Soutis knockdown scales explicitly with DPA (Olsson-predicted).
semi_analytical: the lower of Soutis and the stress at which the weakest
buckled sublaminate grows its delamination (thin-film energy release rate
against a mixed-mode toughness).

On this 8-ply panel the two tiers agree from about 10 J up (the Soutis
bound governs); below that the delamination growth stress is lower and
semi_analytical gives the smaller knockdown (about 0.53 against 0.63 at
3 J).

fe3d is excluded because each energy point takes about 20-30 s at the
default 5 mm mesh (it solves the damaged and the undamaged panel), so the
12-point sweep would take several minutes. Run fe3d spot checks at
individual energies instead (see 02_tier_comparison.py).
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
        "Knockdown-vs-energy: empirical vs. semi_analytical " "(IM7/8552 quasi-iso, 150x100)"
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
