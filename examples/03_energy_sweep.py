"""Sweep impact energy 5..40 J and plot the knockdown curve."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from bvidfe.analysis import AnalysisConfig  # noqa: E402
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry  # noqa: E402
from bvidfe.impact.mapping import ImpactEvent  # noqa: E402
from bvidfe.sweep.parametric_sweep import sweep_energies  # noqa: E402


OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
        impact=ImpactEvent(
            energy_J=10.0,  # placeholder; sweep_energies overrides
            impactor=ImpactorGeometry(diameter_mm=16.0),
            mass_kg=5.5,
        ),
        loading="compression",
        tier="empirical",
    )
    energies = [5, 10, 15, 20, 25, 30, 35, 40]

    csv_path = OUTPUT_DIR / "energy_sweep.csv"
    df = sweep_energies(cfg, energies_J=energies, csv_path=csv_path)
    print(df.to_string(index=False))
    print(f"\nWrote {csv_path}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df["energy_J"], df["knockdown"], marker="o", linewidth=2)
    ax.set_xlabel("Impact energy [J]")
    ax.set_ylabel("Knockdown (residual / pristine)")
    ax.set_title("BVID-FE empirical knockdown curve, IM7/8552 quasi-iso, 150x100 mm")
    ax.set_ylim(0, 1.0)
    ax.grid(True, linestyle="--", alpha=0.3)

    out = OUTPUT_DIR / "knockdown_curve.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
