"""2D matplotlib plots for BVID-FE: damage map, knockdown curve, tier comparison."""

from __future__ import annotations

from typing import Dict, Sequence

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from bvidfe.core.geometry import PanelGeometry
from bvidfe.damage.state import DamageState
from bvidfe.viz.style import COLORS, ELLIPSE_CMAP, FIGSIZE_LANDSCAPE, FIGSIZE_SQUARE


def plot_damage_map(damage: DamageState, panel: PanelGeometry, title: str | None = None):
    """Top-down plan view of ellipse delamination footprints + panel outline.

    Ellipses are color-coded by interface index. Returns the matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    ax.set_xlim(0, panel.Lx_mm)
    ax.set_ylim(0, panel.Ly_mm)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(title or f"BVID damage map (dent {damage.dent_depth_mm:.2f} mm)")

    # Panel outline
    ax.add_patch(
        mpatches.Rectangle(
            (0, 0),
            panel.Lx_mm,
            panel.Ly_mm,
            fill=False,
            edgecolor="black",
            linewidth=1.5,
        )
    )

    if damage.delaminations:
        ifaces = sorted({d.interface_index for d in damage.delaminations})
        cmap = plt.get_cmap(ELLIPSE_CMAP)
        norm = plt.Normalize(vmin=min(ifaces), vmax=max(ifaces) + 1)
        for d in damage.delaminations:
            color = cmap(norm(d.interface_index))
            ax.add_patch(
                mpatches.Ellipse(
                    xy=d.centroid_mm,
                    width=2 * d.major_mm,
                    height=2 * d.minor_mm,
                    angle=d.orientation_deg,
                    facecolor=color,
                    alpha=0.35,
                    edgecolor=color,
                    linewidth=1.0,
                )
            )
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="Interface index")
    else:
        ax.text(
            panel.Lx_mm / 2,
            panel.Ly_mm / 2,
            "no damage",
            ha="center",
            va="center",
            color="grey",
            fontsize=14,
        )

    fig.tight_layout()
    return fig


def plot_knockdown_curve(
    energies_J: Sequence[float],
    knockdowns: Sequence[float],
    tier_label: str = "",
    title: str | None = None,
):
    """Line plot of knockdown vs impact energy."""
    fig, ax = plt.subplots(figsize=FIGSIZE_LANDSCAPE)
    color = COLORS.get(tier_label, COLORS["knockdown"])
    ax.plot(energies_J, knockdowns, "-o", color=color, label=tier_label or "knockdown")
    ax.set_xlabel("Impact energy [J]")
    ax.set_ylabel("Strength retention (knockdown)")
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    if tier_label:
        ax.legend()
    ax.set_title(title or "BVID knockdown curve")
    fig.tight_layout()
    return fig


def plot_tier_comparison(
    energies_J: Sequence[float],
    results_per_tier: Dict[str, Sequence[float]],
    title: str | None = None,
):
    """Overlaid knockdown curves for multiple tiers."""
    fig, ax = plt.subplots(figsize=FIGSIZE_LANDSCAPE)
    for tier, kd in results_per_tier.items():
        color = COLORS.get(tier, None)
        ax.plot(energies_J, kd, "-o", color=color, label=tier)
    ax.set_xlabel("Impact energy [J]")
    ax.set_ylabel("Strength retention (knockdown)")
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()
    ax.set_title(title or "BVID tier comparison")
    fig.tight_layout()
    return fig
