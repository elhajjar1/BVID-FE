"""2D matplotlib plots for BVID-FE: damage map, knockdown curve, tier comparison
and the (experimental) biaxial failure envelope.

All plot functions accept an optional ``fig`` argument. When ``fig`` is ``None``
(default) they create a fresh ``Figure`` via ``plt.subplots`` — convenient for
standalone scripts and tests that want to call ``fig.savefig(...)``. When
``fig`` is an existing ``Figure`` (for example one owned by an embedding UI
canvas) the function clears it in place and draws into it, so the caller keeps
its own figure and canvas instead of swapping in a pyplot-managed one.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, Optional, Sequence

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from bvidfe.core.geometry import PanelGeometry
from bvidfe.damage.state import DamageState
from bvidfe.viz.style import COLORS, ELLIPSE_CMAP, FIGSIZE_LANDSCAPE, FIGSIZE_SQUARE

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd


def _prepare_axes(fig: Optional[Figure], figsize):
    """Return ``(fig, ax)``: if ``fig`` is None make a fresh one, else clear and
    add a single axes. Used by all plot_* functions so the embedded-figure path
    never relies on ``plt.subplots``."""
    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)
        return fig, ax
    fig.clear()
    ax = fig.add_subplot(111)
    return fig, ax


def plot_damage_map(
    damage: DamageState,
    panel: PanelGeometry,
    title: str | None = None,
    fig: Optional[Figure] = None,
):
    """Top-down plan view of ellipse delamination footprints + panel outline.

    Ellipses are color-coded by interface index. Returns the matplotlib Figure.
    If ``fig`` is provided it is cleared and reused (safe for embedded UI
    canvases); otherwise a new Figure is created via ``plt.subplots``.
    """
    fig, ax = _prepare_axes(fig, FIGSIZE_SQUARE)
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

        # Impact-centroid marker: all ellipses in an impact-driven DamageState
        # share the same centroid (the impact point). For inspection-driven
        # DamageStates the centroids can differ; we plot each unique one.
        unique_centroids = {d.centroid_mm for d in damage.delaminations}
        for cx, cy in unique_centroids:
            ax.plot(
                cx,
                cy,
                marker="x",
                markersize=12,
                markeredgewidth=2,
                color="black",
                zorder=10,
                label="impact" if (cx, cy) == next(iter(unique_centroids)) else None,
            )
            # Fiber-break core on the back face: small filled red circle
            if damage.fiber_break_radius_mm > 0:
                ax.add_patch(
                    mpatches.Circle(
                        (cx, cy),
                        damage.fiber_break_radius_mm,
                        facecolor="red",
                        alpha=0.6,
                        edgecolor="darkred",
                        linewidth=1.0,
                        zorder=9,
                        label="fiber break" if (cx, cy) == next(iter(unique_centroids)) else None,
                    )
                )
        ax.legend(loc="upper right", framealpha=0.85)
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
    fig: Optional[Figure] = None,
):
    """Line plot of knockdown vs impact energy.

    If ``fig`` is provided it is cleared and reused; otherwise a new Figure
    is created.
    """
    fig, ax = _prepare_axes(fig, FIGSIZE_LANDSCAPE)
    color = COLORS.get(tier_label, COLORS["knockdown"])
    ax.plot(energies_J, knockdowns, "-o", color=color, label=tier_label or "knockdown")
    ax.set_xlabel("Impact energy [J]")
    ax.set_ylabel("Strength retention (knockdown) [-]")
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
    fig: Optional[Figure] = None,
):
    """Overlaid knockdown curves for multiple tiers.

    If ``fig`` is provided it is cleared and reused; otherwise a new Figure
    is created.
    """
    fig, ax = _prepare_axes(fig, FIGSIZE_LANDSCAPE)
    for tier, kd in results_per_tier.items():
        color = COLORS.get(tier, None)
        ax.plot(energies_J, kd, "-o", color=color, label=tier)
    ax.set_xlabel("Impact energy [J]")
    ax.set_ylabel("Strength retention (knockdown) [-]")
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()
    ax.set_title(title or "BVID tier comparison")
    fig.tight_layout()
    return fig


def _envelope_points(df: "pd.DataFrame", s: float, strength_col: str, sigma_y_col: str):
    """``(ratio, sigma_x, sigma_y)`` arrays of the successful rows of a
    ``sweep_transverse_ratios`` frame, ordered by ratio; ``sigma_x`` is the
    signed axial stress ``s * |sigma_x|``."""
    ok = df.dropna(subset=[strength_col, sigma_y_col]).sort_values("transverse_stress_ratio")
    return (
        ok["transverse_stress_ratio"].to_numpy(dtype=float),
        s * ok[strength_col].to_numpy(dtype=float),
        ok[sigma_y_col].to_numpy(dtype=float),
    )


def biaxial_envelope_plot(
    df_compression: Optional["pd.DataFrame"] = None,
    df_tension: Optional["pd.DataFrame"] = None,
    title: str | None = None,
    fig: Optional[Figure] = None,
):
    """Pristine and damaged (sigma_x, sigma_y) failure points of a biaxial
    ratio sweep. **Experimental, not validated** (issue #115): no biaxial
    residual-strength data exist, so the envelope is not an allowable.

    ``df_compression`` and ``df_tension`` are
    :func:`~bvidfe.sweep.parametric_sweep.sweep_transverse_ratios` frames
    run with ``loading="compression"`` (sigma_x < 0, the left half) and
    ``"tension"`` (sigma_x > 0, the right half); either may be ``None``,
    not both. Each point is the failure stress along one load ray
    ``sigma_y / sigma_x = transverse_stress_ratio``: the pristine envelope
    (``pristine_MPa``, ``pristine_sigma_y_MPa``) in grey and the damaged one
    (``residual_MPa``, ``residual_sigma_y_MPa``) per ``tier_used``, in the
    tier's colour. The pristine envelope is the same for every tier and is
    drawn once. Failed rows (NaN) are left out. If ``fig`` is provided it
    is cleared and reused; otherwise a new Figure is created.
    """
    frames = [(df, s) for df, s in ((df_compression, -1.0), (df_tension, 1.0)) if df is not None]
    if not frames:
        raise ValueError("biaxial_envelope_plot needs df_compression, df_tension or both")
    fig, ax = _prepare_axes(fig, FIGSIZE_SQUARE)
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.6)
    ax.axvline(0.0, color="black", linewidth=0.8, alpha=0.6)

    n_points = 0
    labelled: set[str] = set()
    for df, s in frames:
        _, sx, sy = _envelope_points(df, s, "pristine_MPa", "pristine_sigma_y_MPa")
        if len(sx):
            ax.plot(
                sx,
                sy,
                "--o",
                color=COLORS["pristine"],
                linewidth=2.0,
                markersize=8,
                markerfacecolor="white",
                label=None if "pristine" in labelled else "pristine",
            )
            labelled.add("pristine")
        tiers = [t for t in dict.fromkeys(df["tier_used"].dropna()) if t]
        for tier in tiers:
            _, sx, sy = _envelope_points(
                df[df["tier_used"] == tier], s, "residual_MPa", "residual_sigma_y_MPa"
            )
            n_points += len(sx)
            label = f"damaged ({tier})"
            ax.plot(
                sx,
                sy,
                "-o",
                color=COLORS.get(tier, COLORS["knockdown"]),
                linewidth=2.0,
                markersize=8,
                label=None if label in labelled else label,
            )
            labelled.add(label)

    if n_points == 0:
        ax.text(
            0.5,
            0.5,
            "no successful sweep points",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="grey",
            fontsize=14,
        )
    else:
        ax.legend(loc="best", framealpha=0.85)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("σx at failure [MPa]")
    ax.set_ylabel("σy at failure [MPa]")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_title(title or "Biaxial failure envelope\n(experimental, not validated)")
    fig.tight_layout()
    return fig
