"""Buckling Mode tab — plots buckling eigenvalues when available.

For tiers that produce buckling eigenvalues (`semi_analytical` and `fe3d`)
this tab renders a bar chart showing the eigenvalues and the resulting
critical buckling stresses. For `empirical` it displays a static note
explaining that the empirical tier does not perform a buckling solve.
"""

from __future__ import annotations

from typing import Any, Optional

from PyQt6.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class BucklingTab(QWidget):
    """Buckling-mode visualization (matplotlib, main-thread safe)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = FigureCanvas(Figure(figsize=(8, 5)))
        layout.addWidget(self.canvas)
        self._draw_placeholder()

    def _draw_placeholder(self) -> None:
        fig = self.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.text(
            0.5,
            0.5,
            "Run a semi_analytical or fe3d analysis to see buckling eigenvalues.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="grey",
        )
        ax.set_axis_off()
        self.canvas.draw()

    def update(self, results: Any) -> None:  # type: ignore[override]
        eigs: Optional[list] = results.buckling_eigenvalues
        tier = results.tier_used
        fig = self.canvas.figure
        fig.clear()

        if not eigs:
            ax = fig.add_subplot(111)
            note = (
                f"Tier '{tier}' does not produce buckling eigenvalues.\n\n"
                "The empirical tier uses closed-form strength knockdown and skips\n"
                "the eigenvalue step entirely. Switch the tier to semi_analytical\n"
                "(Rayleigh-Ritz sublaminate buckling) or fe3d (geometric-stiffness\n"
                "FE buckling) to populate this tab."
            )
            ax.text(
                0.5,
                0.5,
                note,
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=11,
                family="monospace",
            )
            ax.set_axis_off()
            self.canvas.draw()
            return

        # Bar chart of the first (up to 6) eigenvalues
        eigs_to_plot = eigs[:6]
        ax = fig.add_subplot(111)
        idx = list(range(1, len(eigs_to_plot) + 1))
        colors = ["#1f77b4"] + ["#cccccc"] * (len(eigs_to_plot) - 1)
        ax.bar(idx, eigs_to_plot, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Mode number")
        ax.set_ylabel("Buckling eigenvalue / load factor")
        ax.set_title(
            f"Buckling eigenvalues ({tier} tier)\n" f"mode 1 critical = {eigs_to_plot[0]:.3g}"
        )
        ax.grid(axis="y", linestyle="--", alpha=0.3)
        for i, v in enumerate(eigs_to_plot):
            ax.text(idx[i], v, f"{v:.3g}", ha="center", va="bottom", fontsize=9)
        fig.tight_layout()
        self.canvas.draw()
