"""Knockdown curve tab: embedded matplotlib knockdown-vs-energy plot."""

from __future__ import annotations

from typing import Dict, Sequence

from PyQt6.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from bvidfe.viz.plots_2d import plot_knockdown_curve, plot_tier_comparison


class KnockdownTab(QWidget):
    """Embedded knockdown-curve tab. Populated via ``update_series`` or
    ``update_tier_comparison``."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.canvas = FigureCanvas(Figure(figsize=(8, 5)))
        lay = QVBoxLayout(self)
        lay.addWidget(self.canvas)

    def update_series(
        self,
        energies_J: Sequence[float],
        knockdowns: Sequence[float],
        tier_label: str = "",
    ) -> None:
        fig = plot_knockdown_curve(list(energies_J), list(knockdowns), tier_label=tier_label)
        self.canvas.figure = fig
        self.canvas.draw()

    def update_tier_comparison(
        self,
        energies_J: Sequence[float],
        knockdowns_by_tier: Dict[str, Sequence[float]],
    ) -> None:
        """Plot multiple tiers' knockdown curves overlaid on the same axes."""
        fig = plot_tier_comparison(
            list(energies_J), {k: list(v) for k, v in knockdowns_by_tier.items()}
        )
        self.canvas.figure = fig
        self.canvas.draw()
