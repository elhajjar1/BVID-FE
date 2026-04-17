"""Knockdown curve tab: embedded matplotlib knockdown-vs-energy plot."""

from __future__ import annotations

from typing import Sequence

from PyQt6.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from bvidfe.viz.plots_2d import plot_knockdown_curve


class KnockdownTab(QWidget):
    """Embedded knockdown-curve tab. Populated via ``update_series``."""

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
