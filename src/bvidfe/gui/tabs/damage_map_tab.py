"""Damage map tab: embedded matplotlib damage-map plot."""

from __future__ import annotations

from PyQt6.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from bvidfe.analysis import AnalysisResults
from bvidfe.core.geometry import PanelGeometry
from bvidfe.viz.plots_2d import plot_damage_map


class DamageMapTab(QWidget):
    """Embedded damage-map plot tab."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.canvas = FigureCanvas(Figure(figsize=(6, 6)))
        lay = QVBoxLayout(self)
        lay.addWidget(self.canvas)

    def update(self, results: AnalysisResults, panel: PanelGeometry) -> None:  # type: ignore[override]
        fig = plot_damage_map(results.damage, panel)
        self.canvas.figure = fig
        self.canvas.draw()
