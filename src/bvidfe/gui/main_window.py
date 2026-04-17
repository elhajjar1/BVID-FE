"""BVID-FE main window (QMainWindow)."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDockWidget,
    QLabel,
    QMainWindow,
    QScrollArea,
    QStatusBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from bvidfe.gui.panels.analysis_panel import AnalysisPanel
from bvidfe.gui.panels.damage_panel import DamagePanel
from bvidfe.gui.panels.impact_panel import ImpactPanel
from bvidfe.gui.panels.input_mode_panel import InputModePanel
from bvidfe.gui.panels.material_panel import MaterialPanel
from bvidfe.gui.panels.panel_panel import PanelPanel
from bvidfe.gui.panels.sweep_panel import SweepPanel


class BvidMainWindow(QMainWindow):
    """Main application window for BVID-FE."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("BVID-FE")
        self.resize(1200, 800)

        # Status bar
        self.setStatusBar(QStatusBar(self))
        self.statusBar().showMessage("Ready")

        # --- Input panels ---
        self.material_panel = MaterialPanel(self)
        self.panel_panel = PanelPanel(self)
        self.input_mode_panel = InputModePanel(self)
        self.impact_panel = ImpactPanel(self)
        self.damage_panel = DamagePanel(self)
        self.analysis_panel = AnalysisPanel(self)
        self.sweep_panel = SweepPanel(self)

        # Wire input mode toggle to enable/disable the relevant data panels
        self.input_mode_panel.configChanged.connect(self._on_mode_changed)
        self._on_mode_changed()  # set initial enabled state

        panel_container = QWidget()
        layout = QVBoxLayout(panel_container)
        for p in (
            self.material_panel,
            self.panel_panel,
            self.input_mode_panel,
            self.impact_panel,
            self.damage_panel,
            self.analysis_panel,
            self.sweep_panel,
        ):
            layout.addWidget(p)
        layout.addStretch(1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(panel_container)

        dock = QDockWidget("Inputs", self)
        dock.setAllowedAreas(
            Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea
        )
        dock.setWidget(scroll)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)

        # Central tabbed results area (placeholder for future result views)
        self.results_tabs = QTabWidget(self)
        placeholder = QLabel("Run an analysis to see results.")
        placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.results_tabs.addTab(placeholder, "Summary")
        self.setCentralWidget(self.results_tabs)

    def _on_mode_changed(self) -> None:
        """Toggle impact / damage panels based on the selected input mode."""
        mode = self.input_mode_panel.current_mode()
        self.impact_panel.setEnabled(mode == "impact")
        self.damage_panel.setEnabled(mode == "damage")
