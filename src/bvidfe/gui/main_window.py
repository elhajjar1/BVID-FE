"""BVID-FE main window (QMainWindow)."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDockWidget,
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

        # Central tabbed results area
        from bvidfe.gui.tabs.damage_map_tab import DamageMapTab
        from bvidfe.gui.tabs.knockdown_tab import KnockdownTab
        from bvidfe.gui.tabs.placeholder_tab import PlaceholderTab
        from bvidfe.gui.tabs.summary_tab import SummaryTab

        self.results_tabs = QTabWidget(self)
        self.summary_tab = SummaryTab(self)
        self.damage_map_tab = DamageMapTab(self)
        self.knockdown_tab = KnockdownTab(self)
        self.mesh_tab = PlaceholderTab("3D Mesh viewer — available in v0.2.0", self)
        self.buckling_tab = PlaceholderTab("Buckling mode shape — available in v0.2.0", self)
        self.stress_tab = PlaceholderTab("Stress field contour — available in v0.2.0", self)

        self.results_tabs.addTab(self.summary_tab, "Summary")
        self.results_tabs.addTab(self.damage_map_tab, "Damage Map")
        self.results_tabs.addTab(self.knockdown_tab, "Knockdown Curve")
        self.results_tabs.addTab(self.mesh_tab, "3D Mesh")
        self.results_tabs.addTab(self.buckling_tab, "Buckling Mode")
        self.results_tabs.addTab(self.stress_tab, "Stress Field")
        self.setCentralWidget(self.results_tabs)

    def _on_mode_changed(self) -> None:
        """Toggle impact / damage panels based on the selected input mode."""
        mode = self.input_mode_panel.current_mode()
        self.impact_panel.setEnabled(mode == "impact")
        self.damage_panel.setEnabled(mode == "damage")
