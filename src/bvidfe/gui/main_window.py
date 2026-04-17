"""BVID-FE main window (QMainWindow skeleton)."""

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

from bvidfe.gui.panels.material_panel import MaterialPanel


class BvidMainWindow(QMainWindow):
    """Main application window for BVID-FE."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("BVID-FE")
        self.resize(1200, 800)

        # Status bar
        self.setStatusBar(QStatusBar(self))
        self.statusBar().showMessage("Ready")

        # Input panels in a left dock (scrollable)
        self.material_panel = MaterialPanel(self)
        panel_container = QWidget()
        layout = QVBoxLayout(panel_container)
        layout.addWidget(self.material_panel)
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

        # Central tabbed results area (empty placeholders for now)
        self.results_tabs = QTabWidget(self)
        placeholder = QLabel("Run an analysis to see results.")
        placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.results_tabs.addTab(placeholder, "Summary")
        self.setCentralWidget(self.results_tabs)
