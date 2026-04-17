"""BVID-FE main window (QMainWindow)."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PyQt6.QtCore import Qt

if TYPE_CHECKING:
    from bvidfe.gui.workers import AnalysisWorker, SweepWorker
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

        # Keep a reference to workers to prevent garbage-collection during run
        self._analysis_worker: AnalysisWorker | None = None
        self._sweep_worker: SweepWorker | None = None
        self.analysis_panel.runRequested.connect(self._run_analysis)
        self.sweep_panel.sweepRequested.connect(self._run_sweep)

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

    def _build_config(self):
        """Assemble an AnalysisConfig from current panel state."""
        from bvidfe.analysis import AnalysisConfig, MeshParams
        from bvidfe.core.geometry import PanelGeometry

        panel = PanelGeometry(
            Lx_mm=self.panel_panel.get_Lx_mm(),
            Ly_mm=self.panel_panel.get_Ly_mm(),
            boundary=self.panel_panel.get_boundary(),
        )
        mode = self.input_mode_panel.current_mode()
        impact = self.impact_panel.get_impact_event() if mode == "impact" else None
        damage = self.damage_panel.get_damage_state() if mode == "damage" else None
        mesh_params = None
        if self.analysis_panel.get_tier() == "fe3d":
            mesh_params = MeshParams(
                elements_per_ply=self.analysis_panel.get_elements_per_ply(),
                in_plane_size_mm=self.analysis_panel.get_in_plane_size_mm(),
            )
        return AnalysisConfig(
            material=self.material_panel.get_material_name(),
            layup_deg=self.material_panel.get_layup_deg(),
            ply_thickness_mm=self.material_panel.get_ply_thickness_mm(),
            panel=panel,
            loading=self.analysis_panel.get_loading(),
            tier=self.analysis_panel.get_tier(),
            impact=impact,
            damage=damage,
            mesh=mesh_params,
        )

    def _run_analysis(self) -> None:
        from bvidfe.gui.workers import AnalysisWorker

        try:
            cfg = self._build_config()
        except (ValueError, AssertionError) as exc:
            self.statusBar().showMessage(f"Invalid config: {exc}", 10000)
            return

        self.statusBar().showMessage("Running analysis\u2026")
        worker = AnalysisWorker(cfg)
        worker.resultReady.connect(self._on_analysis_ready)
        worker.error.connect(self._on_worker_error)
        worker.progress.connect(self._on_progress)
        worker.finished.connect(lambda: worker.deleteLater())
        self._analysis_worker = worker
        worker.start()

    def _run_sweep(self) -> None:
        from bvidfe.gui.workers import SweepWorker

        try:
            cfg = self._build_config()
        except (ValueError, AssertionError) as exc:
            self.statusBar().showMessage(f"Invalid config: {exc}", 10000)
            return
        energies = self.sweep_panel.get_energies_J()
        if not energies:
            self.statusBar().showMessage("No energies specified for sweep", 5000)
            return
        csv_path = self.sweep_panel.get_csv_path() or None

        self.statusBar().showMessage("Running sweep\u2026")
        worker = SweepWorker(cfg, energies_J=energies, csv_path=csv_path)
        worker.resultReady.connect(self._on_sweep_ready)
        worker.error.connect(self._on_worker_error)
        worker.progress.connect(self._on_progress)
        worker.finished.connect(lambda: worker.deleteLater())
        self._sweep_worker = worker
        worker.start()

    def _on_analysis_ready(self, result) -> None:
        self.summary_tab.update(result)
        self.damage_map_tab.update(result, panel=self.panel_panel_as_geometry())
        self.statusBar().showMessage(
            f"Analysis complete: knockdown = {result.knockdown:.3f}", 10000
        )

    def _on_sweep_ready(self, df) -> None:
        energies = df["energy_J"].tolist() if "energy_J" in df.columns else list(range(len(df)))
        knockdowns = df["knockdown"].tolist()
        self.knockdown_tab.update_series(
            energies,
            knockdowns,
            tier_label=self.analysis_panel.get_tier(),
        )
        self.statusBar().showMessage(f"Sweep complete: {len(df)} points", 10000)

    def _on_worker_error(self, tb: str) -> None:
        from PyQt6.QtWidgets import QMessageBox

        QMessageBox.critical(self, "Analysis error", tb)
        self.statusBar().showMessage("Analysis failed", 5000)

    def _on_progress(self, percent: int) -> None:
        self.statusBar().showMessage(f"Running\u2026 {percent}%")

    def panel_panel_as_geometry(self):
        from bvidfe.core.geometry import PanelGeometry

        return PanelGeometry(
            Lx_mm=self.panel_panel.get_Lx_mm(),
            Ly_mm=self.panel_panel.get_Ly_mm(),
            boundary=self.panel_panel.get_boundary(),
        )
