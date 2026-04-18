"""BVID-FE main window (QMainWindow)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction

if TYPE_CHECKING:
    from bvidfe.gui.workers import AnalysisWorker, SweepWorker
from PyQt6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QMainWindow,
    QMessageBox,
    QScrollArea,
    QStatusBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from bvidfe.gui.config_io import config_from_dict, config_to_dict

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
        from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab
        from bvidfe.gui.tabs.placeholder_tab import PlaceholderTab
        from bvidfe.gui.tabs.summary_tab import SummaryTab

        self.results_tabs = QTabWidget(self)
        self.summary_tab = SummaryTab(self)
        self.damage_map_tab = DamageMapTab(self)
        self.knockdown_tab = KnockdownTab(self)
        self.mesh_tab = Mesh3DTab(self)
        self.buckling_tab = PlaceholderTab("Buckling mode shape — available in v0.2.0", self)
        self.stress_tab = PlaceholderTab("Stress field contour — available in v0.2.0", self)

        self.results_tabs.addTab(self.summary_tab, "Summary")
        self.results_tabs.addTab(self.damage_map_tab, "Damage Map")
        self.results_tabs.addTab(self.knockdown_tab, "Knockdown Curve")
        self.results_tabs.addTab(self.mesh_tab, "Damage View")
        self.results_tabs.addTab(self.buckling_tab, "Buckling Mode")
        self.results_tabs.addTab(self.stress_tab, "Stress Field")
        self.setCentralWidget(self.results_tabs)

        self._last_result = None
        self._last_config = None
        self._build_file_menu()

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

        # fe3d mesh-size sanity check. Oversized problems can both hang the
        # GUI and, more seriously, crash the process outright via SIGSEGV from
        # scipy's native sparse solvers when memory is exhausted. The solver
        # itself enforces a hard cap (FE3D_MAX_DOF) via FE3DSizeError; this
        # GUI dialog is a friendly early check.
        if cfg.tier == "fe3d":
            from bvidfe.analysis.fe_mesh import estimate_fe_mesh_size
            from bvidfe.analysis.fe_tier import FE3D_MAX_DOF

            stats = estimate_fe_mesh_size(cfg)
            if stats["n_dof"] > FE3D_MAX_DOF:
                # Past the hard cap: block run outright, explain clearly
                QMessageBox.critical(
                    self,
                    "Mesh too large for fe3d tier",
                    f"The requested fe3d mesh has {stats['n_elements']:,} elements "
                    f"({stats['n_dof']:,} DOFs), which exceeds the safe-size cap "
                    f"of {FE3D_MAX_DOF:,} DOFs.\n\n"
                    f"At this size the pure-Python FE assembler and scipy sparse "
                    f"solvers can exhaust memory and crash the process.\n\n"
                    f"Please do one of:\n"
                    f"  \u2022 increase 'In-plane size (mm)' in the Analysis panel\n"
                    f"  \u2022 decrease 'Elements per ply'\n"
                    f"  \u2022 switch tier to empirical or semi_analytical",
                )
                self.statusBar().showMessage("Run cancelled: mesh too large", 5000)
                return

            if stats["n_elements"] > 50_000 or stats["n_dof"] > 150_000:
                msg = (
                    f"The requested fe3d mesh has {stats['n_elements']:,} elements "
                    f"({stats['n_dof']:,} DOFs).\n\n"
                    f"Python FE in this version is single-threaded; expect multi-minute "
                    f"wall time at this size.\n\n"
                    f"Suggested tweaks:\n"
                    f"  \u2022 increase 'In-plane size (mm)' in the Analysis panel\n"
                    f"  \u2022 decrease 'Elements per ply'\n"
                    f"  \u2022 switch tier to empirical or semi_analytical\n\n"
                    f"Run anyway?"
                )
                result = QMessageBox.question(
                    self,
                    "Large mesh",
                    msg,
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if result != QMessageBox.StandardButton.Yes:
                    self.statusBar().showMessage("Run cancelled", 3000)
                    return

        self._last_config = cfg
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

        # Same fe3d mesh-size guard as single-run. A sweep is N x single-run
        # cost, so the soft warning threshold here is stricter.
        if cfg.tier == "fe3d":
            from bvidfe.analysis.fe_mesh import estimate_fe_mesh_size
            from bvidfe.analysis.fe_tier import FE3D_MAX_DOF

            stats = estimate_fe_mesh_size(cfg)
            n_runs = len(energies)
            if stats["n_dof"] > FE3D_MAX_DOF:
                QMessageBox.critical(
                    self,
                    "Mesh too large for fe3d tier",
                    f"Sweep mesh has {stats['n_elements']:,} elements / "
                    f"{stats['n_dof']:,} DOFs, exceeding the safe-size cap "
                    f"of {FE3D_MAX_DOF:,}.\n\n"
                    f"Increase 'In-plane size (mm)', decrease 'Elements per ply', "
                    f"or switch to tier='empirical' / 'semi_analytical'.",
                )
                self.statusBar().showMessage("Sweep cancelled: mesh too large", 5000)
                return
            # Soft threshold: N runs * ~single-run cost. Warn at 10k elements * N.
            if stats["n_elements"] * n_runs > 50_000 or n_runs > 20:
                msg = (
                    f"Sweep: {n_runs} fe3d runs, each with {stats['n_elements']:,} "
                    f"elements ({stats['n_dof']:,} DOFs).\n\n"
                    f"Total projected cost \u2248 {n_runs * stats['n_elements']:,} "
                    f"element-solves. Expect multi-minute wall time.\n\n"
                    f"Suggested: start with tier='empirical' for a quick sweep, "
                    f"then use fe3d for spot-checks at key energies.\n\n"
                    f"Run anyway?"
                )
                result = QMessageBox.question(
                    self,
                    "Large fe3d sweep",
                    msg,
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if result != QMessageBox.StandardButton.Yes:
                    self.statusBar().showMessage("Sweep cancelled", 3000)
                    return

        self.statusBar().showMessage("Running sweep\u2026")
        worker = SweepWorker(cfg, energies_J=energies, csv_path=csv_path)
        worker.resultReady.connect(self._on_sweep_ready)
        worker.error.connect(self._on_worker_error)
        worker.progress.connect(self._on_progress)
        worker.finished.connect(lambda: worker.deleteLater())
        self._sweep_worker = worker
        worker.start()

    def _on_analysis_ready(self, result) -> None:
        self._last_result = result
        self.summary_tab.update(result)
        self.damage_map_tab.update(result, panel=self.panel_panel_as_geometry())
        if self._last_config is not None:
            self.mesh_tab.update(self._last_config, result)
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

    # ------------------------------------------------------------------
    # File menu
    # ------------------------------------------------------------------

    def _build_file_menu(self) -> None:
        menu = self.menuBar().addMenu("&File")

        save_cfg = QAction("Save Config\u2026", self)
        save_cfg.triggered.connect(self._save_config)
        menu.addAction(save_cfg)

        load_cfg = QAction("Load Config\u2026", self)
        load_cfg.triggered.connect(self._load_config)
        menu.addAction(load_cfg)

        menu.addSeparator()

        export_json = QAction("Export Results JSON\u2026", self)
        export_json.triggered.connect(self._export_results_json)
        menu.addAction(export_json)

        export_png = QAction("Export Damage Map PNG\u2026", self)
        export_png.triggered.connect(self._export_damage_png)
        menu.addAction(export_png)

    def _save_config(self) -> None:
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save AnalysisConfig", "bvidfe_config.json", "JSON (*.json)"
        )
        if not path_str:
            return
        try:
            cfg = self._build_config()
        except (ValueError, AssertionError) as exc:
            QMessageBox.warning(self, "Invalid config", str(exc))
            return
        Path(path_str).write_text(json.dumps(config_to_dict(cfg), indent=2))
        self.statusBar().showMessage(f"Saved config to {path_str}", 5000)

    def _load_config(self) -> None:
        path_str, _ = QFileDialog.getOpenFileName(self, "Load AnalysisConfig", "", "JSON (*.json)")
        if not path_str:
            return
        try:
            d = json.loads(Path(path_str).read_text())
            cfg = config_from_dict(d)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Failed to load", str(exc))
            return
        self._apply_config_to_panels(cfg)
        self.statusBar().showMessage(f"Loaded config from {path_str}", 5000)

    def _apply_config_to_panels(self, cfg) -> None:
        """Push config values back into the input panels."""
        # Material panel
        idx = self.material_panel.material_combo.findText(cfg.material)
        if idx >= 0:
            self.material_panel.material_combo.setCurrentIndex(idx)
        self.material_panel.layup_edit.setText(", ".join(f"{a:g}" for a in cfg.layup_deg))
        self.material_panel.thickness_spin.setValue(cfg.ply_thickness_mm)
        # Panel panel
        self.panel_panel.lx_spin.setValue(cfg.panel.Lx_mm)
        self.panel_panel.ly_spin.setValue(cfg.panel.Ly_mm)
        bi = self.panel_panel.boundary_combo.findText(cfg.panel.boundary)
        if bi >= 0:
            self.panel_panel.boundary_combo.setCurrentIndex(bi)
        # Input mode + impact/damage
        if cfg.impact is not None:
            self.input_mode_panel.impact_radio.setChecked(True)
            self.impact_panel.energy_spin.setValue(cfg.impact.energy_J)
            self.impact_panel.diameter_spin.setValue(cfg.impact.impactor.diameter_mm)
            si = self.impact_panel.shape_combo.findText(cfg.impact.impactor.shape)
            if si >= 0:
                self.impact_panel.shape_combo.setCurrentIndex(si)
            self.impact_panel.mass_spin.setValue(cfg.impact.mass_kg)
            self.impact_panel.location_x.setValue(cfg.impact.location_xy_mm[0])
            self.impact_panel.location_y.setValue(cfg.impact.location_xy_mm[1])
        elif cfg.damage is not None:
            self.input_mode_panel.damage_radio.setChecked(True)
            self.damage_panel.table.setRowCount(0)
            for d in cfg.damage.delaminations:
                self.damage_panel.add_delamination_row(
                    d.interface_index,
                    d.centroid_mm[0],
                    d.centroid_mm[1],
                    d.major_mm,
                    d.minor_mm,
                    d.orientation_deg,
                )
            self.damage_panel.dent_spin.setValue(cfg.damage.dent_depth_mm)
            self.damage_panel.fb_spin.setValue(cfg.damage.fiber_break_radius_mm)
        # Analysis panel
        ti = self.analysis_panel.tier_combo.findText(cfg.tier)
        if ti >= 0:
            self.analysis_panel.tier_combo.setCurrentIndex(ti)
        li = self.analysis_panel.loading_combo.findText(cfg.loading)
        if li >= 0:
            self.analysis_panel.loading_combo.setCurrentIndex(li)

    def _export_results_json(self) -> None:
        if self._last_result is None:
            QMessageBox.information(self, "No results", "Run an analysis first.")
            return
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Export results JSON", "bvidfe_results.json", "JSON (*.json)"
        )
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(self._last_result.to_dict(), indent=2, default=str))
        self.statusBar().showMessage(f"Exported results to {path_str}", 5000)

    def _export_damage_png(self) -> None:
        if self._last_result is None:
            QMessageBox.information(self, "No results", "Run an analysis first.")
            return
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Export damage map", "damage_map.png", "PNG (*.png)"
        )
        if not path_str:
            return
        self.damage_map_tab.canvas.figure.savefig(path_str, dpi=150)
        self.statusBar().showMessage(f"Saved PNG to {path_str}", 5000)
