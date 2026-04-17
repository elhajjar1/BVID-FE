"""3D mesh viewer tab powered by pyvistaqt + PyVista.

Renders the FE mesh built for the current analysis config, coloured by
damage factor, so the user can see where delaminations drop element
stiffness.

Why lazy initialization
-----------------------
``pyvistaqt.QtInteractor`` creates a VTK render window at widget
construction time.  On macOS this triggers full OpenGL context setup,
which can block the Qt main thread for several seconds.  Any
tab-switch to a tab containing a live QtInteractor also forces a
re-render — on a 10k-element mesh this adds another few seconds of
main-thread work each time.

To keep the app snappy, the Mesh3DTab defers all VTK work until the
user explicitly clicks a "Render 3D mesh" button.  The placeholder
button is cheap, lives in the tab by default, and the expensive
QtInteractor is only instantiated on demand.

Headless note
-------------
Under ``QT_QPA_PLATFORM=offscreen`` (CI / test runs) ``QtInteractor``
segfaults inside ``vtkRenderWindow.SetWindowInfo()``.  The tab
substitutes a lightweight ``_StubPlotter`` that exposes the same
surface used by tests (``clear``, ``add_mesh``, ``add_axes``,
``reset_camera``, ``actors``).
"""

from __future__ import annotations

import logging
import os
import time
from typing import Any, Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QLabel, QPushButton, QVBoxLayout, QWidget

from bvidfe.analysis import AnalysisConfig, AnalysisResults, MeshParams
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.viz.plots_3d import mesh_to_pyvista

_HEADLESS = os.environ.get("QT_QPA_PLATFORM") == "offscreen"
_log = logging.getLogger("bvidfe.gui")


class _StubPlotter:
    """Minimal stand-in for ``QtInteractor`` used in headless test runs."""

    def __init__(self) -> None:
        self.actors: dict[str, Any] = {}
        self._actor_counter = 0

    def clear(self) -> None:
        self.actors.clear()
        self._actor_counter = 0

    def add_mesh(self, mesh: Any, **kwargs: Any) -> None:
        key = f"actor_{self._actor_counter}"
        self.actors[key] = mesh
        self._actor_counter += 1

    def add_axes(self, **kwargs: Any) -> None:
        key = f"actor_{self._actor_counter}"
        self.actors[key] = "axes"
        self._actor_counter += 1

    def reset_camera(self) -> None:
        pass

    def close(self) -> None:
        pass


class Mesh3DTab(QWidget):
    """Lazy 3D mesh viewer. VTK initialization is deferred until the
    user clicks the render button so that switching to this tab never
    freezes the GUI with an OpenGL context creation on the main thread.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

        # Latest analysis state — stored here so that clicking Render
        # after an analysis completes immediately renders the latest mesh.
        self._pending_config: Optional[AnalysisConfig] = None
        self._pending_results: Optional[AnalysisResults] = None
        self._initialized: bool = False
        self.plotter: Any = None

        # Placeholder UI
        self._placeholder_label = QLabel(
            "3D mesh viewer is not initialized.\n\n"
            "Click below to render. Rendering a 10k-element mesh takes a\n"
            "few seconds on macOS because VTK needs to set up an OpenGL\n"
            "context on first use.",
            self,
        )
        self._placeholder_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._layout.addWidget(self._placeholder_label)

        self._render_button = QPushButton("Render 3D mesh", self)
        self._render_button.clicked.connect(self._on_render_clicked)
        self._layout.addWidget(self._render_button)

        if _HEADLESS:
            # In headless mode we auto-initialize the stub so tests don't
            # have to click a button.
            self._init_plotter()

    # ------------------------------------------------------------------
    # Public API — called from BvidMainWindow on analysis completion
    # ------------------------------------------------------------------

    def update(self, config: AnalysisConfig, results: AnalysisResults) -> None:  # type: ignore[override]
        """Store the latest config + results so the tab can render when asked.

        If the viewer has already been initialized (user has clicked render
        at least once), re-render immediately with the new result.
        """
        self._pending_config = config
        self._pending_results = results
        if self._initialized:
            self._render_pending()

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _on_render_clicked(self) -> None:
        if self._pending_results is None:
            self._placeholder_label.setText(
                "No analysis result yet. Run an analysis first,\nthen return "
                "here and click the render button."
            )
            return
        if not self._initialized:
            self._init_plotter()
        self._render_pending()

    def _init_plotter(self) -> None:
        """Create the QtInteractor (or stub) and swap it in for the placeholder."""
        t0 = time.time()
        _log.info("Mesh3DTab: initializing VTK plotter (first render)")
        # Remove placeholder widgets
        self._placeholder_label.hide()
        self._render_button.hide()
        self._layout.removeWidget(self._placeholder_label)
        self._layout.removeWidget(self._render_button)
        self._placeholder_label.deleteLater()
        self._render_button.deleteLater()

        if _HEADLESS:
            self.plotter = _StubPlotter()
        else:
            from pyvistaqt import QtInteractor

            self.plotter = QtInteractor(self, off_screen=False)
            self._layout.addWidget(self.plotter)
        self._initialized = True
        _log.info("Mesh3DTab: VTK plotter ready (%.2fs)", time.time() - t0)

    def _render_pending(self) -> None:
        assert self._pending_config is not None
        assert self._pending_results is not None
        t0 = time.time()
        config = self._pending_config
        results = self._pending_results
        if config.mesh is None:
            config.mesh = MeshParams()
        _log.info(
            "Mesh3DTab: rendering mesh (tier=%s, n_delam=%d)",
            config.tier,
            len(results.damage.delaminations),
        )
        fe_mesh = build_fe_mesh(config, results.damage)
        _log.info(
            "Mesh3DTab: mesh built (%d elements, %.2fs)",
            fe_mesh.n_elements,
            time.time() - t0,
        )
        grid = mesh_to_pyvista(fe_mesh)
        self.plotter.clear()
        self.plotter.add_mesh(
            grid,
            scalars="damage_factor",
            cmap="RdYlGn",
            clim=(0.0, 1.0),
            show_edges=False,
            scalar_bar_args={"title": "Damage factor"},
        )
        self.plotter.add_axes()
        self.plotter.reset_camera()
        _log.info("Mesh3DTab: render complete (%.2fs total)", time.time() - t0)

    def closeEvent(self, event):  # type: ignore[override]
        """Release the VTK render window cleanly when the tab is destroyed."""
        if self.plotter is not None:
            try:
                self.plotter.close()
            except Exception:  # noqa: BLE001
                pass
        super().closeEvent(event)
