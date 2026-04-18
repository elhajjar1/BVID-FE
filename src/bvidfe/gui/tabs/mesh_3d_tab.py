"""3D mesh viewer tab.

v0.2.0-dev note: embedding a ``pyvistaqt.QtInteractor`` inside a Qt
tab freezes the main window on macOS in several configurations — the
OpenGL context creation blocks the Qt main thread, and once created
the interactor can fail to paint the viewport until external events
force it.  Rather than fight that, this tab now opens a **separate
window** (``pyvistaqt.BackgroundPlotter``) when the user clicks the
Render button.  The tab itself stays lightweight and never freezes
the main GUI.

Headless note
-------------
Under ``QT_QPA_PLATFORM=offscreen`` (CI / test runs) ``BackgroundPlotter``
also creates an OpenGL context and can segfault.  The tab substitutes a
lightweight ``_StubPlotter`` that exposes the same surface used by
tests (``clear``, ``add_mesh``, ``add_axes``, ``reset_camera``,
``actors``).
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
    """Minimal stand-in for ``BackgroundPlotter`` used in headless test runs."""

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

    def render(self) -> None:
        pass

    def show(self) -> None:
        pass

    def close(self) -> None:
        pass


class Mesh3DTab(QWidget):
    """3D mesh viewer launcher. Opens a separate PyVista window on demand."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

        self._pending_config: Optional[AnalysisConfig] = None
        self._pending_results: Optional[AnalysisResults] = None
        self.plotter: Any = None

        self._info_label = QLabel(
            "3D mesh viewer\n\n"
            "Runs an analysis first, then click the button below to open\n"
            "the 3D view in a separate window. The separate-window approach\n"
            "avoids a known pyvistaqt/Qt embedding issue on macOS that can\n"
            "freeze the main window.",
            self,
        )
        self._info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._layout.addWidget(self._info_label)

        self._render_button = QPushButton("Open 3D mesh in separate window", self)
        self._render_button.clicked.connect(self._on_render_clicked)
        self._layout.addWidget(self._render_button)

    # ------------------------------------------------------------------
    # Public API — called from BvidMainWindow on analysis completion
    # ------------------------------------------------------------------

    def update(self, config: AnalysisConfig, results: AnalysisResults) -> None:  # type: ignore[override]
        """Cache latest state. If an external window is already open and
        the user wants auto-refresh, we re-render into it."""
        self._pending_config = config
        self._pending_results = results
        if self.plotter is not None and not _HEADLESS:
            try:
                self._render_pending()
            except Exception as exc:  # noqa: BLE001
                # If the external window was closed by the user, discard it.
                _log.info("Mesh3DTab: external plotter unavailable (%s); discarded", exc)
                self.plotter = None

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _on_render_clicked(self) -> None:
        if self._pending_results is None:
            self._info_label.setText(
                "No analysis result yet. Run an analysis first,\nthen click this button."
            )
            return
        try:
            self._render_pending()
        except Exception:  # noqa: BLE001
            _log.exception("Mesh3DTab: failed to open 3D viewer")
            self._info_label.setText(
                "Failed to open 3D viewer.\n"
                "See terminal for the traceback. (This usually means a\n"
                "VTK/OpenGL driver issue that is specific to macOS.)"
            )

    def _get_or_create_plotter(self) -> Any:
        """Return the current plotter, creating a fresh one if needed."""
        if self.plotter is not None:
            return self.plotter
        if _HEADLESS:
            self.plotter = _StubPlotter()
            return self.plotter
        from pyvistaqt import BackgroundPlotter

        _log.info("Mesh3DTab: opening BackgroundPlotter window")
        self.plotter = BackgroundPlotter(
            title="BVID-FE — 3D Mesh",
            window_size=(900, 700),
            show=True,
        )
        return self.plotter

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
        plotter = self._get_or_create_plotter()
        plotter.clear()
        plotter.add_mesh(
            grid,
            scalars="damage_factor",
            cmap="RdYlGn",
            clim=(0.0, 1.0),
            show_edges=False,
        )
        plotter.add_axes()
        plotter.reset_camera()
        _log.info("Mesh3DTab: render complete (%.2fs total)", time.time() - t0)

    def closeEvent(self, event):  # type: ignore[override]
        """Release the VTK render window cleanly when the tab is destroyed."""
        if self.plotter is not None:
            try:
                self.plotter.close()
            except Exception:  # noqa: BLE001
                pass
        super().closeEvent(event)
