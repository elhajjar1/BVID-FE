"""3D mesh viewer tab powered by pyvistaqt + PyVista.

Renders the FE mesh built for the current analysis config, colored by
damage factor, so the user can see where delaminations drop element
stiffness. Available for every tier — for empirical / semi-analytical
runs we still build the mesh here (with default MeshParams) to visualize
the damage.

Headless / offscreen note
-------------------------
``pyvistaqt.QtInteractor`` calls ``vtkRenderWindow.SetWindowInfo()`` during
``__init__``, which segfaults on macOS under the ``offscreen`` Qt platform
because VTK's OpenGL integration cannot create a context against the
virtual window ID.  To allow the rest of the test suite to keep running
headless, ``Mesh3DTab`` detects the offscreen platform and substitutes a
lightweight ``_StubPlotter`` that exposes the same surface used by tests
(``clear``, ``add_mesh``, ``add_axes``, ``reset_camera``, ``actors``).
"""

from __future__ import annotations

import os
from typing import Any

from PyQt6.QtWidgets import QLabel, QVBoxLayout, QWidget

from bvidfe.analysis import AnalysisConfig, AnalysisResults, MeshParams
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.viz.plots_3d import mesh_to_pyvista

_HEADLESS = os.environ.get("QT_QPA_PLATFORM") == "offscreen"


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
    """PyVistaQt-embedded 3D mesh viewer.

    Falls back to a ``_StubPlotter`` when running under the offscreen Qt
    platform (headless CI / test runs) so that the tab can be constructed
    without segfaulting inside VTK's window-handle initialisation.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        if _HEADLESS:
            # Headless environment — use stub so we don't segfault in VTK
            self.plotter: Any = _StubPlotter()
            label = QLabel("3D viewer (headless mode — no render output)", self)
            layout.addWidget(label)
        else:
            from pyvistaqt import QtInteractor

            self.plotter = QtInteractor(self, off_screen=False)
            layout.addWidget(self.plotter)

    def update(self, config: AnalysisConfig, results: AnalysisResults) -> None:  # type: ignore[override]
        """Rebuild the mesh from the current config + result and render it."""
        # Ensure a MeshParams exists for empirical/semi_analytical runs
        cfg_for_mesh = config
        if cfg_for_mesh.mesh is None:
            # Attach a default MeshParams just for visualization
            cfg_for_mesh.mesh = MeshParams()

        fe_mesh = build_fe_mesh(cfg_for_mesh, results.damage)
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

    def closeEvent(self, event):  # type: ignore[override]
        """Release the VTK render window cleanly when the tab is destroyed."""
        self.plotter.close()
        super().closeEvent(event)
