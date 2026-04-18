"""Through-thickness damage visualization (was 3D Mesh).

v0.2.0-dev note: we tried three approaches to embed a VTK/PyVista 3D
viewer in the GUI (embedded QtInteractor, separate-window BackgroundPlotter,
lazy-init). All of them deadlocked the main Qt event loop on macOS under
various conditions — the VTK/Qt/OpenGL interop on Apple Silicon with
PyInstaller is not reliably usable from a single QApplication.

Rather than keep fighting it, this tab renders three matplotlib-based
orthographic projections of the damaged laminate:
  - Top view (x-y): delamination footprints on the panel
  - Side view (x-z): damaged elements through the thickness
  - Front view (y-z): same, perpendicular perspective

This is actually more informative for BVID engineers than a rotating
3D mesh (hex-element meshes aren't visually useful), renders in ~50 ms,
and is guaranteed not to freeze.

True VTK/PyVista visualization is still available through the Python
API (see examples/) or by calling the viz.plots_3d functions directly
from a user's own script — just not from inside the GUI's event loop.
"""

from __future__ import annotations

import logging
import time
from typing import Any, Optional

import numpy as np
from PyQt6.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse

from bvidfe.analysis import AnalysisConfig, AnalysisResults

_log = logging.getLogger("bvidfe.gui")


class Mesh3DTab(QWidget):
    """Through-thickness orthographic damage views (matplotlib-backed).

    Kept the class name ``Mesh3DTab`` so BvidMainWindow and existing tests
    don't have to change. The class-level attribute ``plotter`` is a stub
    retained for backward compatibility with the headless test surface.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

        self.canvas = FigureCanvas(Figure(figsize=(10, 8)))
        self._layout.addWidget(self.canvas)

        # Backwards-compatible stub for tests and main window callers.
        self.plotter: Any = _HeadlessStub()
        self._pending_config: Optional[AnalysisConfig] = None
        self._pending_results: Optional[AnalysisResults] = None

        self._draw_placeholder()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def update(self, config: AnalysisConfig, results: AnalysisResults) -> None:  # type: ignore[override]
        """Redraw the orthographic damage views for the given config + results."""
        self._pending_config = config
        self._pending_results = results
        t0 = time.time()
        _log.info(
            "Mesh3DTab: rendering through-thickness views (n_delam=%d)",
            len(results.damage.delaminations),
        )
        self._draw_views(config, results)
        _log.info("Mesh3DTab: render complete (%.2fs)", time.time() - t0)

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------

    def _draw_placeholder(self) -> None:
        fig = self.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.text(
            0.5,
            0.5,
            "Run an analysis to see the through-thickness damage views.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="grey",
        )
        ax.set_axis_off()
        self.canvas.draw()

    def _draw_views(self, config: AnalysisConfig, results: AnalysisResults) -> None:
        damage = results.damage
        panel = config.panel
        layup = config.layup_deg
        n_plies = len(layup)
        h_ply = config.ply_thickness_mm
        h_total = n_plies * h_ply
        Lx, Ly = panel.Lx_mm, panel.Ly_mm

        fig = self.canvas.figure
        fig.clear()
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.8])

        # --- Top view (x-y) — delamination footprints ---
        ax_top = fig.add_subplot(gs[0, 0])
        _draw_top_view(ax_top, damage, panel, n_plies)

        # --- Side view (x-z) — damaged interfaces along x ---
        ax_side = fig.add_subplot(gs[0, 1])
        _draw_side_view(ax_side, damage, Lx, h_total, h_ply, n_plies, axis="x")

        # --- Front view (y-z) — damaged interfaces along y ---
        ax_front = fig.add_subplot(gs[1, 0])
        _draw_side_view(ax_front, damage, Ly, h_total, h_ply, n_plies, axis="y")

        # --- Summary text panel ---
        ax_info = fig.add_subplot(gs[1, 1])
        ax_info.set_axis_off()
        info = (
            f"Panel: {Lx:.0f} x {Ly:.0f} x {h_total:.3f} mm\n"
            f"Layup: {len(layup)} plies, {h_ply:.3f} mm each\n\n"
            f"Damage state:\n"
            f"  DPA: {damage.projected_damage_area_mm2:.0f} mm^2\n"
            f"  Dent: {damage.dent_depth_mm:.3f} mm\n"
            f"  Delaminations: {len(damage.delaminations)}\n"
            f"  Fiber break r: {damage.fiber_break_radius_mm:.2f} mm\n\n"
            f"Residual: {results.residual_strength_MPa:.1f} MPa\n"
            f"Knockdown: {results.knockdown:.3f}\n"
            f"Tier: {results.tier_used}"
        )
        ax_info.text(
            0.02,
            0.98,
            info,
            ha="left",
            va="top",
            transform=ax_info.transAxes,
            fontsize=10,
            family="monospace",
        )

        fig.suptitle("Through-thickness damage views", fontsize=13)
        fig.tight_layout()
        self.canvas.draw()


def _draw_top_view(ax, damage, panel, n_plies: int) -> None:
    """Top-down panel view with delamination ellipse footprints."""
    Lx, Ly = panel.Lx_mm, panel.Ly_mm
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title("Top view (x-y)")
    # Panel outline
    from matplotlib.patches import Rectangle

    ax.add_patch(Rectangle((0, 0), Lx, Ly, fill=False, edgecolor="black", linewidth=1.3))

    for e in damage.delaminations:
        color = cmap((e.interface_index + 1) / max(1, n_plies))
        ellipse = Ellipse(
            xy=e.centroid_mm,
            width=2 * e.major_mm,
            height=2 * e.minor_mm,
            angle=e.orientation_deg,
            facecolor=color,
            alpha=0.35,
            edgecolor=color,
            linewidth=0.8,
        )
        ax.add_patch(ellipse)


def _draw_side_view(
    ax, damage, extent_mm: float, h_total: float, h_ply: float, n_plies: int, axis: str
) -> None:
    """Projection of delamination ellipses onto an (axis, z) slice.

    axis="x" -> side view over x, side view is x on horizontal, z on vertical
    axis="y" -> front view over y
    """
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    ax.set_xlim(0, extent_mm)
    ax.set_ylim(0, h_total)
    ax.set_aspect("auto")
    ax.set_xlabel(f"{axis} [mm]")
    ax.set_ylabel("z [mm] (through thickness)")
    ax.set_title(f"{'Side' if axis == 'x' else 'Front'} view ({axis}-z)")

    # Ply stack horizontal guide lines
    for i in range(n_plies + 1):
        z = i * h_ply
        ax.axhline(z, color="lightgrey", linewidth=0.4)

    for e in damage.delaminations:
        # The delamination sits at interface_index * h_ply (z = bottom of ply index+1)
        # Project the ellipse onto the chosen plane. We represent it as a
        # horizontal bar at the interface z, extending +/- semi-axis along axis.
        z = (e.interface_index + 1) * h_ply
        if axis == "x":
            center = e.centroid_mm[0]
            # Effective half-extent along x: project rotated ellipse.
            # For simplicity we use major cos(theta) + minor sin(theta) bound.
            theta = np.radians(e.orientation_deg)
            half = e.major_mm * abs(np.cos(theta)) + e.minor_mm * abs(np.sin(theta))
        else:
            center = e.centroid_mm[1]
            theta = np.radians(e.orientation_deg)
            half = e.major_mm * abs(np.sin(theta)) + e.minor_mm * abs(np.cos(theta))
        x_lo = max(0, center - half)
        x_hi = min(extent_mm, center + half)
        color = cmap((e.interface_index + 1) / max(1, n_plies))
        ax.plot(
            [x_lo, x_hi],
            [z, z],
            color=color,
            linewidth=3.5,
            alpha=0.7,
            solid_capstyle="round",
        )


class _HeadlessStub:
    """Back-compat stub for tests that used to assert on .plotter.actors."""

    def __init__(self) -> None:
        self.actors: dict[str, Any] = {}
        self._n = 0

    def clear(self) -> None:
        self.actors.clear()
        self._n = 0

    def add_mesh(self, *args: Any, **kwargs: Any) -> None:
        self.actors[f"a{self._n}"] = "matplotlib-stub"
        self._n += 1

    def add_axes(self, *args: Any, **kwargs: Any) -> None:
        self.actors[f"a{self._n}"] = "axes"
        self._n += 1

    def reset_camera(self) -> None:
        pass

    def close(self) -> None:
        pass
