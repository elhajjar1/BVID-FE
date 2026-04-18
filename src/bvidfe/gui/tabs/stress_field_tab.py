"""Stress Field tab — shows a damage-severity heatmap of the panel.

For fe3d runs we have per-element damage factors; this tab projects them
down to 2D (x-y) by summing over the through-thickness direction, giving
a "damage severity" map analogous to what a C-scan operator sees. For
other tiers we build the mesh from the current config+damage and show
the same map (the mesh is cheap to build for visualization).
"""

from __future__ import annotations

from PyQt6.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from bvidfe.analysis import AnalysisConfig, AnalysisResults, MeshParams
from bvidfe.analysis.fe_mesh import build_fe_mesh


class StressFieldTab(QWidget):
    """Damage-severity heatmap tab (through-thickness-summed)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = FigureCanvas(Figure(figsize=(8, 6)))
        layout.addWidget(self.canvas)
        self._draw_placeholder()

    def _draw_placeholder(self) -> None:
        fig = self.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.text(
            0.5,
            0.5,
            "Run an analysis to see the damage-severity heatmap.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="grey",
        )
        ax.set_axis_off()
        self.canvas.draw()

    def update(self, config: AnalysisConfig, results: AnalysisResults) -> None:  # type: ignore[override]
        """Build the mesh (if needed) and plot a through-thickness-summed
        damage severity map."""
        if config.mesh is None:
            config.mesh = MeshParams()
        mesh = build_fe_mesh(config, results.damage)

        # Sum (1 - damage_factor) over the through-thickness direction for
        # each (i, j) column to get a "damage depth" metric: 0 = no
        # delamination at any interface, n_plies-1 = every interface
        # delaminated at this (x, y).
        import math

        in_plane_size = config.mesh.in_plane_size_mm
        nx = max(1, math.ceil(config.panel.Lx_mm / in_plane_size))
        ny = max(1, math.ceil(config.panel.Ly_mm / in_plane_size))
        n_plies = len(config.layup_deg)
        nz = n_plies * config.mesh.elements_per_ply
        damage = 1.0 - mesh.damage_factors
        damage = damage.reshape((nz, ny, nx))
        damage_columns = damage.sum(axis=0)  # (ny, nx) — summed through thickness

        fig = self.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        if damage_columns.max() == 0:
            im = ax.imshow(
                damage_columns,
                origin="lower",
                extent=(0, config.panel.Lx_mm, 0, config.panel.Ly_mm),
                cmap="hot_r",
                vmin=0.0,
                vmax=1.0,
            )
        else:
            im = ax.imshow(
                damage_columns,
                origin="lower",
                extent=(0, config.panel.Lx_mm, 0, config.panel.Ly_mm),
                cmap="hot_r",
                vmin=0.0,
                vmax=max(1.0, damage_columns.max()),
            )
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title(
            f"Damage severity map (through-thickness sum)\n"
            f"tier={results.tier_used}  residual={results.residual_strength_MPa:.1f} MPa "
            f"KD={results.knockdown:.3f}"
        )
        cbar = fig.colorbar(im, ax=ax, shrink=0.85)
        cbar.set_label("Damage depth (stacked interfaces)")
        ax.set_aspect("equal")
        fig.tight_layout()
        self.canvas.draw()
