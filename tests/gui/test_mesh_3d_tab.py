"""Tests for the through-thickness damage view (was 3D Mesh tab).

v0.2.0-dev: VTK in the GUI was abandoned after three embedding approaches
all deadlocked the Qt main loop on macOS. The tab is now a matplotlib
FigureCanvas with three orthographic damage projections.
"""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


@pytest.fixture
def sample_cfg_and_result():
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(60, 40),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(20.0, ImpactorGeometry(), mass_kg=5.5),
    )
    result = BvidAnalysis(cfg).run()
    return cfg, result


def test_mesh_3d_tab_constructs(qtbot):
    """Tab constructs without any VTK/OpenGL touch."""
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    # Matplotlib canvas is always available; VTK-style stub also present
    assert tab.canvas is not None
    assert tab.plotter is not None


def test_mesh_3d_tab_update_renders_canvas(qtbot, sample_cfg_and_result):
    """update() replaces the placeholder with the three damage views."""
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    cfg, result = sample_cfg_and_result
    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    tab.update(cfg, result)
    # After updating the canvas should have the expected 4-axes layout
    # (top view, side view, front view, info panel).
    fig = tab.canvas.figure
    assert len(fig.axes) == 4


def test_mesh_3d_tab_update_with_empty_damage(qtbot):
    """Empty damage should still render without crashing."""
    from bvidfe.damage.state import DamageState
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20, 10),
        loading="compression",
        tier="empirical",
        damage=DamageState([], dent_depth_mm=0.0),
    )
    result = BvidAnalysis(cfg).run()

    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    tab.update(cfg, result)  # must not raise
    assert len(tab.canvas.figure.axes) == 4
