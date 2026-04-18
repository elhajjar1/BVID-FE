"""Tests for the 3D mesh viewer tab (Task: wire pyvistaqt interactor)."""

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
    """The tab constructs without touching VTK. Plotter is lazily created."""
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    # Lazy design: no plotter until the user clicks the Open button.
    assert tab.plotter is None
    assert tab._render_button is not None


def test_mesh_3d_tab_caches_state_without_rendering(qtbot, sample_cfg_and_result):
    """update(config, result) before the user clicks Open should cache only."""
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    cfg, result = sample_cfg_and_result
    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    tab.update(cfg, result)
    # No plotter created yet — but the state is cached for later rendering.
    assert tab.plotter is None
    assert tab._pending_results is result


def test_mesh_3d_tab_render_after_click_creates_stub(qtbot, sample_cfg_and_result):
    """Clicking Open with a cached result creates the plotter (stub headless)."""
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    cfg, result = sample_cfg_and_result
    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    tab.update(cfg, result)
    tab._on_render_clicked()  # simulates button click
    assert tab.plotter is not None
    assert len(tab.plotter.actors) >= 1


def test_mesh_3d_tab_update_with_empty_damage(qtbot):
    """Empty damage shouldn't crash the render path."""
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
    tab.update(cfg, result)
    tab._on_render_clicked()
    assert tab.plotter is not None
