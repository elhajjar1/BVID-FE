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
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    assert tab.plotter is not None  # QtInteractor attached


def test_mesh_3d_tab_updates_with_result(qtbot, sample_cfg_and_result):
    from bvidfe.gui.tabs.mesh_3d_tab import Mesh3DTab

    cfg, result = sample_cfg_and_result
    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    tab.update(cfg, result)  # must not raise
    # After updating, the plotter should have at least one actor
    assert len(tab.plotter.actors) >= 1


def test_mesh_3d_tab_update_with_empty_damage(qtbot):
    """Empty damage shouldn't crash the viewer."""
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
    # Build a minimal AnalysisResults-ish object via actually running the analysis
    result = BvidAnalysis(cfg).run()

    tab = Mesh3DTab()
    qtbot.addWidget(tab)
    tab.update(cfg, result)  # must not raise
