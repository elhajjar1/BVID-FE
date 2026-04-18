"""Smoke tests for the BucklingTab and StressFieldTab (formerly placeholders)."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


@pytest.fixture
def result_semi_analytical():
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(60, 40),
        loading="compression",
        tier="semi_analytical",
        impact=ImpactEvent(20.0, ImpactorGeometry(), mass_kg=5.5),
    )
    return cfg, BvidAnalysis(cfg).run()


@pytest.fixture
def result_empirical():
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(40, 30),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(10.0, ImpactorGeometry(), mass_kg=5.5),
    )
    return cfg, BvidAnalysis(cfg).run()


def test_buckling_tab_constructs(qtbot):
    from bvidfe.gui.tabs.buckling_tab import BucklingTab

    tab = BucklingTab()
    qtbot.addWidget(tab)
    assert tab.canvas is not None
    assert len(tab.canvas.figure.axes) == 1  # placeholder axis


def test_buckling_tab_update_with_eigenvalues(qtbot, result_semi_analytical):
    from bvidfe.gui.tabs.buckling_tab import BucklingTab

    _, result = result_semi_analytical
    tab = BucklingTab()
    qtbot.addWidget(tab)
    tab.update(result)
    assert len(tab.canvas.figure.axes) == 1
    # For semi_analytical, eigenvalues are populated
    assert result.buckling_eigenvalues is not None


def test_buckling_tab_update_empirical_shows_note(qtbot, result_empirical):
    from bvidfe.gui.tabs.buckling_tab import BucklingTab

    _, result = result_empirical
    tab = BucklingTab()
    qtbot.addWidget(tab)
    tab.update(result)
    # No eigenvalues for empirical tier; the tab should still render.
    assert result.buckling_eigenvalues is None
    assert len(tab.canvas.figure.axes) == 1


def test_stress_field_tab_constructs(qtbot):
    from bvidfe.gui.tabs.stress_field_tab import StressFieldTab

    tab = StressFieldTab()
    qtbot.addWidget(tab)
    assert tab.canvas is not None


def test_stress_field_tab_updates(qtbot, result_empirical):
    from bvidfe.gui.tabs.stress_field_tab import StressFieldTab

    cfg, result = result_empirical
    tab = StressFieldTab()
    qtbot.addWidget(tab)
    tab.update(cfg, result)
    # After updating, the figure should have the heatmap + colorbar
    assert len(tab.canvas.figure.axes) >= 1
