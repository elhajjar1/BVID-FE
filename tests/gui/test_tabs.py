"""Tests for GUI result tabs (Task 13.3)."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


@pytest.fixture
def sample_result():
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5),
    )
    return BvidAnalysis(cfg).run()


def test_summary_tab_updates(qtbot, sample_result):
    from bvidfe.gui.tabs.summary_tab import SummaryTab

    tab = SummaryTab()
    qtbot.addWidget(tab)
    tab.update(sample_result)
    # Summary should display the knockdown value
    assert f"{sample_result.knockdown:.3f}" in tab.text_area.toPlainText()


def test_damage_map_tab_updates(qtbot, sample_result):
    from bvidfe.gui.tabs.damage_map_tab import DamageMapTab

    tab = DamageMapTab()
    qtbot.addWidget(tab)
    tab.update(sample_result, panel=PanelGeometry(150, 100))
    # Canvas should exist and have a figure
    assert tab.canvas.figure is not None


def test_knockdown_tab_updates_with_series(qtbot, sample_result):
    from bvidfe.gui.tabs.knockdown_tab import KnockdownTab

    tab = KnockdownTab()
    qtbot.addWidget(tab)
    # Provide a pre-computed series
    tab.update_series([5, 10, 20, 30], [0.95, 0.80, 0.55, 0.35], tier_label="empirical")
    assert tab.canvas.figure is not None


def test_placeholder_tab_shows_message(qtbot):
    from bvidfe.gui.tabs.placeholder_tab import PlaceholderTab

    tab = PlaceholderTab("3D Mesh (v0.2.0)")
    qtbot.addWidget(tab)
    assert "v0.2.0" in tab.label.text()


def test_main_window_has_tabs(qtbot):
    from bvidfe.gui.main_window import BvidMainWindow

    w = BvidMainWindow()
    qtbot.addWidget(w)
    tab_titles = [w.results_tabs.tabText(i) for i in range(w.results_tabs.count())]
    assert "Summary" in tab_titles
    assert "Damage Map" in tab_titles
    assert "Knockdown Curve" in tab_titles
