"""Smoke tests for all six new GUI input panels and the updated BvidMainWindow."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_panel_panel_constructs(qtbot):
    from bvidfe.gui.panels.panel_panel import PanelPanel

    w = PanelPanel()
    qtbot.addWidget(w)
    assert w.get_Lx_mm() > 0 and w.get_Ly_mm() > 0
    assert w.get_boundary() in ("clamped", "simply_supported", "free")


def test_panel_panel_emits_config_changed(qtbot):
    from bvidfe.gui.panels.panel_panel import PanelPanel

    w = PanelPanel()
    qtbot.addWidget(w)
    with qtbot.waitSignal(w.configChanged, timeout=1000):
        w.lx_spin.setValue(200.0)


def test_input_mode_panel_default_is_impact(qtbot):
    from bvidfe.gui.panels.input_mode_panel import InputModePanel

    w = InputModePanel()
    qtbot.addWidget(w)
    assert w.current_mode() == "impact"


def test_input_mode_panel_toggle_to_damage(qtbot):
    from bvidfe.gui.panels.input_mode_panel import InputModePanel

    w = InputModePanel()
    qtbot.addWidget(w)
    with qtbot.waitSignal(w.configChanged, timeout=1000):
        w.damage_radio.setChecked(True)
    assert w.current_mode() == "damage"


def test_impact_panel_defaults(qtbot):
    from bvidfe.gui.panels.impact_panel import ImpactPanel

    w = ImpactPanel()
    qtbot.addWidget(w)
    ev = w.get_impact_event()
    assert ev.energy_J > 0
    assert ev.impactor.diameter_mm == 16.0
    assert ev.impactor.shape == "hemispherical"


def test_impact_panel_shape_dropdown(qtbot):
    from bvidfe.gui.panels.impact_panel import ImpactPanel

    w = ImpactPanel()
    qtbot.addWidget(w)
    items = [w.shape_combo.itemText(i) for i in range(w.shape_combo.count())]
    assert "hemispherical" in items
    assert "flat" in items
    assert "conical" in items


def test_damage_panel_empty_default(qtbot):
    from bvidfe.gui.panels.damage_panel import DamagePanel

    w = DamagePanel()
    qtbot.addWidget(w)
    ds = w.get_damage_state()
    assert ds.dent_depth_mm == 0.0
    assert ds.delaminations == []


def test_damage_panel_add_row_creates_delamination(qtbot):
    from bvidfe.gui.panels.damage_panel import DamagePanel

    w = DamagePanel()
    qtbot.addWidget(w)
    w.add_delamination_row(interface_index=2, cx=50, cy=40, major=12, minor=8, angle=30)
    ds = w.get_damage_state()
    assert len(ds.delaminations) == 1
    d = ds.delaminations[0]
    assert d.interface_index == 2
    assert d.major_mm == 12.0


def test_analysis_panel_defaults(qtbot):
    from bvidfe.gui.panels.analysis_panel import AnalysisPanel

    w = AnalysisPanel()
    qtbot.addWidget(w)
    assert w.get_tier() in ("empirical", "semi_analytical", "fe3d")
    assert w.get_loading() in ("compression", "tension")


def test_analysis_panel_run_clicked_emits_signal(qtbot):
    from bvidfe.gui.panels.analysis_panel import AnalysisPanel

    w = AnalysisPanel()
    qtbot.addWidget(w)
    with qtbot.waitSignal(w.runRequested, timeout=1000):
        w.run_button.click()


def test_sweep_panel_constructs(qtbot):
    from bvidfe.gui.panels.sweep_panel import SweepPanel

    w = SweepPanel()
    qtbot.addWidget(w)
    assert w.get_energies_J() == [5.0, 10.0, 20.0, 30.0, 40.0]


def test_main_window_has_all_panels(qtbot):
    from bvidfe.gui.main_window import BvidMainWindow

    w = BvidMainWindow()
    qtbot.addWidget(w)
    for name in (
        "material_panel",
        "panel_panel",
        "input_mode_panel",
        "impact_panel",
        "damage_panel",
        "analysis_panel",
        "sweep_panel",
    ):
        assert hasattr(w, name), f"missing {name}"
