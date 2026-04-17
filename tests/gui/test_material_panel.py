import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from bvidfe.gui.panels.material_panel import MaterialPanel


def test_material_panel_populates_library(qtbot):
    panel = MaterialPanel()
    qtbot.addWidget(panel)
    # All four presets should be present in the dropdown
    names = [panel.material_combo.itemText(i) for i in range(panel.material_combo.count())]
    assert "IM7/8552" in names
    assert "AS4/3501-6" in names
    assert "T700/2510" in names
    assert "T800/epoxy" in names


def test_material_panel_config_changed_signal_fires_on_selection(qtbot):
    panel = MaterialPanel()
    qtbot.addWidget(panel)
    with qtbot.waitSignal(panel.configChanged, timeout=1000):
        # Switch to the second preset
        panel.material_combo.setCurrentIndex(1)


def test_material_panel_layup_editable(qtbot):
    panel = MaterialPanel()
    qtbot.addWidget(panel)
    panel.layup_edit.setText("0, 45, -45, 90")
    assert panel.get_layup_deg() == [0.0, 45.0, -45.0, 90.0]


def test_material_panel_invalid_layup_returns_empty_list(qtbot):
    panel = MaterialPanel()
    qtbot.addWidget(panel)
    panel.layup_edit.setText("not a layup")
    assert panel.get_layup_deg() == []


def test_material_panel_ply_thickness_default(qtbot):
    panel = MaterialPanel()
    qtbot.addWidget(panel)
    assert abs(panel.get_ply_thickness_mm() - 0.152) < 1e-6
