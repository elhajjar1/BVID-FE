import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import json

from bvidfe.analysis import AnalysisConfig
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.gui.config_io import config_to_dict, config_from_dict
from bvidfe.impact.mapping import ImpactEvent


def _sample_cfg():
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(30.0, ImpactorGeometry(20.0, "flat"), mass_kg=4.0),
    )


def test_config_round_trip():
    cfg = _sample_cfg()
    d = config_to_dict(cfg)
    cfg2 = config_from_dict(d)
    assert cfg2.material == cfg.material
    assert cfg2.layup_deg == cfg.layup_deg
    assert cfg2.panel.Lx_mm == cfg.panel.Lx_mm
    assert cfg2.impact is not None
    assert cfg2.impact.impactor.shape == "flat"
    assert cfg2.impact.impactor.diameter_mm == 20.0
    assert cfg2.impact.mass_kg == 4.0


def test_config_json_serializable():
    cfg = _sample_cfg()
    d = config_to_dict(cfg)
    s = json.dumps(d)  # must not raise
    assert isinstance(s, str)


def test_damage_only_config_round_trip():
    from bvidfe.damage.state import DamageState, DelaminationEllipse

    ds = DamageState(
        [DelaminationEllipse(3, (75, 50), 20, 12, 45)],
        dent_depth_mm=0.4,
        fiber_break_radius_mm=1.5,
    )
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(100, 50),
        loading="tension",
        tier="semi_analytical",
        damage=ds,
    )
    d = config_to_dict(cfg)
    cfg2 = config_from_dict(d)
    assert cfg2.damage is not None
    assert len(cfg2.damage.delaminations) == 1
    assert cfg2.damage.delaminations[0].interface_index == 3


def test_main_window_has_file_menu(qtbot):
    from bvidfe.gui.main_window import BvidMainWindow

    w = BvidMainWindow()
    qtbot.addWidget(w)
    actions = [a.text() for a in w.menuBar().actions()]
    # File menu should exist (title "&File" or "File")
    assert any("File" in t for t in actions)
