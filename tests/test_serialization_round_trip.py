"""Results and config JSON round-trip (issue #170).

``AnalysisResults.to_dict()["damage"]`` is the schema-v1.0 C-scan layout
(with ``schema_version``), so a predicted damage state can be re-run
inspection-driven; ``AnalysisConfig.from_dict`` loads ``config_snapshot``
(and ``dataclasses.asdict``) output back into an equal config.
"""

from __future__ import annotations

import importlib
import json
import logging
from dataclasses import asdict, replace
from pathlib import Path

import pytest

import bvidfe.cli as cli
from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.analysis.bvid import _config_snapshot_dict
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.io import (
    SCHEMA_VERSION,
    CScanSchemaError,
    damage_state_from_dict,
    damage_state_to_dict,
    load_cscan_json,
)
from bvidfe.impact.mapping import ImpactEvent

REPO_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_CSCAN = REPO_ROOT / "examples" / "sample_cscan.json"
SCHEMA_DIR = REPO_ROOT / "docs" / "schemas"
LAYUP = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]

CLI_BASE = [
    "--material",
    "IM7/8552",
    "--layup",
    "0,45,-45,90,90,-45,45,0",
    "--thickness",
    "0.152",
    "--panel",
    "150x100",
]


@pytest.fixture(autouse=True)
def _restore_fe3d_logger():
    """The in-process CLI calls attach a stderr handler; undo it per test."""
    log = logging.getLogger("bvidfe.fe3d")
    handlers, level = list(log.handlers), log.level
    yield
    log.handlers[:] = handlers
    log.setLevel(level)


def _cli(args, capsys):
    try:
        code = cli.main(CLI_BASE + args)
    except SystemExit as exc:
        code = exc.code
    out, err = capsys.readouterr()
    return code, out, err


def _impact_cfg(**kw) -> AnalysisConfig:
    args = dict(
        material="IM7/8552",
        layup_deg=list(LAYUP),
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        impact=ImpactEvent(energy_J=15.0),
    )
    args.update(kw)
    return AnalysisConfig(**args)


def _custom_material():
    return replace(MATERIAL_LIBRARY["IM7/8552"], name="custom", E11=150000.0, Zt=60.0)


CONFIGS = {
    "impact_defaults": lambda: _impact_cfg(),
    "impact_everything": lambda: _impact_cfg(
        ply_thickness_mm=[0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1],
        panel=PanelGeometry(80.0, 50.0, boundary="clamped"),
        impact_panel=PanelGeometry(60.0, 60.0, boundary="clamped"),
        impact=ImpactEvent(
            energy_J=8.0,
            impactor=ImpactorGeometry(diameter_mm=12.7, shape="conical"),
            mass_kg=3.0,
            location_xy_mm=(28.0, 31.0),
        ),
        loading="tension",
        tier="fe3d",
        mesh=MeshParams(elements_per_ply=2, in_plane_size_mm=2.5, cohesive_zone_factor=1.5),
    ),
    "custom_material_damage": lambda: AnalysisConfig(
        material=_custom_material(),
        layup_deg=list(LAYUP),
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0, boundary="free"),
        damage=load_cscan_json(SAMPLE_CSCAN),
        tier="semi_analytical",
    ),
}


# -------------------------------------------------------- results damage


def test_results_damage_is_the_cscan_layout():
    result = BvidAnalysis(_impact_cfg()).run()
    d = result.to_dict()["damage"]
    assert d["schema_version"] == SCHEMA_VERSION
    assert d == damage_state_to_dict(result.damage)


def test_results_damage_loads_back():
    result = BvidAnalysis(_impact_cfg()).run()
    assert result.damage.delaminations
    assert damage_state_from_dict(result.to_dict()["damage"]) == result.damage
    via_json = json.loads(json.dumps(result.to_dict()))["damage"]
    assert damage_state_from_dict(via_json) == result.damage


def test_cli_results_damage_feeds_cscan(tmp_path, capsys):
    """``bvidfe ... --energy 10 > res.json; jq .damage res.json > dmg.json;
    bvidfe ... --cscan dmg.json`` used to fail on the missing schema_version."""
    code, out, err = _cli(["--energy", "10"], capsys)
    assert code == 0, err
    res = json.loads(out)
    dmg = tmp_path / "dmg.json"
    dmg.write_text(json.dumps(res["damage"]))
    code, out, err = _cli(["--cscan", str(dmg)], capsys)
    assert code == 0, err
    rerun = json.loads(out)
    assert rerun["knockdown"] == res["knockdown"]
    assert rerun["damage"] == res["damage"]


def test_cli_save_damage_writes_a_cscan_file(tmp_path, capsys):
    path = tmp_path / "saved.json"
    code, out, err = _cli(["--energy", "10", "--save-damage", str(path)], capsys)
    assert code == 0, err
    res = json.loads(out)
    assert json.loads(path.read_text()) == res["damage"]
    code, out, err = _cli(["--cscan", str(path), "--quick-json"], capsys)
    assert code == 0, err
    assert json.loads(out)["knockdown"] == res["knockdown"]


def test_cli_save_damage_with_cscan_and_several_tiers(tmp_path, capsys):
    path = tmp_path / "saved.json"
    code, out, err = _cli(
        [
            "--cscan",
            str(SAMPLE_CSCAN),
            "--tier",
            "empirical,semi_analytical",
            "--quick",
            "--save-damage",
            str(path),
        ],
        capsys,
    )
    assert code == 0, err
    assert load_cscan_json(path) == load_cscan_json(SAMPLE_CSCAN)
    assert len(out.splitlines()) == 2


def test_cli_save_damage_unwritable_is_a_clean_error(tmp_path, capsys):
    path = tmp_path / "no_such_dir" / "saved.json"
    code, out, err = _cli(["--energy", "10", "--quick", "--save-damage", str(path)], capsys)
    assert code == 1
    assert "bvidfe: error: --save-damage" in err
    assert "Traceback" not in err
    assert out == ""


# ----------------------------------------------------------- from_dict


@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_from_dict_round_trips_config_snapshot(name):
    cfg = CONFIGS[name]()
    snap = _config_snapshot_dict(cfg)
    assert AnalysisConfig.from_dict(snap) == cfg
    assert AnalysisConfig.from_dict(json.loads(json.dumps(snap))) == cfg


def test_from_dict_round_trips_a_run_config_snapshot():
    cfg = CONFIGS["custom_material_damage"]()
    result = BvidAnalysis(cfg).run()
    assert AnalysisConfig.from_dict(result.config_snapshot) == cfg
    as_json = json.loads(json.dumps(result.to_dict()))
    assert AnalysisConfig.from_dict(as_json["config_snapshot"]) == cfg


def test_from_dict_loads_the_cli_config_snapshot(capsys):
    code, out, err = _cli(
        ["--energy", "12", "--boundary", "clamped", "--location", "70", "45"], capsys
    )
    assert code == 0, err
    snap = json.loads(out)["config_snapshot"]
    cfg = AnalysisConfig.from_dict(snap)
    assert cfg.panel.boundary == "clamped"
    assert cfg.impact.location_xy_mm == (70.0, 45.0)
    assert BvidAnalysis(cfg).run().knockdown == json.loads(out)["knockdown"]


def test_from_dict_accepts_damage_with_schema_version():
    """The results' ``damage`` (C-scan layout) can be put into a config."""
    result = BvidAnalysis(_impact_cfg()).run()
    d = _config_snapshot_dict(_impact_cfg())
    d["impact"] = None
    d["damage"] = result.to_dict()["damage"]
    assert AnalysisConfig.from_dict(d).damage == result.damage


def test_from_dict_applies_dataclass_defaults():
    cfg = AnalysisConfig.from_dict(
        {
            "material": "IM7/8552",
            "layup_deg": LAYUP,
            "ply_thickness_mm": 0.152,
            "panel": {"Lx_mm": 150, "Ly_mm": 100},
            "impact": {"energy_J": 15.0},
        }
    )
    assert cfg == _impact_cfg()


def test_from_dict_keeps_instances():
    cfg = _impact_cfg()
    d = {f: getattr(cfg, f) for f in ("material", "layup_deg", "ply_thickness_mm", "panel")}
    d["impact"] = cfg.impact
    assert AnalysisConfig.from_dict(d) == cfg


@pytest.mark.parametrize(
    "mutate, match",
    [
        (lambda d: d.update(colour="red"), "colour"),
        (lambda d: d.pop("panel"), "panel"),
        (lambda d: d["panel"].update(Lz_mm=1.0), "Lz_mm"),
        (lambda d: d["impact"].update(speed=3.0), "speed"),
        (lambda d: d["impact"]["impactor"].update(tip="flat"), "tip"),
        (lambda d: d.update(mesh={"in_plane_size": 2.5}), "in_plane_size"),
        (lambda d: d.update(material={"name": "x"}), "E11"),
        (lambda d: d.update(panel=[150, 100]), "panel"),
        (lambda d: d.update(material=3), "material"),
        (lambda d: d.update(tier="fe3D"), "tier"),
    ],
)
def test_from_dict_rejects_bad_input(mutate, match):
    d = json.loads(json.dumps(_config_snapshot_dict(_impact_cfg())))
    mutate(d)
    with pytest.raises(ValueError, match=match):
        AnalysisConfig.from_dict(d)


def test_from_dict_rejects_a_bad_damage_record():
    d = _config_snapshot_dict(CONFIGS["custom_material_damage"]())
    d["damage"]["delaminations"][0]["major_mm"] = -1.0
    with pytest.raises(CScanSchemaError, match="non-positive"):
        AnalysisConfig.from_dict(d)


def test_component_from_dicts():
    assert ImpactEvent.from_dict({"energy_J": 5}) == ImpactEvent(energy_J=5)
    assert ImpactEvent.from_dict(asdict(ImpactEvent(3.0, mass_kg=2.0))) == ImpactEvent(
        3.0, mass_kg=2.0
    )
    assert MeshParams.from_dict({}) == MeshParams()
    assert MeshParams.from_dict({"elements_per_ply": 2}) == MeshParams(elements_per_ply=2)
    assert PanelGeometry.from_dict({"Lx_mm": 1, "Ly_mm": 2}) == PanelGeometry(1, 2)
    assert ImpactorGeometry.from_dict({"shape": "flat"}) == ImpactorGeometry(shape="flat")


# --------------------------------------------------------- JSON schemas


def _validator(filename):
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((SCHEMA_DIR / filename).read_text())
    return jsonschema.Draft202012Validator(schema)


@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_config_snapshot_matches_the_config_schema(name):
    snap = json.loads(json.dumps(_config_snapshot_dict(CONFIGS[name]())))
    _validator("analysis_config.json").validate(snap)


def test_config_with_schema_versioned_damage_matches_the_config_schema():
    d = json.loads(json.dumps(_config_snapshot_dict(CONFIGS["custom_material_damage"]())))
    d["damage"]["schema_version"] = SCHEMA_VERSION
    _validator("analysis_config.json").validate(d)
    AnalysisConfig.from_dict(d)


def test_custom_material_without_the_zone_factor_still_loads():
    """A custom material written before ``fe3d_zone_inplane_factor``
    existed takes the default, and still matches the schema."""
    d = json.loads(json.dumps(_config_snapshot_dict(CONFIGS["custom_material_damage"]())))
    del d["material"]["fe3d_zone_inplane_factor"]
    _validator("analysis_config.json").validate(d)
    cfg = AnalysisConfig.from_dict(d)
    assert cfg.material.fe3d_zone_inplane_factor == 0.15
    assert cfg == CONFIGS["custom_material_damage"]()


def test_results_json_matches_the_results_schema():
    result = BvidAnalysis(_impact_cfg()).run()
    _validator("analysis_results.json").validate(json.loads(json.dumps(result.to_dict())))


# ------------------------------------------------------------------ app


def test_app_config_from_dict_uses_the_library_loader(monkeypatch):
    """The app's loader used to drop ``impact_panel`` and accept only a
    preset-name material."""
    pytest.importorskip("streamlit", reason="streamlit not installed")
    monkeypatch.setenv("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
    monkeypatch.syspath_prepend(str(REPO_ROOT))
    app_module = importlib.import_module("app")
    for name in sorted(CONFIGS):
        cfg = CONFIGS[name]()
        snap = json.loads(json.dumps(_config_snapshot_dict(cfg)))
        assert app_module._config_from_dict(snap) == cfg
