"""``AnalysisConfig.transverse_stress_ratio`` (issue #115): validation,
JSON round trip and backward compatibility."""

from __future__ import annotations

import dataclasses
import json
import math
from pathlib import Path

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, PanelGeometry
from bvidfe.analysis.bvid import _config_snapshot_dict

SCHEMA = Path(__file__).resolve().parents[2] / "docs" / "schemas" / "analysis_config.json"


def _kw(**kw):
    base = dict(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        impact=ImpactEvent(energy_J=10.0),
    )
    base.update(kw)
    return base


def test_default_is_uniaxial():
    assert AnalysisConfig(**_kw()).transverse_stress_ratio == 0.0


def test_field_is_last_so_positional_construction_is_unchanged():
    names = [f.name for f in dataclasses.fields(AnalysisConfig)]
    assert names[-1] == "transverse_stress_ratio"
    assert names[:-1] == [
        "material",
        "layup_deg",
        "ply_thickness_mm",
        "panel",
        "loading",
        "tier",
        "impact",
        "damage",
        "mesh",
        "impact_panel",
    ]
    panel = PanelGeometry(150.0, 100.0)
    cfg = AnalysisConfig(
        "IM7/8552", [0, 90, 90, 0], 0.152, panel, "tension", "empirical", ImpactEvent(5.0)
    )
    assert (cfg.loading, cfg.tier, cfg.transverse_stress_ratio) == ("tension", "empirical", 0.0)


@pytest.mark.parametrize("beta", [-1.0, -0.5, 0.0, 0.25, 1.0, 1, -1, 0])
def test_ratios_in_range_are_accepted_as_floats(beta):
    cfg = AnalysisConfig(**_kw(transverse_stress_ratio=beta))
    assert cfg.transverse_stress_ratio == beta
    assert type(cfg.transverse_stress_ratio) is float


def test_negative_zero_is_stored_as_zero():
    cfg = AnalysisConfig(**_kw(transverse_stress_ratio=-0.0))
    assert math.copysign(1.0, cfg.transverse_stress_ratio) == 1.0
    assert cfg == AnalysisConfig(**_kw())


@pytest.mark.parametrize("beta", [1.0000001, -1.5, 3.0, -100.0])
def test_out_of_range_names_the_rotation_workaround(beta):
    with pytest.raises(ValueError) as exc:
        AnalysisConfig(**_kw(transverse_stress_ratio=beta))
    msg = str(exc.value)
    assert "transverse_stress_ratio" in msg and "between -1 and 1" in msg
    assert "rotate the layup by 90 deg" in msg
    assert "Lx_mm" in msg and "Ly_mm" in msg
    assert f"{1.0 / beta:g}" in msg  # the inverse ratio to use


@pytest.mark.parametrize("beta", [True, False, math.nan, math.inf, -math.inf, "0.5", None])
def test_non_numbers_bools_and_non_finite_are_rejected(beta):
    with pytest.raises(ValueError, match="transverse_stress_ratio"):
        AnalysisConfig(**_kw(transverse_stress_ratio=beta))


def test_from_dict_without_the_key_loads_as_uniaxial():
    """A config JSON written before the field existed."""
    snap = _config_snapshot_dict(AnalysisConfig(**_kw()))
    assert snap["transverse_stress_ratio"] == 0.0
    del snap["transverse_stress_ratio"]
    assert AnalysisConfig.from_dict(snap) == AnalysisConfig(**_kw())


@pytest.mark.parametrize("beta", [-0.75, 0.5])
def test_from_dict_round_trips_the_ratio(beta):
    cfg = AnalysisConfig(**_kw(transverse_stress_ratio=beta))
    snap = json.loads(json.dumps(_config_snapshot_dict(cfg)))
    assert snap["transverse_stress_ratio"] == beta
    assert AnalysisConfig.from_dict(snap) == cfg


def test_from_dict_rejects_an_out_of_range_ratio():
    snap = _config_snapshot_dict(AnalysisConfig(**_kw()))
    snap["transverse_stress_ratio"] = 2.0
    with pytest.raises(ValueError, match="between -1 and 1"):
        AnalysisConfig.from_dict(snap)


def test_run_snapshot_carries_the_ratio():
    r = BvidAnalysis(AnalysisConfig(**_kw(transverse_stress_ratio=-0.25))).run()
    assert r.config_snapshot["transverse_stress_ratio"] == -0.25
    assert AnalysisConfig.from_dict(r.to_dict()["config_snapshot"]).transverse_stress_ratio == -0.25


def test_schema_states_the_range_and_the_sign_convention():
    prop = json.loads(SCHEMA.read_text())["properties"]["transverse_stress_ratio"]
    assert prop["type"] == "number"
    assert (prop["minimum"], prop["maximum"], prop["default"]) == (-1.0, 1.0, 0.0)
    assert "sigma_y / sigma_x" in prop["description"]
    assert "Experimental" in prop["description"]
    assert "transverse_stress_ratio" not in json.loads(SCHEMA.read_text())["required"]
