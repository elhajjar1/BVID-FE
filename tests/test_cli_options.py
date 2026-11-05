"""CLI options aligned with the Python API and the app (issue #169).

``--loading`` defaults to ``compression``; ``--boundary``,
``--impactor-shape``, ``--location``, ``--in-plane-size`` and
``--elements-per-ply`` reach ``PanelGeometry`` / ``ImpactorGeometry`` /
``ImpactEvent`` / ``MeshParams``; options that would have no effect on the
run are usage errors; and the ply-angle and energy ranges are the API's.
The runs here are in-process and closed-form; the fe3d wiring is checked on
the config handed to ``BvidAnalysis`` without solving.
"""

from __future__ import annotations

import json
import logging
import re
import typing
from dataclasses import replace
from pathlib import Path

import pytest

import bvidfe.cli as cli
from bvidfe.analysis import BvidAnalysis
from bvidfe.core.geometry import BoundaryKind, ImpactorShape

REPO_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_CSCAN = str(REPO_ROOT / "examples" / "sample_cscan.json")

BASE = [
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


def _main(args, capsys):
    try:
        code = cli.main(BASE + args)
    except SystemExit as exc:
        code = exc.code
    out, err = capsys.readouterr()
    return code, out, err


def _json(args, capsys) -> dict:
    code, out, err = _main(args, capsys)
    assert code == 0, err
    return json.loads(out)


@pytest.fixture
def captured_configs(monkeypatch):
    """Record the configs the CLI runs, solving each on the empirical tier so
    an fe3d request is checked for its wiring without an FE solve."""
    seen = []

    class _Recorder:
        def __init__(self, config):
            seen.append(config)
            self._config = config

        def run(self):
            return BvidAnalysis(replace(self._config, tier="empirical")).run()

    monkeypatch.setattr(cli, "BvidAnalysis", _Recorder)
    return seen


# ------------------------------------------------------------------ help


def test_help_does_not_mention_a_json_option():
    """There is no ``--json`` flag; the full JSON is the default output."""
    text = cli._build_parser().format_help()
    assert not re.search(r"(?<![\w-])--json\b", text)


# --------------------------------------------------------------- loading


def test_loading_defaults_to_compression(capsys):
    data = _json(["--energy", "10"], capsys)
    assert data["config_snapshot"]["loading"] == "compression"
    explicit = _json(["--energy", "10", "--loading", "compression"], capsys)
    assert data["knockdown"] == explicit["knockdown"]


def test_loading_default_is_the_api_default():
    from bvidfe.analysis import AnalysisConfig

    (action,) = [a for a in cli._build_parser()._actions if a.dest == "loading"]
    assert action.default == AnalysisConfig.__dataclass_fields__["loading"].default


# --------------------------------------------------- new options' wiring


def test_boundary_choices_are_the_library_literal():
    (action,) = [a for a in cli._build_parser()._actions if a.dest == "boundary"]
    assert sorted(action.choices) == sorted(typing.get_args(BoundaryKind))
    assert action.default == "simply_supported"


def test_impactor_shape_choices_are_the_library_literal():
    (action,) = [a for a in cli._build_parser()._actions if a.dest == "impactor_shape"]
    assert sorted(action.choices) == sorted(typing.get_args(ImpactorShape))


def test_boundary_reaches_the_panel(capsys):
    data = _json(["--energy", "10", "--boundary", "clamped"], capsys)
    assert data["config_snapshot"]["panel"]["boundary"] == "clamped"
    default = _json(["--energy", "10"], capsys)
    assert default["config_snapshot"]["panel"]["boundary"] == "simply_supported"
    # The boundary changes the onset energy and so the predicted damage.
    assert data["dpa_mm2"] != default["dpa_mm2"]


def test_boundary_applies_to_a_cscan_run(capsys):
    """The boundary scales the sublaminate buckling load, so it is not an
    impact-only option."""
    data = _json(["--cscan", SAMPLE_CSCAN, "--boundary", "clamped"], capsys)
    assert data["config_snapshot"]["panel"]["boundary"] == "clamped"


def test_impactor_shape_reaches_the_impactor(capsys):
    flat = _json(["--energy", "20", "--impactor-shape", "flat"], capsys)
    hemi = _json(["--energy", "20"], capsys)
    assert flat["config_snapshot"]["impact"]["impactor"]["shape"] == "flat"
    assert hemi["config_snapshot"]["impact"]["impactor"]["shape"] == "hemispherical"
    assert flat["dpa_mm2"] != hemi["dpa_mm2"]


def test_location_reaches_the_impact_event(capsys):
    data = _json(["--energy", "10", "--location", "60", "45"], capsys)
    assert data["config_snapshot"]["impact"]["location_xy_mm"] == [60.0, 45.0]
    xs = {d["centroid_mm"][0] for d in data["damage"]["delaminations"]}
    assert xs == {60.0}


def test_location_default_is_the_panel_centre(capsys):
    """No ``--location`` is ``location_xy_mm=None``, the panel centre (#171)."""
    data = _json(["--energy", "10"], capsys)
    assert data["config_snapshot"]["impact"]["location_xy_mm"] is None
    assert {tuple(d["centroid_mm"]) for d in data["damage"]["delaminations"]} == {(75.0, 50.0)}


def test_location_origin_is_a_deprecated_centre(capsys):
    """``--location 0 0`` was the centre sentinel; it still runs a centre
    impact, with a deprecation note on stderr (the library's
    ``DeprecationWarning`` is hidden by default outside ``__main__``)."""
    code, out, err = _main(["--energy", "10", "--location", "0", "0"], capsys)
    assert code == 0, err
    assert "bvidfe: warning:" in err and "deprecated" in err and "--location" in err
    data = json.loads(out)
    assert data["config_snapshot"]["impact"]["location_xy_mm"] is None
    assert data == _json(["--energy", "10"], capsys)


@pytest.mark.parametrize("bad", [["nan", "10"], ["10"], ["a", "b"]])
def test_location_rejects_bad_pairs(capsys, bad):
    code, out, err = _main(["--energy", "10", "--location", *bad], capsys)
    assert code == 2
    assert "--location" in err


@pytest.mark.parametrize("loc", [["200", "50"], ["0", "50"]])
def test_location_off_the_panel_is_a_clean_error(capsys, loc):
    code, out, err = _main(["--energy", "10", "--location", *loc, "--quick"], capsys)
    assert code == 1
    assert "bvidfe: error: impact location" in err and "inside the panel" in err
    assert "Traceback" not in err
    assert out == ""


def test_mesh_options_reach_mesh_params(capsys, captured_configs):
    code, out, err = _main(
        [
            "--energy",
            "10",
            "--tier",
            "fe3d",
            "--in-plane-size",
            "2.5",
            "--elements-per-ply",
            "2",
            "--quick",
        ],
        capsys,
    )
    assert code == 0, err
    (cfg,) = captured_configs
    assert cfg.tier == "fe3d"
    assert cfg.mesh.in_plane_size_mm == 2.5
    assert cfg.mesh.elements_per_ply == 2
    assert cfg.mesh.cohesive_zone_factor == 1.0


def test_one_mesh_option_keeps_the_other_default(capsys, captured_configs):
    code, out, err = _main(
        ["--energy", "10", "--tier", "empirical,fe3d", "--elements-per-ply", "3", "--quick"],
        capsys,
    )
    assert code == 0, err
    fe3d = [c for c in captured_configs if c.tier == "fe3d"]
    assert len(fe3d) == 1
    assert fe3d[0].mesh.elements_per_ply == 3
    assert fe3d[0].mesh.in_plane_size_mm == 5.0


def test_fe3d_without_mesh_options_uses_mesh_params_defaults(capsys, captured_configs):
    code, out, err = _main(["--energy", "10", "--tier", "fe3d", "--quick"], capsys)
    assert code == 0, err
    from bvidfe.analysis import MeshParams

    assert captured_configs[0].mesh == MeshParams()


@pytest.mark.parametrize(
    "flag, value",
    [("--elements-per-ply", "0"), ("--elements-per-ply", "1.5"), ("--in-plane-size", "0")],
)
def test_mesh_options_reject_bad_values(capsys, flag, value):
    code, out, err = _main(["--energy", "10", "--tier", "fe3d", flag, value], capsys)
    assert code == 2
    assert flag in err


@pytest.mark.parametrize("opt", [["--in-plane-size", "2.5"], ["--elements-per-ply", "2"]])
def test_mesh_options_without_fe3d_are_an_error(capsys, opt):
    code, out, err = _main(["--energy", "10", "--tier", "empirical", *opt], capsys)
    assert code == 2
    assert opt[0] in err and "fe3d" in err
    assert out == ""


# ------------------------------------------- impact-only flags + --cscan


@pytest.mark.parametrize(
    "opt",
    [
        ["--mass", "2"],
        ["--impactor-diameter", "8"],
        ["--impactor-shape", "flat"],
        ["--location", "50", "40"],
    ],
)
def test_impact_only_options_with_cscan_are_an_error(capsys, opt):
    code, out, err = _main(["--cscan", SAMPLE_CSCAN, *opt, "--quick"], capsys)
    assert code == 2
    assert opt[0] in err and "--cscan" in err
    assert out == ""


def test_cscan_alone_still_runs(capsys):
    code, out, err = _main(["--cscan", SAMPLE_CSCAN, "--quick"], capsys)
    assert code == 0, err
    assert 0 < float(out) <= 1


def test_impact_option_defaults_are_the_api_defaults(capsys):
    from bvidfe.core.geometry import ImpactorGeometry
    from bvidfe.impact.mapping import ImpactEvent

    data = _json(["--energy", "10"], capsys)
    impact = data["config_snapshot"]["impact"]
    assert impact["mass_kg"] == ImpactEvent(energy_J=1.0).mass_kg
    assert impact["impactor"] == {
        "diameter_mm": ImpactorGeometry().diameter_mm,
        "shape": ImpactorGeometry().shape,
    }


# ------------------------------------------------------------- ranges


def test_layup_accepts_any_finite_angle():
    assert cli._parse_layup("0,135,-135,270") == [0.0, 135.0, -135.0, 270.0]


def test_layup_135_is_minus_45(capsys):
    """Ply angles are taken modulo 180 deg, as in the API."""
    a = _json(["--energy", "15", "--quick-json"], capsys)
    code, out, err = _main(
        ["--layup", "0,45,135,90,90,135,45,0", "--energy", "15", "--quick-json"], capsys
    )
    assert code == 0, err
    assert json.loads(out)["knockdown"] == pytest.approx(a["knockdown"], rel=1e-9)


def test_energy_zero_matches_impact_event(capsys):
    """``ImpactEvent(energy_J=0)`` is accepted and gives no damage."""
    data = _json(["--energy", "0"], capsys)
    assert data["knockdown"] == 1.0
    assert data["dpa_mm2"] == 0.0
    assert data["damage"]["delaminations"] == []


@pytest.mark.parametrize("bad, match", [("-1", ">= 0"), ("nan", "finite"), ("inf", "finite")])
def test_energy_rejects_negative_and_non_finite(capsys, bad, match):
    code, out, err = _main(["--energy", bad], capsys)
    assert code == 2
    assert match in err


# ------------------------------------------------- transverse stress ratio


def test_transverse_stress_ratio_defaults_to_uniaxial(capsys):
    data = _json(["--energy", "10"], capsys)
    assert data["config_snapshot"]["transverse_stress_ratio"] == 0.0
    assert data["biaxial_interaction_factor"] is None
    explicit = _json(["--energy", "10", "--transverse-stress-ratio", "0"], capsys)
    assert explicit == data


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_transverse_stress_ratio_reaches_the_config(loading, capsys):
    data = _json(
        ["--energy", "10", "--loading", loading, "--transverse-stress-ratio", "-0.5"], capsys
    )
    assert data["config_snapshot"]["transverse_stress_ratio"] == -0.5
    assert data["biaxial_interaction_factor"] is not None
    assert data["notes"][0].startswith("Biaxial sigma_y/sigma_x = -0.50: experimental")
    quick = _main(
        ["--energy", "10", "--loading", loading, "--transverse-stress-ratio", "-0.5", "--quick"],
        capsys,
    )
    assert quick[0] == 0 and float(quick[1]) == pytest.approx(data["knockdown"], abs=1e-6)


@pytest.mark.parametrize("value", ["1.5", "-2", "nan", "inf", "x"])
def test_transverse_stress_ratio_out_of_range_is_a_usage_error(value, capsys):
    code, _, err = _main(["--energy", "10", "--transverse-stress-ratio", value], capsys)
    assert code == 2
    assert "--transverse-stress-ratio" in err


def test_transverse_stress_ratio_help_says_experimental():
    text = " ".join(cli._build_parser().format_help().split())
    assert "--transverse-stress-ratio" in text and "EXPERIMENTAL" in text


@pytest.mark.parametrize("tier", ["fe3d", "empirical,fe3d"])
def test_transverse_stress_ratio_on_an_unsupported_tier_is_an_error(tier, capsys):
    code, out, err = _main(
        ["--energy", "10", "--tier", tier, "--transverse-stress-ratio", "0.5"], capsys
    )
    assert code == 1 and out == ""
    assert "not yet supported" in err and "#115" in err


def test_transverse_stress_ratio_runs_on_semi_analytical(capsys):
    data = _json(
        ["--energy", "10", "--tier", "semi_analytical", "--transverse-stress-ratio", "0.5"], capsys
    )
    assert data["tier_used"] == "semi_analytical"
    assert data["biaxial_interaction_factor"] is not None
    assert "semi_analytical puts the film's lateral stress" in data["notes"][0]
