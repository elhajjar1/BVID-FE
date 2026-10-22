"""Input validation, warnings capture, CLI error handling and logging setup.

Non-finite and out-of-domain inputs used to pass construction (``NaN <= 0``
is False) and fail later, or not at all; the CLI printed raw tracebacks;
``AnalysisResults.warnings`` was always empty; and importing the library
attached a stderr log handler.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import warnings

import pytest

from bvidfe._types import _TIER_NAMES
from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.analysis.bvid import DamageOutsidePanelWarning
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import DPACapClipWarning, ImpactEvent, impact_to_damage

NAN = float("nan")
INF = float("inf")
LAYUP = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]


@pytest.fixture(autouse=True)
def _restore_fe3d_logger():
    """The in-process CLI calls attach a stderr handler; undo it per test."""
    log = logging.getLogger("bvidfe.fe3d")
    handlers, level = list(log.handlers), log.level
    yield
    for h in list(log.handlers):
        if h not in handlers:
            log.removeHandler(h)
    log.setLevel(level)


def _cfg(**kw):
    base = dict(
        material="IM7/8552",
        layup_deg=list(LAYUP),
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        impact=ImpactEvent(energy_J=20.0),
    )
    base.update(kw)
    return AnalysisConfig(**base)


# ---------------------------------------------------------------- dataclasses


@pytest.mark.parametrize("bad", [NAN, INF, -INF, 0.0, -1.0, "150", True])
def test_panel_geometry_rejects_bad_dimensions(bad):
    with pytest.raises(ValueError, match="Lx_mm"):
        PanelGeometry(bad, 100.0)
    with pytest.raises(ValueError, match="Ly_mm"):
        PanelGeometry(100.0, bad)


@pytest.mark.parametrize("bad", [NAN, INF, 0.0, -16.0])
def test_impactor_geometry_rejects_bad_diameter(bad):
    with pytest.raises(ValueError, match="diameter_mm"):
        ImpactorGeometry(diameter_mm=bad)


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"major_mm": NAN}, "major_mm"),
        ({"minor_mm": INF}, "minor_mm"),
        ({"major_mm": -1.0}, "positive"),
        ({"orientation_deg": NAN}, "orientation_deg"),
        ({"centroid_mm": (NAN, 0.0)}, r"centroid_mm\[0\]"),
        ({"centroid_mm": (0.0, 1.0, 2.0)}, "centroid_mm"),
        ({"interface_index": True}, "interface_index"),
        ({"interface_index": -1}, "interface_index"),
    ],
)
def test_delamination_ellipse_rejects_bad_fields(kw, match):
    args = dict(
        interface_index=0, centroid_mm=(75.0, 50.0), major_mm=10.0, minor_mm=5.0, orientation_deg=0
    )
    args.update(kw)
    with pytest.raises(ValueError, match=match):
        DelaminationEllipse(**args)


def test_delamination_centroid_list_becomes_tuple():
    e = DelaminationEllipse(0, [75.0, 50.0], 10.0, 5.0, 0.0)
    assert e.centroid_mm == (75.0, 50.0)


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"dent_depth_mm": -1.0}, "dent_depth_mm"),
        ({"dent_depth_mm": NAN}, "dent_depth_mm"),
        ({"fiber_break_radius_mm": -5.0}, "fiber_break_radius_mm"),
        ({"fiber_break_radius_mm": INF}, "fiber_break_radius_mm"),
        ({"delaminations": [(0, (1, 1), 1, 1, 0)]}, r"delaminations\[0\]"),
    ],
)
def test_damage_state_rejects_bad_fields(kw, match):
    with pytest.raises(ValueError, match=match):
        DamageState(**kw)


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"elements_per_ply": True}, "elements_per_ply"),
        ({"elements_per_ply": 0}, "elements_per_ply"),
        ({"in_plane_size_mm": INF}, "in_plane_size_mm"),
        ({"in_plane_size_mm": NAN}, "in_plane_size_mm"),
        ({"cohesive_zone_factor": INF}, "cohesive_zone_factor"),
    ],
)
def test_mesh_params_rejects_bad_fields(kw, match):
    with pytest.raises(ValueError, match=match):
        MeshParams(**kw)


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"energy_J": -5.0}, "energy_J"),
        ({"energy_J": NAN}, "energy_J"),
        ({"energy_J": INF}, "energy_J"),
        ({"mass_kg": 0.0}, "mass_kg"),
        ({"mass_kg": -1.0}, "mass_kg"),
        ({"mass_kg": NAN}, "mass_kg"),
        ({"location_xy_mm": (NAN, 0.0)}, "location_xy_mm"),
        ({"location_xy_mm": (1.0,)}, "location_xy_mm"),
        ({"impactor": 16.0}, "impactor"),
    ],
)
def test_impact_event_rejects_bad_fields(kw, match):
    args = {"energy_J": 20.0}
    args.update(kw)
    with pytest.raises(ValueError, match=match):
        ImpactEvent(**args)


def test_impact_event_zero_energy_is_allowed():
    assert ImpactEvent(energy_J=0.0).energy_J == 0.0


def test_impact_event_list_location_keeps_panel_centre_sentinel():
    """JSON gives ``[0.0, 0.0]``; it must mean the panel centre like the tuple."""
    ev_list = ImpactEvent(energy_J=20.0, location_xy_mm=[0.0, 0.0])
    assert ev_list.location_xy_mm == (0.0, 0.0)
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], LAYUP, 0.152)
    panel = PanelGeometry(150.0, 100.0)
    a = impact_to_damage(ev_list, lam, panel)
    b = impact_to_damage(ImpactEvent(energy_J=20.0), lam, panel)
    assert a.projected_damage_area_mm2 == b.projected_damage_area_mm2
    assert [e.centroid_mm for e in a.delaminations] == [e.centroid_mm for e in b.delaminations]


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"ply_thickness_mm": INF}, "ply_thickness_mm"),
        ({"ply_thickness_mm": NAN}, "ply_thickness_mm"),
        ({"ply_thickness_mm": [0.152] * 7 + [INF]}, r"ply_thickness_mm\[7\]"),
        ({"layup_deg": []}, "at least one ply"),
        ({"layup_deg": "0,45"}, "not a string"),
        ({"layup_deg": [0.0, NAN, 90.0]}, r"layup_deg\[1\]"),
        ({"layup_deg": [0.0, INF]}, r"layup_deg\[1\]"),
        ({"layup_deg": 45.0}, "layup_deg"),
        ({"panel": (150.0, 100.0)}, "PanelGeometry"),
        ({"material": 42}, "material"),
        ({"mesh": {"in_plane_size_mm": 5.0}}, "MeshParams"),
    ],
)
def test_analysis_config_rejects_bad_fields(kw, match):
    with pytest.raises(ValueError, match=match):
        _cfg(**kw)


@pytest.mark.parametrize("name", ["im7/8552", "IM7-8552", "nope"])
def test_analysis_config_rejects_unknown_material_listing_valid_names(name):
    with pytest.raises(ValueError, match="MATERIAL_LIBRARY") as exc:
        _cfg(material=name)
    assert "IM7/8552" in str(exc.value)
    assert "case-sensitive" in str(exc.value)


def test_valid_config_still_constructs():
    cfg = _cfg(layup_deg=tuple(LAYUP), ply_thickness_mm=[0.152] * 8)
    assert cfg.panel.Lx_mm == 150.0


# ------------------------------------------------- damage outside the panel


def _damage_cfg(centroid):
    return _cfg(
        impact=None,
        damage=DamageState([DelaminationEllipse(3, centroid, 10.0, 6.0, 0.0)], 0.1),
    )


@pytest.mark.parametrize("centroid", [(500.0, 50.0), (75.0, -20.0), (-1.0, -1.0)])
def test_delamination_outside_panel_warns_and_is_recorded(centroid):
    with pytest.warns(DamageOutsidePanelWarning, match="outside the 150 x 100 mm panel"):
        r = BvidAnalysis(_damage_cfg(centroid)).run()
    assert any(w.startswith("DamageOutsidePanelWarning:") for w in r.warnings)


@pytest.mark.parametrize("centroid", [(75.0, 50.0), (0.0, 0.0), (150.0, 100.0)])
def test_delamination_inside_panel_does_not_warn(centroid):
    with warnings.catch_warnings():
        warnings.simplefilter("error", DamageOutsidePanelWarning)
        r = BvidAnalysis(_damage_cfg(centroid)).run()
    assert r.warnings == []


# --------------------------------------------------- results.warnings capture


def _clipping_cfg():
    # 1500 mm^2 of predicted damage on a 30 x 30 mm panel hits the 80% cap.
    return _cfg(panel=PanelGeometry(30.0, 30.0), impact=ImpactEvent(energy_J=30.0))


def test_run_records_warnings_and_still_emits_them():
    with pytest.warns(DPACapClipWarning):
        r = BvidAnalysis(_clipping_cfg()).run()
    clip = [w for w in r.warnings if w.startswith("DPACapClipWarning: ")]
    assert len(clip) == 1 and "80% of panel" in clip[0]
    assert r.to_dict()["warnings"] == r.warnings
    assert len(r.warnings) == len(set(r.warnings))


def test_run_records_warnings_even_when_the_caller_ignores_them():
    with warnings.catch_warnings(record=True) as shown:
        warnings.simplefilter("ignore")
        r = BvidAnalysis(_clipping_cfg()).run()
    assert shown == []
    assert any(w.startswith("DPACapClipWarning") for w in r.warnings)


def test_run_still_raises_under_an_error_filter():
    with warnings.catch_warnings():
        warnings.simplefilter("error", DPACapClipWarning)
        with pytest.raises(DPACapClipWarning):
            BvidAnalysis(_clipping_cfg()).run()


def test_quiet_run_has_no_warnings():
    assert BvidAnalysis(_cfg()).run().warnings == []


# ---------------------------------------------------------------------- CLI

CLI_BASE = [
    "--material",
    "IM7/8552",
    "--layup",
    "0,45,-45,90,90,-45,45,0",
    "--thickness",
    "0.152",
    "--loading",
    "compression",
]


def _main(args, capsys):
    from bvidfe.cli import main

    try:
        code = main(CLI_BASE + args)
    except SystemExit as exc:
        code = exc.code
    out, err = capsys.readouterr()
    return code, out, err


@pytest.mark.parametrize(
    "spec, match",
    [("inf", "finite"), ("nan", "finite"), ("-inf", "finite"), ("0", "> 0")],
)
def test_positive_float_rejects_non_finite(spec, match):
    from bvidfe.cli import _positive_float

    with pytest.raises(argparse.ArgumentTypeError, match=match):
        _positive_float(spec)


@pytest.mark.parametrize(
    "spec, match",
    [("0,nan,90", "finite"), ("0,inf", "finite"), ("0,135", r"\[-90, 90\]"), ("-91", "-90")],
)
def test_parse_layup_rejects_non_finite_and_out_of_range(spec, match):
    from bvidfe.cli import _parse_layup

    with pytest.raises(argparse.ArgumentTypeError, match=match):
        _parse_layup(spec)


def test_parse_layup_accepts_the_range_ends():
    from bvidfe.cli import _parse_layup

    assert _parse_layup("-90,90") == [-90.0, 90.0]


@pytest.mark.parametrize("spec", ["0.1,inf", "0.1,nan", "inf"])
def test_parse_thickness_rejects_non_finite(spec):
    from bvidfe.cli import _parse_thickness

    with pytest.raises(argparse.ArgumentTypeError, match="finite"):
        _parse_thickness(spec)


@pytest.mark.parametrize("spec", ["infx100", "100xnan"])
def test_parse_panel_rejects_non_finite(spec):
    from bvidfe.cli import _parse_panel

    with pytest.raises(argparse.ArgumentTypeError, match="finite"):
        _parse_panel(spec)


def test_cli_tier_list_is_the_shared_tier_set():
    from bvidfe.cli import _VALID_TIERS

    assert set(_VALID_TIERS) == _TIER_NAMES


@pytest.mark.parametrize(
    "content, match",
    [
        ('{"dent_depth_mm": 0, "delaminations": []}', "schema_version"),
        ("{not json", "invalid JSON"),
    ],
)
def test_cli_bad_cscan_is_a_clean_error(tmp_path, capsys, content, match):
    p = tmp_path / "scan.json"
    p.write_text(content)
    code, out, err = _main(["--panel", "150x100", "--cscan", str(p)], capsys)
    assert code == 2
    assert "invalid C-scan file" in err and match in err
    assert "Traceback" not in err


def test_cli_config_value_error_is_a_clean_error(tmp_path, capsys):
    """A C-scan interface beyond the 8-ply layup fails AnalysisConfig."""
    p = tmp_path / "scan.json"
    p.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "dent_depth_mm": 0.1,
                "delaminations": [
                    {
                        "interface_index": 12,
                        "centroid_mm": [75, 50],
                        "major_mm": 10,
                        "minor_mm": 5,
                        "orientation_deg": 0,
                    }
                ],
            }
        )
    )
    code, out, err = _main(["--panel", "150x100", "--cscan", str(p)], capsys)
    assert code == 2
    assert "interface_index 12" in err
    assert "Traceback" not in err


def test_cli_fe3d_size_error_is_a_clean_error(capsys):
    code, out, err = _main(
        ["--panel", "1500x1000", "--energy", "20", "--tier", "fe3d", "--quick"], capsys
    )
    assert code == 1
    assert "bvidfe: error: fe3d problem too large" in err
    assert out == ""


def test_cli_json_includes_run_warnings(capsys):
    code, out, err = _main(["--panel", "30x30", "--energy", "30"], capsys)
    assert code == 0
    data = json.loads(out)
    assert any(w.startswith("DPACapClipWarning") for w in data["warnings"])


def test_cli_help_shows_defaults_once_and_never_none(capsys):
    from bvidfe.cli import _build_parser

    text = _build_parser().format_help()
    flat = " ".join(text.split())
    assert "(default: None)" not in flat
    assert "(default: False)" not in flat
    assert "Impactor diameter in mm (default: 16.0)" in flat
    assert "Impactor mass in kg (default: 5.5)" in flat
    assert "(default 16.0)" not in flat and "(default 5.5)" not in flat


def test_cli_closed_stdout_pipe_exits_without_traceback():
    proc = subprocess.Popen(
        [sys.executable, "-m", "bvidfe.cli", *CLI_BASE, "--panel", "150x100", "--energy", "20"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    proc.stdout.close()  # the reader goes away before any output
    err = proc.stderr.read().decode()
    proc.stderr.close()
    code = proc.wait(timeout=120)
    assert "Traceback" not in err and "BrokenPipeError" not in err, err
    assert code == 1


# ------------------------------------------------------------------ logging


def test_library_import_adds_only_a_null_handler():
    code = (
        "import logging, bvidfe, bvidfe.analysis.fe_tier;"
        "lg = logging.getLogger('bvidfe.fe3d');"
        "print([type(h).__name__ for h in lg.handlers], lg.level)"
    )
    env_free = {k: v for k, v in os.environ.items() if k != "BVIDFE_LOG_LEVEL"}
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=env_free, check=True
    )
    assert out.stdout.strip() == "['NullHandler'] 0"
    assert out.stderr == ""


def test_enable_stderr_logging_is_idempotent_and_follows_stderr(capsys):
    from bvidfe.analysis import fe_tier

    log = logging.getLogger("bvidfe.fe3d")
    for h in [h for h in log.handlers if isinstance(h, fe_tier._StderrHandler)]:
        log.removeHandler(h)  # restored by the autouse fixture
    fe_tier.enable_stderr_logging()
    fe_tier.enable_stderr_logging()
    assert sum(isinstance(h, fe_tier._StderrHandler) for h in log.handlers) == 1
    log.warning("hello from fe3d")
    assert "hello from fe3d" in capsys.readouterr().err
