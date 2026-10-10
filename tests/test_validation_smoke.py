"""Smoke tests for the validation harness."""

import subprocess
import sys
from pathlib import Path


def test_validator_module_imports():
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location(
        "validate_bvid_public",
        Path(__file__).parent.parent / "validation" / "validate_bvid_public.py",
    )
    mod = importlib.util.module_from_spec(spec)
    # Register in sys.modules before exec so @dataclass can resolve cls.__module__
    # on Python 3.9 (needed when loading a file outside any package).
    sys.modules["validate_bvid_public"] = mod
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop("validate_bvid_public", None)
    assert hasattr(mod, "DatasetCase")
    assert hasattr(mod, "run_dataset")
    assert hasattr(mod, "main")


def test_validator_runs_synthetic_dataset():
    res = subprocess.run(
        [
            sys.executable,
            "validation/validate_bvid_public.py",
            "--dataset",
            "synthetic_selfcheck",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, res.stderr
    assert "synthetic_selfcheck" in res.stdout
    assert "MAE" in res.stdout


def test_validator_gate_passes_on_synthetic_dataset():
    res = subprocess.run(
        [
            sys.executable,
            "validation/validate_bvid_public.py",
            "--dataset",
            "synthetic_selfcheck",
            "--gate",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, res.stderr


def _load_validator():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "validate_bvid_public",
        Path(__file__).parent.parent / "validation" / "validate_bvid_public.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["validate_bvid_public"] = mod
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop("validate_bvid_public", None)
    return mod


_CASE = {
    "material": "AS4/8552",
    "layup_deg": [45, 0, -45, 90, 90, -45, 0, 45],
    "ply_thickness_mm": 0.188,
    "panel_Lx_mm": 150.0,
    "panel_Ly_mm": 100.0,
    "impactor_diameter_mm": 16.0,
    "impactor_mass_kg": 5.5,
    "impact_energy_J": 10.0,
    "measured_strength_MPa": 200.0,
}


def test_case_requires_impactor_diameter_and_mass():
    import pytest

    v = _load_validator()
    for key in ("impactor_diameter_mm", "impactor_mass_kg"):
        d = {k: val for k, val in _CASE.items() if k != key}
        with pytest.raises(ValueError, match=key):
            v.case_from_dict(d)


def test_case_boundary_and_impactor_shape_reach_the_analysis(monkeypatch):
    v = _load_validator()
    seen = {}
    real = v.BvidAnalysis

    class Spy(real):
        def __init__(self, config):
            seen["config"] = config
            super().__init__(config)

    monkeypatch.setattr(v, "BvidAnalysis", Spy)
    v.run_case(v.case_from_dict(_CASE), "empirical")
    assert seen["config"].panel.boundary == "simply_supported"
    assert seen["config"].impact.impactor.shape == "hemispherical"

    case = v.case_from_dict({**_CASE, "boundary": "clamped", "impactor_shape": "flat"})
    v.run_case(case, "empirical")
    assert seen["config"].panel.boundary == "clamped"
    assert seen["config"].impact.impactor.shape == "flat"


def test_ncamp_as4_8552_cai_gate_passes_for_the_closed_form_tiers():
    for tier in ("empirical", "semi_analytical"):
        res = subprocess.run(
            [
                sys.executable,
                "validation/validate_bvid_public.py",
                "--dataset",
                "ncamp_as4_8552_cai",
                "--tier",
                tier,
                "--gate",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0, (tier, res.stdout, res.stderr)
        assert "ncamp_as4_8552_cai" in res.stdout


def test_ungated_dataset_reports_its_error_without_failing_the_gate(tmp_path, monkeypatch, capsys):
    import json

    v = _load_validator()
    wrong = {**_CASE, "measured_strength_MPa": 1.0}  # any prediction is far off
    for name, extra in (("gated", {}), ("advisory", {"gate": False, "gate_note": "known miss"})):
        (tmp_path / f"{name}.json").write_text(
            json.dumps({"name": name, "target_mae_pct": 10.0, "cases": [wrong], **extra})
        )
    monkeypatch.setattr(v, "DATASET_DIR", tmp_path)

    assert v.main(["--dataset", "advisory", "--gate"]) == 0
    assert "ADVISORY: advisory is not gated. known miss" in capsys.readouterr().out
    assert v.main(["--dataset", "gated", "--gate"]) == 1


def test_lovejoy_scotti_dataset_is_advisory():
    res = subprocess.run(
        [
            sys.executable,
            "validation/validate_bvid_public.py",
            "--dataset",
            "lovejoy_scotti_im7_8552_cai",
            "--gate",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, res.stderr
    assert "ADVISORY: lovejoy_scotti_im7_8552_cai is not gated" in res.stdout
