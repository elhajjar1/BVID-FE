"""Smoke tests for the validation harness."""

import subprocess
import sys
from pathlib import Path

import pytest


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


@pytest.mark.parametrize("dataset", ["ncamp_as4_8552_cai", "lovejoy_scotti_im7_8552_cai"])
def test_measured_cai_gates_pass_for_the_closed_form_tiers(dataset):
    for tier in ("empirical", "semi_analytical"):
        res = subprocess.run(
            [
                sys.executable,
                "validation/validate_bvid_public.py",
                "--dataset",
                dataset,
                "--tier",
                tier,
                "--gate",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0, (tier, res.stdout, res.stderr)
        assert dataset in res.stdout
        assert "ADVISORY" not in res.stdout


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


def test_validator_runs_several_named_datasets():
    v = _load_validator()
    paths = v._select_datasets(["ncamp_as4_8552_cai", "lovejoy_scotti_im7_8552_cai"])
    assert [p.stem for p in paths] == ["lovejoy_scotti_im7_8552_cai", "ncamp_as4_8552_cai"]


def test_unknown_dataset_name_is_an_error(capsys):
    v = _load_validator()
    assert v.main(["--dataset", "ncamp_as4_8552_cai", "no_such_dataset"]) == 1
    assert "no_such_dataset" in capsys.readouterr().err


def test_measured_pristine_scores_the_knockdown_on_the_measured_scale():
    v = _load_validator()
    plain = v.run_case(v.case_from_dict(_CASE), "empirical")
    anchored = v.run_case(v.case_from_dict({**_CASE, "measured_pristine_MPa": 300.0}), "empirical")
    assert anchored["predicted_knockdown"] == pytest.approx(plain["predicted_knockdown"])
    assert anchored["predicted_MPa"] == pytest.approx(300.0 * plain["predicted_knockdown"])
    assert anchored["error_pct"] == pytest.approx(
        100 * abs(anchored["predicted_MPa"] - 200.0) / 200.0
    )


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_case_impact_panel_reaches_the_analysis(monkeypatch):
    v = _load_validator()
    seen = {}
    real = v.BvidAnalysis

    class Spy(real):
        def __init__(self, config):
            seen["config"] = config
            super().__init__(config)

    monkeypatch.setattr(v, "BvidAnalysis", Spy)
    v.run_case(v.case_from_dict(_CASE), "empirical")
    assert seen["config"].impact_panel is None

    case = {**_CASE, "impact_panel_Lx_mm": 60.0, "impact_panel_Ly_mm": 60.0}
    v.run_case(v.case_from_dict({**case, "impact_boundary": "clamped"}), "empirical")
    p = seen["config"].impact_panel
    assert (p.Lx_mm, p.Ly_mm, p.boundary) == (60.0, 60.0, "clamped")
    v.run_case(v.case_from_dict(case), "empirical")  # boundary defaults to the test panel's
    assert seen["config"].impact_panel.boundary == "simply_supported"


def test_hasebe_dataset_is_advisory_and_runs():
    v = _load_validator()
    name, cases, meta = v.load_dataset(v.DATASET_DIR / "hasebe_t800s_3900_2b_cai.json")
    assert not meta["gate"] and meta["gate_note"]
    assert len(cases) == 343
    assert {c.material for c in cases} == {"T800S/3900-2B"}
    assert all(c.measured_pristine_MPa and c.impact_panel_Lx_mm == 60.0 for c in cases)
    df = v.run_dataset(cases[:5], "empirical")
    assert df["error_pct"].notna().all()


@pytest.mark.parametrize("dataset", ["ncamp_as4_8552_cai", "lovejoy_scotti_im7_8552_cai"])
def test_astm_datasets_impact_on_the_d7136_window(dataset):
    v = _load_validator()
    _, cases, _ = v.load_dataset(v.DATASET_DIR / f"{dataset}.json")
    for c in cases:
        assert (c.impact_panel_Lx_mm, c.impact_panel_Ly_mm) == (125.0, 75.0)
        assert c.impact_boundary == "simply_supported"


def test_dent_summary_is_reported_not_gated(capsys):
    """The validator reports the median predicted / measured dent over the
    cases with a positive measured dent; it does not enter the gate."""
    import pandas as pd

    v = _load_validator()
    df = pd.DataFrame(
        {
            "predicted_dent_mm": [0.2, 0.3, 0.4, 0.1],
            "measured_dent_mm": [0.1, 0.3, 0.8, 0.0],
        }
    )
    s = v.dent_summary(df)
    assert (s["n"], s["n_nonpositive"]) == (3, 1)
    assert s["median"] == pytest.approx(1.0)
    no_dent = pd.DataFrame({"predicted_dent_mm": [0.2], "measured_dent_mm": [None]})
    assert v.dent_summary(no_dent) is None

    row = v.run_case(v.case_from_dict({**_CASE, "measured_dent_mm": 0.25}), "empirical")
    assert row["measured_dent_mm"] == 0.25 and row["predicted_dent_mm"] >= 0.0

    assert v.main(["--dataset", "lovejoy_scotti_im7_8552_cai", "--gate"]) == 0
    out = capsys.readouterr().out
    assert "dent (reported, not gated): median predicted / measured = 1.00" in out
