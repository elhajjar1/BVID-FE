"""Smoke tests for validation/plot_validation.py (predicted-vs-measured plots)."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent


@pytest.fixture
def pv():
    spec = importlib.util.spec_from_file_location(
        "plot_validation", ROOT / "validation" / "plot_validation.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["plot_validation"] = mod
    try:
        spec.loader.exec_module(mod)
        yield mod
    finally:
        sys.modules.pop("plot_validation", None)
        sys.modules.pop("validate_bvid_public", None)


def test_ncamp_unnotched_points(pv):
    un = pv.ncamp_unnotched()
    assert len(un) == 12
    assert set(un["loading"]) == {"tension", "compression"}
    ratio = un["predicted_MPa"] / un["measured_MPa"]
    # The same bounds as tests/validation/test_unnotched_strength_ncamp.py.
    assert ratio.between(0.8, 1.2).all()


def test_limits_cover_every_point(pv):
    import pandas as pd

    df = pd.DataFrame({"a": [0.5, 2.0, None, -1.0], "b": [1.0, 3.0, 0.2, 0.0]})
    lo, hi = pv._limits(df, ("a", "b"), log=True)
    assert lo < 0.2 and hi > 3.0
    lo, hi = pv._limits(df, ("a", "b"), log=False)
    assert lo == 0.0 and hi > 3.0


def test_main_writes_every_panel(pv, tmp_path, monkeypatch):
    """End to end on the two smallest measured datasets (fast)."""
    v = pv.v
    keep = {"ncamp_as4_8552_cai", "korbelin_m21_t800s_tai"}
    real = v._select_datasets
    monkeypatch.setattr(
        v,
        "_select_datasets",
        lambda names, kind=None: [p for p in real(names, kind) if p.stem in keep],
    )
    assert pv.main(["--out", str(tmp_path / "out"), "--cache", str(tmp_path / "cache")]) == 0
    out = tmp_path / "out"
    for name in (
        "strength",
        "knockdown",
        "damage_area",
        "dent_depth",
        "unnotched_strength",
        "validation_overview",
    ):
        assert (out / f"{name}.png").stat().st_size > 10_000, name
    import pandas as pd

    points = pd.read_csv(out / "validation_points.csv")
    assert set(points["dataset"]) == keep
    assert {"predicted_MPa", "measured_MPa", "measured_dpa_mm2"} <= set(points.columns)
    # The cache is reused on a second run.
    assert sorted(p.name for p in (tmp_path / "cache").iterdir()) == [
        "empirical_korbelin_m21_t800s_tai.csv",
        "empirical_ncamp_as4_8552_cai.csv",
    ]
