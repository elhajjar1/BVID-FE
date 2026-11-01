import warnings

import pandas as pd
import pytest

from bvidfe.analysis import AnalysisConfig
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import DPACapClipWarning, ImpactEvent
from bvidfe.impact.olsson import onset_energy
from bvidfe.sweep.parametric_sweep import (
    sweep_energies,
    sweep_layups,
    sweep_thicknesses,
)


def _base_impact_cfg(**overrides):
    kw = dict(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90] * 4,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(10.0, ImpactorGeometry(), mass_kg=5.5),
    )
    kw.update(overrides)
    return AnalysisConfig(**kw)


def test_sweep_energies_returns_dataframe_with_expected_columns():
    cfg = _base_impact_cfg()
    df = sweep_energies(cfg, energies_J=[5, 10, 20, 30])
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    for col in ["energy_J", "knockdown", "residual_MPa", "dpa_mm2", "dent_mm"]:
        assert col in df.columns


def test_sweep_energies_writes_csv(tmp_path):
    cfg = _base_impact_cfg()
    csv_path = tmp_path / "sweep.csv"
    sweep_energies(cfg, energies_J=[5, 10], csv_path=csv_path)
    assert csv_path.exists()
    loaded = pd.read_csv(csv_path)
    assert len(loaded) == 2


def test_sweep_layups():
    cfg = _base_impact_cfg()
    layups = [
        [0, 45, -45, 90] * 4,
        [0, 90] * 8,
        [0, 60, -60] * 5 + [0],
    ]
    df = sweep_layups(cfg, layups=layups)
    assert len(df) == 3
    assert "layup" in df.columns
    assert "knockdown" in df.columns


def test_sweep_thicknesses():
    cfg = _base_impact_cfg()
    df = sweep_thicknesses(cfg, ply_thicknesses_mm=[0.125, 0.152, 0.2])
    assert len(df) == 3
    assert "ply_thickness_mm" in df.columns


def test_sweep_energies_default_raises_on_iteration_failure(monkeypatch):
    """Issue #8 baseline: with on_error='raise' (default) a failed iteration
    aborts the sweep — preserving the legacy behaviour for callers that want
    a hard fail."""
    cfg = _base_impact_cfg()
    from bvidfe.sweep import parametric_sweep as ps

    real_run_one = ps._run_one

    def _patched(cfg):
        if cfg.impact.energy_J == 10.0:
            raise RuntimeError("synthetic failure at 10J")
        return real_run_one(cfg)

    monkeypatch.setattr(ps, "_run_one", _patched)
    import pytest

    with pytest.raises(RuntimeError, match="synthetic failure"):
        sweep_energies(cfg, energies_J=[5, 10, 20])


def test_sweep_energies_skip_preserves_partial_results(monkeypatch, tmp_path):
    """Issue #8: with on_error='skip' a per-iteration failure must not
    abort the whole sweep — the failed energy is recorded with NaN numerics
    and an 'error' column, the surrounding successful iterations land in
    the DataFrame, and the CSV is still written."""
    import math

    cfg = _base_impact_cfg()
    from bvidfe.sweep import parametric_sweep as ps

    real_run_one = ps._run_one

    def _patched(cfg):
        if cfg.impact.energy_J == 10.0:
            raise RuntimeError("synthetic failure at 10J")
        return real_run_one(cfg)

    monkeypatch.setattr(ps, "_run_one", _patched)
    csv_path = tmp_path / "partial.csv"
    df = sweep_energies(cfg, energies_J=[5, 10, 20], csv_path=csv_path, on_error="skip")
    assert len(df) == 3
    bad = df.loc[df["energy_J"] == 10.0].iloc[0]
    assert math.isnan(bad["knockdown"])
    assert "error" in df.columns
    assert "synthetic failure" in str(bad["error"])
    # Surrounding iterations still produced valid knockdowns
    good = df.loc[df["energy_J"] != 10.0]
    assert all(0.0 < kd <= 1.0 for kd in good["knockdown"])
    # Partial CSV is preserved
    assert csv_path.exists()
    loaded = pd.read_csv(csv_path)
    assert len(loaded) == 3


def test_sweep_energies_warn_emits_userwarning(monkeypatch):
    """Issue #8: on_error='warn' is identical to 'skip' but adds a
    UserWarning per failure so callers running interactively see the cause."""
    import warnings as _warnings

    cfg = _base_impact_cfg()
    from bvidfe.sweep import parametric_sweep as ps

    real_run_one = ps._run_one

    def _patched(cfg):
        if cfg.impact.energy_J == 10.0:
            raise RuntimeError("synthetic warn-mode failure")
        return real_run_one(cfg)

    monkeypatch.setattr(ps, "_run_one", _patched)
    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        df = sweep_energies(cfg, energies_J=[5, 10, 20], on_error="warn")
    assert any("synthetic warn-mode failure" in str(w.message) for w in caught)
    assert len(df) == 3


def test_sweep_energies_progress_callback_is_invoked():
    """Issue #8: progress_callback receives (i_done, n_total) after each
    iteration, allowing the GUI worker / a long-running script to drive a
    progress bar without re-implementing the loop."""
    cfg = _base_impact_cfg()
    seen: list[tuple[int, int]] = []
    sweep_energies(
        cfg,
        energies_J=[5, 10, 20],
        progress_callback=lambda i, n: seen.append((i, n)),
    )
    assert seen == [(1, 3), (2, 3), (3, 3)]


def test_sweep_energies_progress_callback_three_arg_receives_inputs():
    """Issue #112: the Streamlit UI passes a 3-arg callback so it can
    render a live ``E=…J`` label. Existing 2-arg callbacks must keep
    working (covered by ``..._progress_callback_is_invoked`` above) —
    this test pins the 3-arg dispatch and the inputs-dict contract."""
    cfg = _base_impact_cfg()
    seen: list[tuple[int, int, dict]] = []
    sweep_energies(
        cfg,
        energies_J=[5, 20],
        progress_callback=lambda i, n, inp: seen.append((i, n, dict(inp))),
    )
    assert [(i, n) for i, n, _ in seen] == [(1, 2), (2, 2)]
    assert seen[0][2] == {"energy_J": 5.0}
    assert seen[1][2] == {"energy_J": 20.0}


def test_sweep_progress_callback_exception_does_not_break_sweep():
    """Issue #112: a buggy UI callback (e.g. a stale Streamlit handle on
    a torn-down session) must not abort the whole sweep. The exception
    is swallowed and logged at DEBUG; the DataFrame still finalises."""
    cfg = _base_impact_cfg()
    calls = {"n": 0}

    def _bad_cb(i, n, inputs):
        calls["n"] += 1
        raise RuntimeError("UI callback exploded")

    df = sweep_energies(cfg, energies_J=[5, 10], progress_callback=_bad_cb)
    assert calls["n"] == 2
    assert len(df) == 2
    # Sweep itself succeeded — knockdowns are populated, error column NaN/empty.
    assert df["knockdown"].notna().all()


def test_sweep_layups_progress_callback_inputs_carries_layup():
    """The 3-arg callback for ``sweep_layups`` carries the layup string."""
    cfg = _base_impact_cfg()
    seen: list[dict] = []
    sweep_layups(
        cfg,
        layups=[[0, 90, 0, 90], [45, -45, 45, -45]],
        progress_callback=lambda i, n, inp: seen.append(dict(inp)),
    )
    assert seen == [{"layup": "0/90/0/90"}, {"layup": "45/-45/45/-45"}]


def test_sweep_thicknesses_progress_callback_inputs_carries_thickness():
    """The 3-arg callback for ``sweep_thicknesses`` carries the thickness."""
    cfg = _base_impact_cfg()
    seen: list[dict] = []
    sweep_thicknesses(
        cfg,
        ply_thicknesses_mm=[0.125, 0.200],
        progress_callback=lambda i, n, inp: seen.append(dict(inp)),
    )
    assert seen == [{"ply_thickness_mm": 0.125}, {"ply_thickness_mm": 0.200}]


def test_sweep_energies_rejects_unknown_on_error():
    cfg = _base_impact_cfg()
    import pytest

    with pytest.raises(ValueError, match="on_error must be one of"):
        sweep_energies(cfg, energies_J=[5], on_error="bogus")


def test_sweep_column_schema_is_deterministic():
    """Issue #38: column order is fixed (swept var first, `error` always
    present) regardless of whether any iteration failed."""
    expected_results = [
        "knockdown",
        "residual_MPa",
        "pristine_MPa",
        "dpa_mm2",
        "dent_mm",
        "n_delaminations",
        "tier_used",
        "error",
    ]
    cfg = _base_impact_cfg()

    # All-success run with on_error='raise': `error` column still present.
    df_e = sweep_energies(cfg, energies_J=[5, 10, 20])
    assert list(df_e.columns) == ["energy_J", *expected_results]
    assert df_e["error"].isna().all()

    df_l = sweep_layups(cfg, layups=[[0, 90] * 8, [0, 45, -45, 90] * 4])
    assert list(df_l.columns) == ["layup", *expected_results]

    df_t = sweep_thicknesses(cfg, ply_thicknesses_mm=[0.125, 0.152])
    assert list(df_t.columns) == ["ply_thickness_mm", *expected_results]


def test_sweep_dedupes_dpa_cap_warning_to_once_per_sweep():
    """Issue #55: the DPA-cap diagnostic must fire once for the whole sweep,
    not once per point. A small panel + high energies makes every iteration
    exceed the 80% DPA clip; each warning's message embeds the per-point DPA
    value, so Python's 'once' filter (keyed on message text) cannot collapse
    them — _dedupe_impact_warnings() collapses by category instead."""
    cfg = _base_impact_cfg(panel=PanelGeometry(40, 30))
    energies = [80.0, 100.0, 120.0, 140.0, 160.0]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        df = sweep_energies(cfg, energies_J=energies)
    assert len(df) == len(energies)
    dpa_warnings = [w for w in caught if issubclass(w.category, DPACapClipWarning)]
    # Without per-category dedup this equals len(energies) (one per point).
    assert len(dpa_warnings) == 1


# --- Invalid sweep points honour on_error ---------------------------------
# A point whose config fails AnalysisConfig / ImpactEvent validation used to
# raise from the sweep loop itself, before on_error was consulted.

_RESULT_NUMERIC = ["knockdown", "residual_MPa", "pristine_MPa", "dpa_mm2", "dent_mm"]


def _assert_failed_row(row, match):
    import math

    for col in _RESULT_NUMERIC:
        assert math.isnan(row[col])
    assert row["n_delaminations"] == 0
    assert row["tier_used"] == ""
    assert str(row["error"]).startswith("ValueError: ")
    assert match in str(row["error"])


def _assert_good_rows(df):
    assert df["error"].isna().all()
    assert all(0.0 < kd <= 1.0 for kd in df["knockdown"])


@pytest.mark.parametrize("on_error", ["skip", "warn"])
def test_sweep_energies_invalid_point_becomes_failed_row(on_error):
    cfg = _base_impact_cfg()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        df = sweep_energies(cfg, energies_J=[5, -1, 20], on_error=on_error)
    assert df["energy_J"].tolist() == [5.0, -1.0, 20.0]
    _assert_failed_row(df.iloc[1], "energy_J")
    _assert_good_rows(df.iloc[[0, 2]])
    failures = [w for w in caught if "sweep iteration energy_J=-1 failed" in str(w.message)]
    assert len(failures) == (1 if on_error == "warn" else 0)


def test_sweep_energies_invalid_point_raises_by_default():
    with pytest.raises(ValueError, match="energy_J"):
        sweep_energies(_base_impact_cfg(), energies_J=[5, -1, 20])


def test_sweep_layups_thickness_list_mismatch_skips_point(tmp_path):
    # A 4-entry per-ply thickness list fits the 4-ply layups only.
    cfg = _base_impact_cfg(layup_deg=[0, 45, -45, 90], ply_thickness_mm=[0.152] * 4)
    csv_path = tmp_path / "layups.csv"
    df = sweep_layups(
        cfg,
        layups=[[0, 90, 90, 0], [0, 45, -45, 90] * 2, [45, -45, -45, 45]],
        csv_path=csv_path,
        on_error="skip",
    )
    assert df["layup"].tolist() == ["0/90/90/0", "0/45/-45/90/0/45/-45/90", "45/-45/-45/45"]
    _assert_failed_row(df.iloc[1], "ply_thickness_mm sequence length")
    _assert_good_rows(df.iloc[[0, 2]])
    assert len(pd.read_csv(csv_path)) == 3


def test_sweep_layups_interface_past_shorter_layup_skips_point():
    # A given damage state on interface 5 needs at least 7 plies.
    damage = DamageState([DelaminationEllipse(5, (75, 50), 20, 12, 45)])
    cfg = _base_impact_cfg(impact=None, damage=damage)
    df = sweep_layups(
        cfg, layups=[[0, 45, -45, 90] * 2, [0, 90, 90, 0], [0, 90] * 4], on_error="skip"
    )
    _assert_failed_row(df.iloc[1], "interface_index 5")
    assert df.loc[[0, 2], "error"].isna().all()
    assert df.loc[[0, 2], "knockdown"].notna().all()


def test_sweep_layups_invalid_point_raises_by_default():
    cfg = _base_impact_cfg(layup_deg=[0, 45, -45, 90], ply_thickness_mm=[0.152] * 4)
    with pytest.raises(ValueError, match="ply_thickness_mm sequence length"):
        sweep_layups(cfg, layups=[[0, 90, 90, 0], [0, 45, -45, 90] * 2])


def test_sweep_thicknesses_invalid_point_becomes_failed_row():
    df = sweep_thicknesses(
        _base_impact_cfg(), ply_thicknesses_mm=[0.125, 0.0, 0.2], on_error="skip"
    )
    assert df["ply_thickness_mm"].tolist() == [0.125, 0.0, 0.2]
    _assert_failed_row(df.iloc[1], "ply_thickness_mm")
    _assert_good_rows(df.iloc[[0, 2]])


def test_sweep_thicknesses_invalid_point_raises_by_default():
    with pytest.raises(ValueError, match="ply_thickness_mm"):
        sweep_thicknesses(_base_impact_cfg(), ply_thicknesses_mm=[0.125, 0.0])


@pytest.mark.parametrize(
    "sweep, kwargs",
    [
        (sweep_energies, {"energies_J": []}),
        (sweep_layups, {"layups": []}),
        (sweep_thicknesses, {"ply_thicknesses_mm": []}),
    ],
)
def test_sweep_rejects_unknown_on_error_with_empty_list(sweep, kwargs):
    with pytest.raises(ValueError, match="on_error must be one of"):
        sweep(_base_impact_cfg(), on_error="bogus", **kwargs)


def test_sweep_energies_fine_grid_across_onset_has_no_failed_rows():
    """Issue #163: a fine energy grid straddling the onset energy used to
    raise from distribute_damage a few mJ above onset. Every point runs,
    the damage area is zero at and below onset and non-decreasing above."""
    layup = [0, 45, -45, 90, 90, -45, 45, 0]
    cfg = _base_impact_cfg(layup_deg=layup)
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], layup, 0.152)
    E_onset = onset_energy(lam, cfg.panel, cfg.impact.impactor)
    energies = [E_onset + k * 1e-4 for k in range(-5, 6)] + [E_onset + 1e-3, E_onset + 1e-2]
    df = sweep_energies(cfg, energies_J=energies)
    assert len(df) == len(energies)
    assert df["dpa_mm2"].notna().all()
    assert df["knockdown"].notna().all()
    at_or_below = df["energy_J"] <= E_onset
    assert (df.loc[at_or_below, "dpa_mm2"] == 0.0).all()
    assert (df["dpa_mm2"].diff().dropna() >= 0.0).all()
    assert (df["knockdown"] <= 1.0).all()
