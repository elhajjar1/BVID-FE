"""``sweep_transverse_ratios``: the biaxial ratio sweep (issue #115, PR 4)."""

from __future__ import annotations

import math
from dataclasses import replace

import pandas as pd
import pytest

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.sweep import parametric_sweep as ps
from bvidfe.sweep.parametric_sweep import (
    sweep_energies,
    sweep_layups,
    sweep_thicknesses,
    sweep_transverse_ratios,
)

RATIOS = [-1.0, -0.5, 0.0, 0.5, 1.0]

_STANDARD_COLS = [
    "knockdown",
    "residual_MPa",
    "pristine_MPa",
    "dpa_mm2",
    "dent_mm",
    "n_delaminations",
    "tier_used",
    "error",
]
_RATIO_COLS = [
    "transverse_stress_ratio",
    *_STANDARD_COLS[:-1],
    "residual_sigma_y_MPa",
    "pristine_sigma_y_MPa",
    "biaxial_interaction_factor",
    "error",
]


def _cfg(**overrides) -> AnalysisConfig:
    kw: dict = dict(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(15.0, ImpactorGeometry(), mass_kg=5.5),
    )
    kw.update(overrides)
    return AnalysisConfig(**kw)


def test_columns_and_order():
    df = sweep_transverse_ratios(_cfg(), RATIOS)
    assert list(df.columns) == _RATIO_COLS
    assert df["transverse_stress_ratio"].tolist() == RATIOS
    assert df["error"].isna().all()
    assert set(df["tier_used"]) == {"empirical"}


def test_other_sweeps_columns_unchanged():
    cfg = _cfg()
    assert list(sweep_energies(cfg, [10.0]).columns) == ["energy_J", *_STANDARD_COLS]
    assert list(sweep_layups(cfg, [cfg.layup_deg]).columns) == ["layup", *_STANDARD_COLS]
    assert list(sweep_thicknesses(cfg, [0.152]).columns) == ["ply_thickness_mm", *_STANDARD_COLS]


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_rows_equal_looping_bvid_analysis(loading):
    """Each row is exactly what a ``BvidAnalysis`` run at that ratio gives."""
    base = _cfg(loading=loading)
    df = sweep_transverse_ratios(base, RATIOS)
    for beta, (_, row) in zip(RATIOS, df.iterrows()):
        r = BvidAnalysis(replace(base, transverse_stress_ratio=beta)).run()
        assert row["knockdown"] == r.knockdown
        assert row["residual_MPa"] == r.residual_strength_MPa
        assert row["pristine_MPa"] == r.pristine_strength_MPa
        assert row["dpa_mm2"] == r.dpa_mm2
        assert row["dent_mm"] == r.damage.dent_depth_mm
        assert row["n_delaminations"] == len(r.damage.delaminations)
        assert row["tier_used"] == r.tier_used
        if beta == 0:
            # Uniaxial: the result fields are None; by definition sigma_y = 0
            # and the interaction factor is 1.
            assert r.residual_transverse_stress_MPa is None
            assert r.biaxial_interaction_factor is None
            assert (row["residual_sigma_y_MPa"], row["pristine_sigma_y_MPa"]) == (0.0, 0.0)
            assert row["biaxial_interaction_factor"] == 1.0
        else:
            assert row["residual_sigma_y_MPa"] == r.residual_transverse_stress_MPa
            assert row["pristine_sigma_y_MPa"] == r.pristine_transverse_stress_MPa
            assert row["biaxial_interaction_factor"] == r.biaxial_interaction_factor


@pytest.mark.parametrize("loading, s", [("compression", -1.0), ("tension", 1.0)])
def test_sigma_y_is_the_signed_ray(loading, s):
    df = sweep_transverse_ratios(_cfg(loading=loading), RATIOS)
    for _, row in df.iterrows():
        beta = row["transverse_stress_ratio"]
        assert row["pristine_sigma_y_MPa"] == pytest.approx(s * beta * row["pristine_MPa"])
        assert row["residual_sigma_y_MPa"] == pytest.approx(s * beta * row["residual_MPa"])


def test_ratio_zero_row_equals_uniaxial_energy_sweep_row():
    cfg = _cfg()
    ratio_row = sweep_transverse_ratios(cfg, [0.0]).iloc[0]
    energy_row = sweep_energies(cfg, [cfg.impact.energy_J]).iloc[0]
    for col in _STANDARD_COLS[:-1]:
        assert ratio_row[col] == energy_row[col], col


def test_inspection_driven_config_works():
    """Unlike ``sweep_energies``, the ratio sweep needs no impact."""
    damage = BvidAnalysis(_cfg()).run().damage
    df = sweep_transverse_ratios(_cfg(impact=None, damage=damage), [0.0, 0.5])
    assert df["error"].isna().all()
    assert df["dpa_mm2"].tolist() == [damage.projected_damage_area_mm2] * 2


def test_csv_round_trip(tmp_path):
    path = tmp_path / "ratios.csv"
    df = sweep_transverse_ratios(_cfg(), [0.0, 0.5], csv_path=path)
    loaded = pd.read_csv(path)
    assert list(loaded.columns) == _RATIO_COLS
    pd.testing.assert_series_equal(
        loaded["biaxial_interaction_factor"], df["biaxial_interaction_factor"]
    )


def test_out_of_range_ratio_is_a_failed_point():
    with pytest.raises(ValueError, match="between -1 and 1"):
        sweep_transverse_ratios(_cfg(), [0.0, 1.5])
    df = sweep_transverse_ratios(_cfg(), [0.0, 1.5], on_error="skip")
    bad = df.iloc[1]
    assert bad["transverse_stress_ratio"] == 1.5
    assert "ValueError" in bad["error"]
    for col in ("knockdown", "residual_sigma_y_MPa", "pristine_sigma_y_MPa"):
        assert math.isnan(bad[col])
    assert math.isnan(bad["biaxial_interaction_factor"])
    assert (bad["n_delaminations"], bad["tier_used"]) == (0, "")


def test_unsupported_tier_is_recorded_not_crashing(monkeypatch):
    """A tier without a biaxial model raises ``NotImplementedError`` at a
    nonzero ratio: re-raised by default, a failed row otherwise. Simulated
    so the test does not depend on which tiers support the ratio yet."""
    real = BvidAnalysis.run

    def _run(self):
        if self.config.transverse_stress_ratio != 0:
            raise NotImplementedError("tier='fake' has no biaxial model")
        return real(self)

    monkeypatch.setattr(BvidAnalysis, "run", _run)
    with pytest.raises(NotImplementedError, match="no biaxial model"):
        sweep_transverse_ratios(_cfg(), [0.0, 0.5])
    with pytest.warns(UserWarning, match="transverse_stress_ratio=0.5 failed"):
        df = sweep_transverse_ratios(_cfg(), [0.0, 0.5], on_error="warn")
    assert pd.isna(df["error"].iloc[0])
    assert df["error"].iloc[1].startswith("NotImplementedError: ")
    assert math.isnan(df["knockdown"].iloc[1])


def test_progress_callback_inputs():
    seen: list = []
    sweep_transverse_ratios(
        _cfg(), [0.0, -0.5], progress_callback=lambda i, n, d: seen.append((i, n, d))
    )
    assert seen == [
        (1, 2, {"transverse_stress_ratio": 0.0}),
        (2, 2, {"transverse_stress_ratio": -0.5}),
    ]


def test_on_error_typo_rejected_before_running():
    with pytest.raises(ValueError, match="on_error"):
        sweep_transverse_ratios(_cfg(), [], on_error="ignore")  # type: ignore[arg-type]


def test_empty_sweep():
    df = sweep_transverse_ratios(_cfg(), [])
    assert df.empty and list(df.columns) == _RATIO_COLS


def test_point_runner_is_the_per_point_seam(monkeypatch):
    """Each point goes through the runner ``_ratio_point_runner`` returns
    (where fe3d's basis-solve cache plugs in)."""
    calls: list = []

    def _runner(base_cfg):
        def run(cfg):
            calls.append(cfg.transverse_stress_ratio)
            return ps._run_one_ratio(cfg)

        return run

    monkeypatch.setattr(ps, "_ratio_point_runner", _runner)
    df = sweep_transverse_ratios(_cfg(), [0.0, 0.5])
    assert calls == [0.0, 0.5]
    assert df["error"].isna().all()


def _count_fe3d(monkeypatch) -> dict:
    """Count fe3d's linear solves and mesh builds."""
    from bvidfe.analysis import fe_tier

    calls = {"solves": 0, "meshes": 0}
    solve, preflight = fe_tier.solve_linear_static, fe_tier._fe3d_preflight

    def _solve(*a, **k):
        calls["solves"] += 1
        return solve(*a, **k)

    def _preflight(*a, **k):
        calls["meshes"] += 1
        return preflight(*a, **k)

    monkeypatch.setattr(fe_tier, "solve_linear_static", _solve)
    monkeypatch.setattr(fe_tier, "_fe3d_preflight", _preflight)
    return calls


def _small_fe3d_cfg(loading: str) -> AnalysisConfig:
    """A fast fe3d config: a 60 x 40 mm panel at 5 mm with two delaminations."""
    from bvidfe.analysis import MeshParams
    from bvidfe.damage.state import DamageState, DelaminationEllipse

    damage = DamageState(
        [
            DelaminationEllipse(2, (30.0, 20.0), 10.0, 6.0, 0.0),
            DelaminationEllipse(5, (30.0, 20.0), 9.0, 5.0, 90.0),
        ]
    )
    return _cfg(
        tier="fe3d",
        loading=loading,
        panel=PanelGeometry(60.0, 40.0),
        impact=None,
        damage=damage,
        mesh=MeshParams(in_plane_size_mm=5.0),
    )


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_fe3d_sweep_solves_each_basis_once(loading, monkeypatch):
    """fe3d (issue #115, PR 3): the basis solves of a mesh do not depend on
    the ratio, so the sweep builds the damaged and the undamaged mesh once
    and solves each twice (x and y) for all five ratios, instead of 18
    solves and 10 meshes for the runs one by one; every row equals the
    run's."""
    calls = _count_fe3d(monkeypatch)
    cfg = _small_fe3d_cfg(loading)
    df = sweep_transverse_ratios(cfg, RATIOS)
    assert calls == {"solves": 4, "meshes": 2}
    calls.update(solves=0, meshes=0)
    rows = [ps._run_one_ratio(replace(cfg, transverse_stress_ratio=b)) for b in RATIOS]
    assert calls == {"solves": 18, "meshes": 10}
    expected = pd.DataFrame(rows)
    for col in expected.columns:
        assert df[col].tolist() == expected[col].tolist(), col
    assert df["error"].isna().all()


def test_fe3d_sweep_of_zero_alone_solves_once_per_mesh(monkeypatch):
    cfg = _small_fe3d_cfg("compression")
    calls = _count_fe3d(monkeypatch)
    df = sweep_transverse_ratios(cfg, [0.0, 0.0])
    assert calls == {"solves": 2, "meshes": 2}
    assert df["knockdown"].iloc[0] == BvidAnalysis(cfg).run().knockdown
