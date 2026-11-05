"""``plots_2d.biaxial_envelope_plot`` (issue #115, PR 4)."""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from bvidfe.analysis import AnalysisConfig  # noqa: E402
from bvidfe.core.geometry import PanelGeometry  # noqa: E402
from bvidfe.impact.mapping import ImpactEvent  # noqa: E402
from bvidfe.sweep.parametric_sweep import sweep_transverse_ratios  # noqa: E402
from bvidfe.viz.plots_2d import biaxial_envelope_plot  # noqa: E402
from bvidfe.viz.style import COLORS  # noqa: E402

RATIOS = [-1.0, -0.5, 0.0, 0.5, 1.0]


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _sweep(loading: str) -> pd.DataFrame:
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading=loading,
        impact=ImpactEvent(15.0),
    )
    return sweep_transverse_ratios(cfg, RATIOS)


def _line(ax, label):
    (line,) = [ln for ln in ax.get_lines() if ln.get_label() == label]
    return line


def test_plots_both_halves_and_saves_png(tmp_path):
    df_c, df_t = _sweep("compression"), _sweep("tension")
    fig = biaxial_envelope_plot(df_c, df_t)
    (ax,) = fig.axes
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels == ["pristine", "damaged (empirical)"]
    assert "not validated" in ax.get_title()
    out = tmp_path / "envelope.png"
    fig.savefig(out, dpi=100)
    assert out.exists() and out.stat().st_size > 5000

    # The compression half: sigma_x = -|sigma_x|, sigma_y the signed column.
    pristine = _line(ax, "pristine")
    np.testing.assert_allclose(pristine.get_xdata(), -df_c["pristine_MPa"])
    np.testing.assert_allclose(pristine.get_ydata(), df_c["pristine_sigma_y_MPa"])
    damaged = _line(ax, "damaged (empirical)")
    np.testing.assert_allclose(damaged.get_xdata(), -df_c["residual_MPa"])
    np.testing.assert_allclose(damaged.get_ydata(), df_c["residual_sigma_y_MPa"])
    assert damaged.get_color() == COLORS["empirical"]
    # The tension half is on the sigma_x > 0 side.
    tension_lines = [ln for ln in ax.get_lines() if (np.asarray(ln.get_xdata()) > 0).all()]
    assert len(tension_lines) == 2  # pristine and damaged


def test_one_half_only_and_fig_reuse():
    df_c = _sweep("compression")
    fig = plt.figure()
    out = biaxial_envelope_plot(df_compression=df_c, fig=fig)
    assert out is fig
    (ax,) = fig.axes
    data_lines = [ln for ln in ax.get_lines() if not ln.get_label().startswith("_")]
    assert len(data_lines) == 2
    assert all((np.asarray(ln.get_xdata()) <= 0).all() for ln in data_lines)
    out = biaxial_envelope_plot(df_tension=_sweep("tension"), fig=fig, title="custom")
    (ax,) = fig.axes
    assert ax.get_title() == "custom"


def test_failed_rows_are_dropped_and_tiers_are_separate():
    df = _sweep("compression")
    other = df.copy()
    other["tier_used"] = "semi_analytical"
    other["residual_MPa"] *= 0.9
    other["residual_sigma_y_MPa"] *= 0.9
    other.loc[other["transverse_stress_ratio"] == 1.0, ["residual_MPa", "residual_sigma_y_MPa"]] = (
        math.nan
    )
    fig = biaxial_envelope_plot(pd.concat([df, other], ignore_index=True))
    (ax,) = fig.axes
    semi = _line(ax, "damaged (semi_analytical)")
    assert len(semi.get_xdata()) == len(RATIOS) - 1
    assert semi.get_color() == COLORS["semi_analytical"]
    assert len(_line(ax, "damaged (empirical)").get_xdata()) == len(RATIOS)


def test_all_failed_says_so():
    df = _sweep("compression")
    df[["residual_MPa", "residual_sigma_y_MPa", "pristine_MPa", "pristine_sigma_y_MPa"]] = math.nan
    df["tier_used"] = ""
    fig = biaxial_envelope_plot(df)
    (ax,) = fig.axes
    assert any("no successful sweep points" in t.get_text() for t in ax.texts)
    assert ax.get_legend() is None


def test_needs_a_frame():
    with pytest.raises(ValueError, match="df_compression, df_tension or both"):
        biaxial_envelope_plot()
