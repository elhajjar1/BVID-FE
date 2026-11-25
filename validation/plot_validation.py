#!/usr/bin/env python3
"""Predicted-versus-measured validation plots for the measured datasets.

Each panel plots a model output against its measured value, one point per
coupon, with the 1:1 line (perfect agreement) and a ±15% band (the
validation target):

* residual strength, closed-form tiers (``empirical``; ``semi_analytical``
  is identical on every measured coupon, since the Soutis bound governs);
* residual strength, ``fe3d`` (``--fe3d``; slow, about 15 min);
* knockdown (residual over pristine strength), where the pristine strength
  was measured on the same coupons (Hasebe et al., Körbelin et al.) or is
  NCAMP's unnotched strength (NCAMP, Lovejoy & Scotti);
* projected damage area against the C-scan area;
* dent depth;
* unnotched (pristine) laminate strength against NCAMP's twelve values.

Run from the repository root::

    python validation/plot_validation.py            # closed-form panels
    python validation/plot_validation.py --fe3d     # also the fe3d panel

Figures and the per-coupon CSV go to ``--out`` (default
``docs/validation_plots``). ``--cache DIR`` reuses per-dataset CSVs written
by an earlier run.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import validate_bvid_public as v  # noqa: E402

from bvidfe.core.laminate import Laminate  # noqa: E402
from bvidfe.core.material import MATERIAL_LIBRARY  # noqa: E402
from bvidfe.failure.laminate_strength import unnotched_strength  # noqa: E402

# The light-surface tokens and the categorical slots live in bvidfe.viz.style,
# shared with the animations. These five slots are validated all-pairs for
# five series (CVD ΔE 6.9), so every series also has its own marker shape.
from bvidfe.viz.style import BAND, CATEGORICAL, GRID, SURFACE, TEXT, TEXT_2  # noqa: E402

STYLE = {
    # dataset: (label, colour, marker)
    "ncamp_as4_8552_cai": ("NCAMP AS4/8552", CATEGORICAL[0], "o"),
    "lovejoy_scotti_im7_8552_cai": ("Lovejoy & Scotti IM7/8552", CATEGORICAL[7], "s"),
    "hasebe_t800s_3900_2b_cai": ("Hasebe et al. T800S/3900-2B (not gated)", CATEGORICAL[2], "^"),
    "gonzalez_as4_8552_cai": ("González AS4/8552 (not gated)", CATEGORICAL[5], "D"),
    "korbelin_m21_t800s_tai": ("Körbelin et al. M21/T800S, TAI (not gated)", CATEGORICAL[3], "v"),
}
BAND_PCT = 15.0

# NCAMP RTD normalized mean unnotched strengths, ksi (tension, compression);
# the same values as tests/validation/test_unnotched_strength_ncamp.py.
_KSI = 6.894757
_NCAMP_UNNOTCHED_KSI = {
    ("AS4/8552", "25/50/25"): (88.61, 81.32),
    ("AS4/8552", "10/80/10"): (63.62, 62.42),
    ("AS4/8552", "50/40/10"): (152.32, 131.05),
    ("IM7/8552", "25/50/25"): (104.69, 87.05),
    ("IM7/8552", "10/80/10"): (67.01, 66.44),
    ("IM7/8552", "50/40/10"): (175.63, 120.84),
}
_CPT_MM = {"AS4/8552": 0.0074 * 25.4, "IM7/8552": 0.0072 * 25.4}


def _sym(half: list) -> list:
    return half + half[::-1]


_NCAMP_LAYUPS = {
    "25/50/25": _sym([45, 0, -45, 90] * 2),
    "10/80/10": _sym([45, -45, 0, 45, -45, 90, 45, -45, 45, -45]),
    "50/40/10": _sym([0, 45, 0, 90, 0, -45, 0, 45, 0, -45]),
}


_NUMERIC = (
    "predicted_MPa",
    "measured_MPa",
    "predicted_knockdown",
    "measured_knockdown",
    "predicted_dpa_mm2",
    "measured_dpa_mm2",
    "predicted_dent_mm",
    "measured_dent_mm",
    "predicted_pristine_MPa",
    "measured_pristine_MPa",
)


def run_tier(tier: str, cache: Path | None, names: list[str] | None = None) -> pd.DataFrame:
    """Per-coupon results of ``tier`` on the measured datasets."""
    frames = []
    for path in v._select_datasets(names, kind="measured"):
        name, cases, _ = v.load_dataset(path)
        cached = cache / f"{tier}_{name}.csv" if cache else None
        if cached is not None and cached.exists():
            df = pd.read_csv(cached)
        else:
            df = v.run_dataset(cases, tier)
            df["measured_dpa_mm2"] = [c.measured_dpa_mm2 for c in cases]
            df["dataset"] = name
            if cached is not None:
                cached.parent.mkdir(parents=True, exist_ok=True)
                df.to_csv(cached, index=False)
        df["tier"] = tier
        # A fixed float dtype, so a column that is all-NA in one dataset
        # (NCAMP has no C-scan areas) concatenates without a dtype change.
        for col in _NUMERIC:
            if col in df:
                df[col] = pd.to_numeric(df[col], errors="coerce").astype(float)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def ncamp_unnotched() -> pd.DataFrame:
    rows = []
    for (mat, layup), pair in _NCAMP_UNNOTCHED_KSI.items():
        lam = Laminate(MATERIAL_LIBRARY[mat], _NCAMP_LAYUPS[layup], _CPT_MM[mat])
        for loading, ksi in zip(("tension", "compression"), pair):
            rows.append(
                {
                    "material": mat,
                    "layup": layup,
                    "loading": loading,
                    "measured_MPa": ksi * _KSI,
                    "predicted_MPa": unnotched_strength(lam, loading),
                }
            )
    return pd.DataFrame(rows)


def _mae_pct(pred, meas) -> float:
    pred, meas = np.asarray(pred, float), np.asarray(meas, float)
    return float(np.mean(np.abs(pred / meas - 1.0)) * 100.0)


def _limits(df: pd.DataFrame, cols: tuple, log: bool, floor: float | None = None) -> tuple:
    """Axis range covering every plotted value of ``cols`` (both axes share it)."""
    vals = pd.concat([pd.to_numeric(df[c], errors="coerce") for c in cols]).dropna()
    vals = vals[vals > 0]
    if log:
        return float(vals.min()) / 1.3, float(vals.max()) * 1.3
    return (0.0 if floor is None else floor), float(vals.max()) * 1.06


def _frame(ax, lo: float, hi: float, log: bool) -> None:
    """Equal axes from lo to hi, the ±15% band and the 1:1 line."""
    x = np.geomspace(lo, hi, 200) if log else np.linspace(lo, hi, 200)
    ax.fill_between(
        x, x * (1 - BAND_PCT / 100), x * (1 + BAND_PCT / 100), color=BAND, lw=0, zorder=0
    )
    ax.plot([lo, hi], [lo, hi], color=TEXT_2, lw=1.2, zorder=1)
    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")
        plain = FuncFormatter(lambda val, _: f"{val:g}")
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
            axis.set_major_formatter(plain)
            axis.set_minor_formatter(NullFormatter())
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_2)
    ax.tick_params(colors=TEXT_2, labelsize=8)


def _scatter(ax, df: pd.DataFrame, xcol: str, ycol: str, label_extra: str = "") -> list:
    """One series per dataset; returns the legend handles."""
    handles = []
    for name, (label, colour, marker) in STYLE.items():
        sub = df[df["dataset"] == name][[xcol, ycol]].dropna()
        sub = sub[(sub[xcol] > 0) & (sub[ycol] > 0)]
        if sub.empty:
            continue
        dense = len(sub) > 50
        h = ax.scatter(
            sub[xcol],
            sub[ycol],
            s=16 if dense else 34,
            marker=marker,
            color=colour,
            alpha=0.55 if dense else 0.95,
            edgecolors=SURFACE,
            linewidths=0.6,
            zorder=3 if dense else 4,
            label=f"{label}: n={len(sub)}, MAE {_mae_pct(sub[ycol], sub[xcol]):.1f}%",
        )
        handles.append(h)
    return handles


def _finish(ax, title: str, xlabel: str, ylabel: str, handles: list, loc="upper left") -> None:
    ax.set_title(title, fontsize=10, color=TEXT, loc="left")
    ax.set_xlabel(xlabel, fontsize=9, color=TEXT_2)
    ax.set_ylabel(ylabel, fontsize=9, color=TEXT_2)
    if len(handles) >= 2:
        leg = ax.legend(
            handles=handles, fontsize=6.5, loc=loc, frameon=True, framealpha=0.92, borderpad=0.4
        )
        leg.get_frame().set_edgecolor(GRID)
        for t in leg.get_texts():
            t.set_color(TEXT)


def panel_strength(ax, df: pd.DataFrame, tier_note: str) -> None:
    _frame(ax, *_limits(df, ("measured_MPa", "predicted_MPa"), log=True), log=True)
    h = _scatter(ax, df, "measured_MPa", "predicted_MPa")
    _finish(
        ax,
        f"Residual strength, {tier_note}",
        "Measured residual strength (MPa)",
        "Predicted residual strength (MPa)",
        h,
    )


def panel_knockdown(ax, df: pd.DataFrame) -> None:
    _frame(ax, *_limits(df, ("measured_knockdown", "predicted_knockdown"), log=False), log=False)
    h = _scatter(ax, df, "measured_knockdown", "predicted_knockdown")
    _finish(
        ax,
        "Knockdown (residual / pristine strength)",
        "Measured knockdown",
        "Predicted knockdown",
        h,
    )
    ax.text(
        0.98,
        0.03,
        "pristine: same-coupon (Hasebe, Körbelin)\nor NCAMP unnotched (NCAMP, L&S)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.5,
        color=TEXT_2,
    )


def panel_dpa(ax, df: pd.DataFrame) -> None:
    _frame(ax, *_limits(df, ("measured_dpa_mm2", "predicted_dpa_mm2"), log=True), log=True)
    h = _scatter(ax, df, "measured_dpa_mm2", "predicted_dpa_mm2")
    _finish(
        ax,
        "Projected damage area",
        "Measured C-scan area (mm²)",
        "Predicted damage area (mm²)",
        h,
    )


def panel_dent(ax, df: pd.DataFrame) -> None:
    d = df[pd.to_numeric(df["measured_dent_mm"], errors="coerce") > 0]
    _frame(ax, *_limits(d, ("measured_dent_mm", "predicted_dent_mm"), log=True), log=True)
    h = _scatter(ax, df, "measured_dent_mm", "predicted_dent_mm")
    _finish(ax, "Dent depth", "Measured dent depth (mm)", "Predicted dent depth (mm)", h)


def panel_unnotched(ax, un: pd.DataFrame) -> None:
    _frame(ax, *_limits(un, ("measured_MPa", "predicted_MPa"), log=False), log=False)
    colour = STYLE["ncamp_as4_8552_cai"][1]
    handles = []
    for loading, marker, face in (("tension", "^", colour), ("compression", "v", SURFACE)):
        sub = un[un["loading"] == loading]
        handles.append(
            ax.scatter(
                sub["measured_MPa"],
                sub["predicted_MPa"],
                s=40,
                marker=marker,
                facecolors=face,
                edgecolors=colour,
                linewidths=1.4,
                zorder=4,
                label=(
                    f"NCAMP unnotched {loading}: n={len(sub)}, "
                    f"MAE {_mae_pct(sub['predicted_MPa'], sub['measured_MPa']):.1f}%"
                ),
            )
        )
    # Label the points more than 10% off; the four 10/80/10 points share one label.
    soft = un[un["layup"] == "10/80/10"]
    ax.annotate(
        "10/80/10, both materials\n(10-17% low)",
        (soft["measured_MPa"].max(), soft["predicted_MPa"].min()),
        xytext=(10, -14),
        textcoords="offset points",
        fontsize=6.5,
        color=TEXT_2,
    )
    for _, r in un[un["layup"] != "10/80/10"].iterrows():
        ratio = r["predicted_MPa"] / r["measured_MPa"]
        if abs(ratio - 1) > 0.10:
            ax.annotate(
                f"{r['material'].split('/')[0]} {r['layup']} {r['loading'][0].upper()} "
                f"({100 * (ratio - 1):+.0f}%)",
                (r["measured_MPa"], r["predicted_MPa"]),
                xytext=(8, -3),
                textcoords="offset points",
                fontsize=6.5,
                color=TEXT_2,
            )
    _finish(
        ax,
        "Unnotched (pristine) laminate strength",
        "NCAMP measured strength (MPa)",
        "Predicted strength (MPa)",
        handles,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=Path("docs/validation_plots"))
    ap.add_argument("--cache", type=Path, default=None)
    ap.add_argument("--fe3d", action="store_true", help="also run fe3d (about 15 min)")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE})

    emp = run_tier("empirical", args.cache)
    frames = [emp]
    fe = None
    if args.fe3d:
        fe = run_tier(
            "fe3d",
            args.cache,
            [
                "ncamp_as4_8552_cai",
                "lovejoy_scotti_im7_8552_cai",
                "gonzalez_as4_8552_cai",
                "korbelin_m21_t800s_tai",
            ],
        )
        frames.append(fe)
    un = ncamp_unnotched()
    pd.concat(frames, ignore_index=True).to_csv(args.out / "validation_points.csv", index=False)
    un.to_csv(args.out / "ncamp_unnotched_points.csv", index=False)

    panels = [
        ("strength", lambda ax: panel_strength(ax, emp, "empirical and semi_analytical")),
        ("knockdown", lambda ax: panel_knockdown(ax, emp)),
        ("damage_area", lambda ax: panel_dpa(ax, emp)),
        ("dent_depth", lambda ax: panel_dent(ax, emp)),
        ("unnotched_strength", lambda ax: panel_unnotched(ax, un)),
    ]
    if fe is not None:
        panels.insert(1, ("strength_fe3d", lambda ax: panel_strength(ax, fe, "fe3d (5 mm mesh)")))

    for key, draw in panels:
        fig, ax = plt.subplots(figsize=(5.2, 5.2), dpi=150)
        draw(ax)
        fig.tight_layout()
        fig.savefig(args.out / f"{key}.png")
        plt.close(fig)

    ncols = 3
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 5.0 * nrows), dpi=150)
    for ax in axes.flat[len(panels) :]:
        ax.set_visible(False)
    for ax, (_, draw) in zip(axes.flat, panels):
        draw(ax)
    fig.suptitle(
        "BVID-FE: predicted against measured (1:1 line; shaded ±15%)",
        fontsize=12,
        color=TEXT,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(args.out / "validation_overview.png")
    plt.close(fig)
    print(f"wrote {len(panels)} panels and validation_overview.png to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
