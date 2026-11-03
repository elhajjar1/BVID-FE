"""Streamlit web app for BVID-FE residual-strength analysis.

Wraps ``bvidfe.analysis.BvidAnalysis`` in a browser UI: sidebar inputs
build an ``AnalysisConfig``, the main area shows results across tabs
(Summary, Damage Map, Knockdown Curve, 3D Damage, Buckling, Damage
Severity, Sweep, Export NCR). The app does not expose
``AnalysisConfig.impact_panel``: the impact is applied on the panel that
is tested.

Run locally:
    streamlit run app.py

Deploy: see DEPLOYMENT_STREAMLIT.md.
"""

from __future__ import annotations

import io
import json
import math
import sys
import time
import warnings
from dataclasses import replace
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend; required on Streamlit Cloud

import matplotlib.pyplot as plt  # noqa: E402  (matplotlib.use must run first)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

# Make the src-layout package importable on Streamlit Cloud, which clones
# the repo but does not necessarily editable-install it before the first
# script run.
_SRC = Path(__file__).resolve().parent / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams  # noqa: E402
from bvidfe.analysis.fe_mesh import (  # noqa: E402
    build_fe_mesh,
    element_grid_shape,
    estimate_fe_mesh_size,
)
from bvidfe.analysis.fe_tier import FE3D_MAX_DOF, enable_stderr_logging  # noqa: E402
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry  # noqa: E402
from bvidfe.core.material import MATERIAL_LIBRARY  # noqa: E402
from bvidfe.damage.io import CScanSchemaError, damage_state_from_dict  # noqa: E402
from bvidfe.damage.state import DamageState  # noqa: E402
from bvidfe.impact.mapping import (  # noqa: E402
    DIAMETER_DPA_EXPONENT,
    DIAMETER_REF_MM,
    SHAPE_DPA_FACTOR,
    ImpactEvent,
    impact_to_damage,
)
from bvidfe.impact.olsson import onset_energy  # noqa: E402
from bvidfe.report.ncr import (  # noqa: E402
    BVID_DENT_THRESHOLD_MM,
    build_ncr,
    recommend_disposition,
)
from bvidfe.core.laminate import Laminate  # noqa: E402
from bvidfe.sweep.parametric_sweep import sweep_energies  # noqa: E402
from bvidfe.viz.plots_2d import plot_damage_map, plot_knockdown_curve  # noqa: E402
from bvidfe.viz.plotly_3d import mesh_damage_figure  # noqa: E402

# fe3d stage logs to the terminal running `streamlit run` (idempotent across reruns).
enable_stderr_logging()

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="BVID-FE",
    page_icon=":airplane:",
    layout="wide",
)

st.title("BVID-FE — Barely Visible Impact Damage analysis")
st.caption(
    "Predicts the compression- and tension-after-impact residual strength "
    "of composite laminates after low-velocity impact. Configure parameters "
    "in the sidebar and click **Run analysis**."
)


# ---------------------------------------------------------------------------
# Constants / helpers
# ---------------------------------------------------------------------------

MATERIAL_NAMES: list[str] = sorted(MATERIAL_LIBRARY.keys())
DEFAULT_MATERIAL = "IM7/8552" if "IM7/8552" in MATERIAL_NAMES else MATERIAL_NAMES[0]
TIERS: list[str] = ["empirical", "semi_analytical", "fe3d"]
BOUNDARIES: list[str] = ["simply_supported", "clamped", "free"]
IMPACTOR_SHAPES: list[str] = ["hemispherical", "flat", "conical"]


MAX_PLIES = 400
MAX_SWEEP_POINTS = 25
DPA_SATURATION_PCT = 79.0

# MRB / NCR disposition support -------------------------------------------
# The BVID dent-detectability threshold (``BVID_DENT_THRESHOLD_MM``, the
# default of a UI input each program must override) and the disposition
# bands are in ``bvidfe.report.ncr``.

NC_TYPES: list[str] = [
    "Impact damage (BVID/VID)",
    "Delamination / disbond",
    "Porosity / voids",
    "Foreign object / inclusion",
    "Resin-rich / resin-starved",
    "Wrinkle / fiber waviness",
    "Out-of-tolerance geometry",
    "Other",
]

FOUND_DURING: list[str] = [
    "In-process inspection",
    "Final inspection",
    "Assembly / installation",
    "In-service / field inspection",
    "Receiving inspection",
    "Other",
]

_DF_FORMATS = {
    "energy_J": "{:g}",
    "knockdown": "{:.3f}",
    "residual_MPa": "{:.1f}",
    "pristine_MPa": "{:.1f}",
    "dpa_mm2": "{:.1f}",
    "dent_mm": "{:.3f}",
}


# ---------------------------------------------------------------------------
# Presets and shareable URL state
# ---------------------------------------------------------------------------
#
# Presets are canonical sidebar configurations. Picking a preset writes its
# values into ``st.session_state`` BEFORE the widgets render (so the widgets
# pick them up via their ``key=`` arg). A "preset-applied" guard
# (``_last_preset_key``) prevents a re-applied preset from clobbering user
# edits on every rerun — only a *change* in preset selection triggers a
# rewrite of the controlled keys.
#
# URL state uses short keys (``mat``, ``lay``, ``lx``, ``ly``, ``E``, ``d``,
# ``t``, ``ld`` and the expert-mode ``tp``, ``bc``, ``sh``, ``m``, ``x0``,
# ``y0``) so a shared link stays well under the 2 kB safe ceiling for
# common browsers. On app load, if the URL carries known keys, we hydrate
# session state from them BEFORE widgets render (same mechanism as presets).

# Expert-mode widget values (session-state key -> default). These widgets
# only render in Expert mode, but their values live in session state either
# way, so toggling Expert mode off and on keeps them and a share link or
# preset carries them. ``loc_centre`` True means "impact the panel centre";
# ``loc_x_mm`` / ``loc_y_mm`` (from the panel corner) are then ignored.
_EXPERT_DEFAULTS: dict = {
    "ply_thickness_mm": 0.152,
    "boundary": "simply_supported",
    "impactor_shape": "hemispherical",
    "mass_kg": 5.5,
    "loc_centre": True,
}

PRESET_CUSTOM = "Custom"
PRESETS: dict[str, dict] = {
    PRESET_CUSTOM: {},  # sentinel — no overrides
    "Empirical quick (IM7/8552 quasi-iso, 20 J)": {
        "material_choice": "IM7/8552",
        "layup_text": "0, 45, -45, 90, 90, -45, 45, 0",
        "Lx_mm": 150.0,
        "Ly_mm": 100.0,
        "energy_J": 20.0,
        "diameter_mm": 16.0,
        "tier": "empirical",
        "loading": "compression",
        "input_mode": "Impact event",
        **_EXPERT_DEFAULTS,
    },
    "Semi-analytical buckling (IM7/8552 quasi-iso, 10 J)": {
        "material_choice": "IM7/8552",
        "layup_text": "0, 45, -45, 90, 90, -45, 45, 0",
        "Lx_mm": 150.0,
        "Ly_mm": 100.0,
        "energy_J": 10.0,
        "diameter_mm": 16.0,
        "tier": "semi_analytical",
        "loading": "compression",
        "input_mode": "Impact event",
        **_EXPERT_DEFAULTS,
    },
    "FE3D detailed (T800/epoxy, 15 J, coarse mesh)": {
        "material_choice": "T800/epoxy",
        "layup_text": "0, 45, -45, 90, 90, -45, 45, 0",
        "Lx_mm": 150.0,
        "Ly_mm": 100.0,
        "energy_J": 15.0,
        "diameter_mm": 16.0,
        "tier": "fe3d",
        "loading": "compression",
        "input_mode": "Impact event",
        **_EXPERT_DEFAULTS,
    },
}

PRESET_NAMES: list[str] = list(PRESETS.keys())

# Widget keys we control via preset / URL hydration. Each maps a sidebar
# widget's ``key`` argument to the type its widget expects in session_state.
_URL_TO_STATE: dict[str, tuple[str, type]] = {
    "mat": ("material_choice", str),
    "lay": ("layup_text", str),
    "lx": ("Lx_mm", float),
    "ly": ("Ly_mm", float),
    "E": ("energy_J", float),
    "d": ("diameter_mm", float),
    "t": ("tier", str),
    "ld": ("loading", str),
    "tp": ("ply_thickness_mm", float),
    "bc": ("boundary", str),
    "sh": ("impactor_shape", str),
    "m": ("mass_kg", float),
    "x0": ("loc_x_mm", float),
    "y0": ("loc_y_mm", float),
}

# Bounds of each numeric widget, so a URL value the widget would reject (and
# Streamlit would silently replace by its minimum) is dropped and reported
# instead. The location bounds are the panel's and are checked separately.
_STATE_BOUNDS: dict[str, tuple[float, float]] = {
    "Lx_mm": (10.0, 2000.0),
    "Ly_mm": (10.0, 2000.0),
    "energy_J": (0.1, 200.0),
    "diameter_mm": (1.0, 100.0),
    "ply_thickness_mm": (0.01, 2.0),
    "mass_kg": (0.01, 100.0),
}

# Allowed values of each choice widget.
_STATE_CHOICES: dict[str, list[str]] = {
    "material_choice": MATERIAL_NAMES,
    "tier": TIERS,
    "loading": ["compression", "tension"],
    "input_mode": ["Impact event", "C-scan inspection"],
    "boundary": BOUNDARIES,
    "impactor_shape": IMPACTOR_SHAPES,
}

_LOCATION_KEYS = ("loc_x_mm", "loc_y_mm")


def _coerce_layup_url(text: str) -> str:
    """Sanity-check a URL-supplied layup string before pushing into a widget.

    Falls back to the default if the string parses to nothing valid; the
    widget itself will surface a friendly error on bad user-edited text but
    we don't want a malformed URL to crash before widgets even render.
    """
    try:
        angles = _parse_layup(text)
    except ValueError:
        return text  # let the widget echo the original; user can fix it
    if not angles:
        return text
    return ", ".join(f"{a:g}" for a in angles)


def _apply_state_overrides(overrides: dict) -> list[str]:
    """Write controlled overrides into ``st.session_state`` with type coercion.

    Unknown keys are skipped. A value that does not cast, is not finite, is
    outside its widget's bounds (``_STATE_BOUNDS``; the impact location must
    lie on the panel) or is not one of a choice widget's options is dropped,
    so its widget keeps its current value, and is returned as a
    ``"key=value (reason)"`` string for the caller to report: a stale shared
    link still loads, but says what it could not reproduce.
    """
    state_to_url = {state_key: url_key for url_key, (state_key, _) in _URL_TO_STATE.items()}
    type_map = {state_key: cast for (_, (state_key, cast)) in _URL_TO_STATE.items()}
    # widgets controlled only by presets (no URL key) — keep validation here.
    extra_types = {"input_mode": str, "loc_centre": bool}
    rejected: list[str] = []

    def _reject(key: str, value, reason: str) -> None:
        rejected.append(f"{state_to_url.get(key, key)}={value} ({reason})")

    # The location is checked against the panel, so after Lx / Ly are applied.
    ordered = [k for k in overrides if k not in _LOCATION_KEYS]
    ordered += [k for k in _LOCATION_KEYS if k in overrides]
    for key in ordered:
        value = overrides[key]
        cast = type_map.get(key, extra_types.get(key))
        if cast is None:
            continue
        if cast is bool:
            if not isinstance(value, bool):
                _reject(key, value, "not true/false")
                continue
            coerced = value
        else:
            try:
                coerced = cast(value)
            except (TypeError, ValueError):
                _reject(key, value, "not a number" if cast is float else "invalid")
                continue
        if cast is float:
            if not math.isfinite(coerced):
                _reject(key, value, "not finite")
                continue
            if key in _LOCATION_KEYS:
                span = "Lx_mm" if key == "loc_x_mm" else "Ly_mm"
                lo, hi = 0.0, float(st.session_state.get(span, 150.0 if span == "Lx_mm" else 100.0))
            else:
                lo, hi = _STATE_BOUNDS[key]
            if not lo <= coerced <= hi:
                _reject(key, value, f"outside {lo:g}-{hi:g}")
                continue
        if key == "layup_text":
            coerced = _coerce_layup_url(str(coerced))
        if key in _STATE_CHOICES and coerced not in _STATE_CHOICES[key]:
            _reject(key, value, "unknown option")
            continue
        st.session_state[key] = coerced
        if key in _LOCATION_KEYS:
            st.session_state["loc_centre"] = False
    return rejected


def _hydrate_from_url() -> None:
    """Populate session state from ``st.query_params`` on first script run.

    Only runs once per session (guarded by ``_url_hydrated``) so user edits
    after load are not overwritten by the original URL on every rerun.
    Values the link carried but that were dropped are kept in
    ``_url_rejected`` for the page to report.
    """
    if st.session_state.get("_url_hydrated"):
        return
    st.session_state["_url_hydrated"] = True
    qp = st.query_params
    overrides: dict = {}
    for url_key, (state_key, _cast) in _URL_TO_STATE.items():
        if url_key in qp:
            overrides[state_key] = qp[url_key]
    if overrides:
        rejected = _apply_state_overrides(overrides)
        if rejected:
            st.session_state["_url_rejected"] = rejected


def _hydrate_from_preset() -> None:
    """Apply preset overrides when the user picks a non-Custom preset.

    Only fires when the selection *changes* (tracked in ``_last_preset_key``);
    re-runs with the same selection are no-ops, so user edits after applying
    a preset survive subsequent script reruns.
    """
    chosen = st.session_state.get("preset_choice", PRESET_CUSTOM)
    last = st.session_state.get("_last_preset_key")
    if chosen == last:
        return
    st.session_state["_last_preset_key"] = chosen
    if chosen == PRESET_CUSTOM:
        return
    _apply_state_overrides(PRESETS[chosen])


def _sync_query_params() -> None:
    """Write the live sidebar config back to ``st.query_params``.

    The shareable subset (material, layup, panel and boundary, ply
    thickness, impact energy, impactor diameter / shape / mass, tier,
    loading, and the impact location when it is not the panel centre) is
    encoded — short keys keep URLs under the ~2 kB safe ceiling. Values are
    stringified; Streamlit handles URL escaping for us.
    """
    ss = st.session_state
    try:
        new_qp = {
            "mat": str(ss.get("material_choice", DEFAULT_MATERIAL)),
            "lay": str(ss.get("layup_text", "")),
            "lx": f"{float(ss.get('Lx_mm', 150.0)):g}",
            "ly": f"{float(ss.get('Ly_mm', 100.0)):g}",
            "E": f"{float(ss.get('energy_J', 30.0)):g}",
            "d": f"{float(ss.get('diameter_mm', 16.0)):g}",
            "t": str(ss.get("tier", "empirical")),
            "ld": str(ss.get("loading", "compression")),
            "tp": f"{float(ss.get('ply_thickness_mm', _EXPERT_DEFAULTS['ply_thickness_mm'])):g}",
            "bc": str(ss.get("boundary", _EXPERT_DEFAULTS["boundary"])),
            "sh": str(ss.get("impactor_shape", _EXPERT_DEFAULTS["impactor_shape"])),
            "m": f"{float(ss.get('mass_kg', _EXPERT_DEFAULTS['mass_kg'])):g}",
        }
        if not ss.get("loc_centre", True):
            new_qp["x0"] = f"{float(ss.get('loc_x_mm', 0.0)):g}"
            new_qp["y0"] = f"{float(ss.get('loc_y_mm', 0.0)):g}"
    except (TypeError, ValueError):
        return
    # Only write if the dict actually changes — avoids needless reruns.
    # (The whole dict, so a dropped ``x0`` / ``y0`` is noticed too.)
    current = st.query_params.to_dict()
    if current != new_qp:
        try:
            st.query_params.from_dict(new_qp)
        except Exception:
            # Streamlit's query-param API is occasionally fussy on stale
            # sessions; never let URL syncing crash the app.
            pass


def _show_df(df: pd.DataFrame) -> None:
    """Render a results DataFrame with sane float formatting."""
    fmt = {c: f for c, f in _DF_FORMATS.items() if c in df.columns}
    st.dataframe(df.style.format(fmt, na_rep="—"), width="stretch")


def _parse_layup(text: str) -> list[float]:
    """Parse a comma- or whitespace-separated layup string into floats.

    Rejects non-finite values (nan/inf), out-of-range angles, and absurdly
    long stacks so a malformed paste can't produce NaN ABD matrices or
    exhaust memory.
    """
    cleaned = text.replace(";", ",").replace("\n", ",")
    tokens = [x.strip() for x in cleaned.split(",") if x.strip()]
    if len(tokens) > MAX_PLIES:
        raise ValueError(f"layup has {len(tokens)} plies; the cap is {MAX_PLIES}")
    angles = [float(t) for t in tokens]
    for a in angles:
        if not math.isfinite(a):
            raise ValueError("ply angles must be finite numbers")
        if not -90.0 <= a <= 90.0:
            raise ValueError("ply angles must be within [-90, 90] degrees")
    return angles


def _knockdown_curve_energies(base_E: float, n: int = 7) -> tuple[float, ...]:
    """Ascending energy grid for the Knockdown Curve tab around ``base_E``.

    From 0.2·E to 2·E, but at least 1 J wide, so the grid never descends
    (the old ``max(1, 0.2·E)`` floor did for E < 1.25 J).
    """
    lo = 0.2 * base_E
    hi = max(2.0 * base_E, lo + 1.0)
    return tuple(np.linspace(lo, hi, n).tolist())


# Shared by the Sweep and Knockdown Curve tabs' tier selectors.
_SWEEP_TIER_HELP = (
    "Use **empirical** or **semi_analytical** for fast scans. "
    "**fe3d** takes about 15-30 s per point at the 5 mm default "
    "mesh (longer for thick laminates or finer meshes); keep "
    "fe3d sweeps short on the 1 GB cloud tier."
)


@st.cache_data(show_spinner=False, max_entries=8)
def _run_analysis_cached(config_dict: dict) -> dict:
    """Run BvidAnalysis from a JSON-safe config dict and return a dict.

    Caching key is the config dict itself, so wiggling unrelated UI
    widgets won't re-trigger a 30-second fe3d solve. ``st.cache_data``
    returns a fresh copy on each access, so the contained objects are
    safe for the result tabs to consume. ``max_entries`` bounds memory
    on the 1 GB Streamlit Cloud tier. The display mesh is built here
    (once per unique config) instead of on every rerun in two tabs.
    """
    cfg = _config_from_dict(config_dict)
    result = BvidAnalysis(cfg).run()
    cfg_mesh = cfg if cfg.mesh is not None else replace(cfg, mesh=MeshParams())
    try:
        mesh = build_fe_mesh(cfg_mesh, result.damage)
    except Exception:
        # A display-mesh failure must not sink the numeric result; the
        # 3D/Severity tabs surface the absence with a friendly message.
        mesh = None
    return {
        "result_dict": result.to_dict(),
        "damage": result.damage,
        "config": cfg,
        "config_dict": config_dict,
        "mesh_config": cfg_mesh,
        "mesh": mesh,
        "buckling_eigenvalues": result.buckling_eigenvalues,
        "notes": list(result.notes),
        "summary_text": result.summary(),
    }


@st.cache_data(show_spinner=False, max_entries=16)
def _sweep_cached(config_dict: dict, tier: str, energies: tuple[float, ...]) -> pd.DataFrame:
    """Energy sweep for a given tier, cached and resilient to per-point failures.

    ``on_error="warn"`` makes the library return NaN rows for points that
    fail instead of aborting the whole sweep.
    """
    cfg = replace(_config_from_dict(config_dict), tier=tier)
    return sweep_energies(cfg, list(energies), on_error="warn", progress_callback=None)


@st.cache_data(show_spinner=False, max_entries=32)
def _onset_preview_cached(
    material: str,
    layup: tuple[float, ...],
    ply_thickness_mm: float,
    Lx_mm: float,
    Ly_mm: float,
    boundary: str,
    energy_J: float,
    diameter_mm: float,
    shape: str,
    mass_kg: float,
    loc_xy: tuple[float, float] | None,
) -> tuple[str | None, str | None]:
    """Cached Olsson onset/DPA preview so it isn't recomputed every keystroke.

    Returns ``(preview, None)``, or ``(None, error text)`` when the preview
    cannot be computed, so the page can say why instead of dropping it.
    """
    try:
        mat = MATERIAL_LIBRARY[material]
        lam = Laminate(mat, list(layup), ply_thickness_mm)
        pan = PanelGeometry(Lx_mm=Lx_mm, Ly_mm=Ly_mm, boundary=boundary)
        impactor = ImpactorGeometry(diameter_mm=diameter_mm, shape=shape)
        ev = ImpactEvent(
            energy_J=energy_J,
            impactor=impactor,
            mass_kg=mass_kg,
            location_xy_mm=loc_xy,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            E_onset = onset_energy(lam, pan, impactor)
            dmg = impact_to_damage(ev, lam, pan)
        A_panel = Lx_mm * Ly_mm
        pct = 100.0 * dmg.projected_damage_area_mm2 / A_panel if A_panel else 0.0
        flag = " ⚠ SATURATED" if pct >= DPA_SATURATION_PCT else ""
        return (
            f"E_onset ≈ **{E_onset:.2f} J** · "
            f"predicted DPA ≈ **{dmg.projected_damage_area_mm2:.0f} mm²** "
            f"({pct:.1f}% of panel){flag}"
        ), None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def _ply_thickness_from_json(t):
    """A uniform thickness (number) or per-ply thicknesses (list)."""
    if isinstance(t, (list, tuple)):
        return [float(x) for x in t]
    return float(t)


def _format_ply_thickness(t) -> str:
    if isinstance(t, (list, tuple)):
        return ", ".join(f"{x:g}" for x in t) + " mm (per ply)"
    return f"{t:g} mm"


def _config_from_dict(d: dict) -> AnalysisConfig:
    """``AnalysisConfig.from_dict``, keeping the app's old leniency for a
    damage record without ``dent_depth_mm`` or ``delaminations``."""
    if isinstance(d.get("damage"), dict):
        d = {**d, "damage": {"dent_depth_mm": 0.0, "delaminations": [], **d["damage"]}}
    return AnalysisConfig.from_dict(d)


# ---------------------------------------------------------------------------
# NCR (Nonconformance Report) export — MRB disposition support
# ---------------------------------------------------------------------------
#
# ``recommend_disposition`` and ``build_ncr`` (pure: dicts in, Markdown and
# dict out) live in ``bvidfe.report.ncr`` with their unit tests (#172); the
# Export NCR tab below only collects the form and renders the result.


# ---------------------------------------------------------------------------
# Sidebar — AnalysisConfig builder
# ---------------------------------------------------------------------------

# Hydrate session_state from URL on first run, then from preset whenever the
# preset selection changes. Both must run BEFORE the widgets render so the
# widgets pick up the new values via their ``key=`` arg.
_hydrate_from_url()
_hydrate_from_preset()

if st.session_state.get("_url_rejected"):
    st.warning(
        "The shared link carried values this page cannot use; they were "
        "ignored and the sidebar keeps its defaults for them: "
        + ", ".join(f"`{r}`" for r in st.session_state["_url_rejected"])
        + "."
    )

# Seed the expert-mode widgets, then re-assign each to itself: Streamlit
# drops the session-state value of a keyed widget that is not rendered in a
# run, and these widgets do not render outside Expert mode.
for _k, _v in _EXPERT_DEFAULTS.items():
    st.session_state.setdefault(_k, _v)
for _k in (*_EXPERT_DEFAULTS, *_LOCATION_KEYS):
    if _k in st.session_state:
        st.session_state[_k] = st.session_state[_k]

with st.sidebar:
    st.markdown("**Configuration preset**")
    st.selectbox(
        "Preset",
        PRESET_NAMES,
        index=(
            PRESET_NAMES.index(st.session_state.get("preset_choice", PRESET_CUSTOM))
            if st.session_state.get("preset_choice", PRESET_CUSTOM) in PRESET_NAMES
            else 0
        ),
        key="preset_choice",
        help=(
            "Pick a canonical configuration to populate the sidebar. "
            "Choose **Custom** to keep your current settings. Selecting a "
            "preset overwrites the affected widgets; subsequent edits stick."
        ),
        on_change=_hydrate_from_preset,
    )

    expert_mode = st.toggle(
        "Expert mode",
        value=False,
        key="expert_mode",
        help=(
            "**Off (default)** — simplified sidebar: material, layup, panel, "
            "impact event, tier, loading.\n\n"
            "**On** — full controls: ply thickness, boundary condition, "
            "impactor shape, mass and impact location, and the fe3d mesh "
            "parameters expanded."
        ),
    )
    if not expert_mode:
        # Expert values set earlier, by a share link or by a preset still
        # apply with Expert mode off; say so.
        _labels = {
            "ply_thickness_mm": "ply thickness {:g} mm",
            "boundary": "boundary {}",
            "impactor_shape": "impactor shape {}",
            "mass_kg": "impactor mass {:g} kg",
        }
        _changed = [
            fmt.format(st.session_state[k])
            for k, fmt in _labels.items()
            if st.session_state[k] != _EXPERT_DEFAULTS[k]
        ]
        if not st.session_state["loc_centre"]:
            _changed.append("an off-centre impact point")
        if _changed:
            st.caption(
                "Expert-mode settings in effect: "
                + ", ".join(_changed)
                + ". Turn on Expert mode to see or change them."
            )

    st.markdown("**Material & layup**")
    if "material_choice" not in st.session_state:
        st.session_state["material_choice"] = DEFAULT_MATERIAL
    material_choice = st.selectbox(
        "Material",
        MATERIAL_NAMES,
        key="material_choice",
        help=(
            "Built-in carbon/epoxy preset from the BVID-FE material library. "
            "CAI is checked against measured coupons for AS4/8552 and "
            "IM7/8552 (calibration) and T800S/3900-2B (reported, a known "
            "miss); AS4/3501-6, T700/2510 and T800/epoxy have no measured "
            "residual-strength comparison."
        ),
    )

    if expert_mode:
        ply_thickness_mm = st.number_input(
            "Ply thickness [mm]",
            min_value=0.01,
            max_value=2.0,
            key="ply_thickness_mm",
            step=0.001,
            format="%.4f",
            help=(
                "Cured ply thickness, the same for every ply. 0.152 mm "
                "(0.006 in) is a common nominal value; NCAMP AS4/8552 tape "
                "cures to about 0.188 mm."
            ),
        )
    else:
        ply_thickness_mm = float(st.session_state["ply_thickness_mm"])

    if "layup_text" not in st.session_state:
        st.session_state["layup_text"] = "0, 45, -45, 90, 90, -45, 45, 0"
    layup_text = st.text_area(
        "Layup [deg]",
        key="layup_text",
        height=80,
        help=(
            "Comma-separated ply angles in degrees, every ply of the full "
            "stack (write symmetric halves out in full), listed from the "
            "impacted face to the back face."
        ),
    )
    try:
        layup_deg = _parse_layup(layup_text)
        if not layup_deg:
            st.error("Layup is empty — enter at least one ply angle.")
    except ValueError as exc:
        layup_deg = []
        st.error(f"Invalid layup: {exc}")

    st.markdown("**Panel geometry**")
    if "Lx_mm" not in st.session_state:
        st.session_state["Lx_mm"] = 150.0
    if "Ly_mm" not in st.session_state:
        st.session_state["Ly_mm"] = 100.0
    col_lx, col_ly = st.columns(2)
    Lx_mm = col_lx.number_input("Lx [mm]", min_value=10.0, max_value=2000.0, step=10.0, key="Lx_mm")
    Ly_mm = col_ly.number_input("Ly [mm]", min_value=10.0, max_value=2000.0, step=10.0, key="Ly_mm")
    if expert_mode:
        boundary = st.selectbox(
            "Boundary condition",
            BOUNDARIES,
            key="boundary",
            help=(
                "Panel edge support. Scales the plate bending stiffness in "
                "the Olsson onset model (clamped 2.5×, free 0.4× simply "
                "supported) and the sublaminate and panel buckling loads "
                "(clamped 1.9×, free 0.5×), and so the delamination growth "
                "stress in semi_analytical and fe3d; it does not change "
                "fe3d's first-ply-failure solve. The Soutis and Whitney-Nuismer "
                "knockdowns see it only through the damage area. 'free' is "
                "a crude approximation — interpret qualitatively."
            ),
        )
    else:
        boundary = st.session_state["boundary"]

    st.markdown("**Input mode**")
    if "input_mode" not in st.session_state:
        st.session_state["input_mode"] = "Impact event"
    input_mode = st.radio(
        "Damage source",
        ["Impact event", "C-scan inspection"],
        horizontal=False,
        key="input_mode",
        help=(
            "**Impact event** — predict damage from impact energy: Olsson's "
            "quasi-static onset threshold, then an energy-balance damage "
            "area spread over the ply interfaces and kept inside the panel. "
            "The impact is applied on the panel that is tested; to impact a "
            "support window (e.g. the ASTM D7136 125 × 75 mm opening) and "
            "test another panel, use `AnalysisConfig.impact_panel` in the "
            "Python API.\n\n"
            "**C-scan inspection** — load measured delamination ellipses "
            "directly from a JSON C-scan file."
        ),
    )

    impact_event: ImpactEvent | None = None
    damage_state: DamageState | None = None

    if input_mode == "Impact event":
        st.markdown("**Impact event**")
        if "energy_J" not in st.session_state:
            st.session_state["energy_J"] = 30.0
        if "diameter_mm" not in st.session_state:
            st.session_state["diameter_mm"] = 16.0
        energy_J = st.number_input(
            "Energy [J]",
            min_value=0.1,
            max_value=200.0,
            step=1.0,
            key="energy_J",
            help="Kinetic energy of the impactor at contact.",
        )
        diameter_mm = st.number_input(
            "Impactor diameter [mm]",
            min_value=1.0,
            max_value=100.0,
            step=1.0,
            key="diameter_mm",
            help=(
                "Impactor tip diameter; ASTM D7136 uses 16 mm. The damage "
                f"area scales with ({DIAMETER_REF_MM:g} mm / d)^{DIAMETER_DPA_EXPONENT:g} "
                "(not for a cone)."
            ),
        )
        if expert_mode:
            impactor_shape = st.selectbox(
                "Impactor shape",
                IMPACTOR_SHAPES,
                key="impactor_shape",
                help=(
                    f"Relative to a {DIAMETER_REF_MM:g} mm hemisphere, a flat-ended "
                    f"impactor leaves {SHAPE_DPA_FACTOR['flat']:g}× and a 60° cone "
                    f"{SHAPE_DPA_FACTOR['conical']:g}× the damage area (fitted to "
                    "Hasebe et al. 2025)."
                ),
            )
            mass_kg = st.number_input(
                "Impactor mass [kg]",
                min_value=0.01,
                max_value=100.0,
                key="mass_kg",
                step=0.1,
                help=(
                    "5.5 kg is the reference drop-weight mass. At the same "
                    "energy a lighter impactor leaves slightly more damage "
                    "((5.5 / m)^0.1); below the plate's own mass the "
                    "quasi-static model may underpredict damage."
                ),
            )
            st.checkbox(
                "Centre of panel",
                key="loc_centre",
                help="Impact the panel centre. Untick to place the impact point.",
            )
            if not st.session_state["loc_centre"]:
                # Seed at the centre; keep a stored point on the panel when
                # Lx / Ly shrink (the inputs are bounded to it).
                for _key, _span in (("loc_x_mm", float(Lx_mm)), ("loc_y_mm", float(Ly_mm))):
                    _cur = float(st.session_state.get(_key, _span / 2.0))
                    st.session_state[_key] = min(max(_cur, 0.0), _span)
                col_lx_loc, col_ly_loc = st.columns(2)
                _loc_help = (
                    "Impact point measured from the panel corner, on the "
                    "panel. The delaminations are kept inside the panel "
                    "around this point."
                )
                col_lx_loc.number_input(
                    "Location X [mm]",
                    min_value=0.0,
                    max_value=float(Lx_mm),
                    step=1.0,
                    key="loc_x_mm",
                    help=_loc_help,
                )
                col_ly_loc.number_input(
                    "Location Y [mm]",
                    min_value=0.0,
                    max_value=float(Ly_mm),
                    step=1.0,
                    key="loc_y_mm",
                    help=_loc_help,
                )
        else:
            impactor_shape = st.session_state["impactor_shape"]
            mass_kg = float(st.session_state["mass_kg"])
        # None is the library's panel centre.
        loc_xy: tuple[float, float] | None = None
        if not st.session_state["loc_centre"]:
            loc_xy = (
                min(max(float(st.session_state.get("loc_x_mm", Lx_mm / 2.0)), 0.0), float(Lx_mm)),
                min(max(float(st.session_state.get("loc_y_mm", Ly_mm / 2.0)), 0.0), float(Ly_mm)),
            )
        if loc_xy == (0.0, 0.0):
            # The library still reads (0, 0) as its deprecated centre
            # sentinel, so it must not get the corner (#171).
            st.error(
                "Invalid impact event: the impact point (0, 0) is the panel "
                "corner; it must lie strictly inside the panel. Move it off "
                "the corner, or tick *Centre of panel*."
            )
        else:
            try:
                impact_event = ImpactEvent(
                    energy_J=energy_J,
                    impactor=ImpactorGeometry(diameter_mm=diameter_mm, shape=impactor_shape),
                    mass_kg=mass_kg,
                    location_xy_mm=loc_xy,
                )
            except (TypeError, ValueError) as exc:
                st.error(f"Invalid impact event: {exc}")
    else:
        st.markdown("**C-scan inspection**")
        uploaded = st.file_uploader(
            "C-scan JSON",
            type=["json"],
            help="JSON conforming to docs/cscan_schema.md.",
        )
        override_dent = st.checkbox(
            "Override file value: dent depth",
            key="override_dent",
            help="Use the dent depth below instead of the C-scan file's.",
        )
        dent_depth_mm = st.number_input(
            "Dent depth [mm]",
            min_value=0.0,
            max_value=10.0,
            value=0.0,
            step=0.05,
            disabled=not override_dent,
            help=(
                "Replaces the file's dent depth (0 included) when "
                "*Override file value* is ticked. Reported and used by the "
                "NCR detectability check; it does not change the predicted "
                "strength."
            ),
        )
        override_fb = st.checkbox(
            "Override file value: fiber-break radius",
            key="override_fb",
            help="Use the fiber-break radius below instead of the C-scan file's.",
        )
        fb_radius_mm = st.number_input(
            "Fiber-break radius [mm]",
            min_value=0.0,
            max_value=50.0,
            value=0.0,
            step=0.5,
            disabled=not override_fb,
            help=(
                "Replaces the file's fiber-break core radius (0 included) "
                "when *Override file value* is ticked. Softens the fe3d mesh "
                "around the delamination centroids; the closed-form tiers "
                "ignore it."
            ),
        )
        if uploaded is not None:
            try:
                # Parse in memory — never write an attacker-controlled
                # filename to disk (path-traversal write primitive).
                data = json.loads(uploaded.getvalue())
                damage_state = damage_state_from_dict(data)
                # Override dent / fiber-break with the sidebar inputs when
                # asked to (0 included); otherwise keep the file's values.
                if override_dent:
                    damage_state = replace(damage_state, dent_depth_mm=float(dent_depth_mm))
                if override_fb:
                    damage_state = replace(damage_state, fiber_break_radius_mm=float(fb_radius_mm))
                st.success(f"Loaded {len(damage_state.delaminations)} delamination(s)")
            except (CScanSchemaError, json.JSONDecodeError, ValueError) as exc:
                st.error(f"Invalid C-scan JSON: {exc}")
                damage_state = None
        else:
            st.info(
                "Upload a C-scan JSON to populate the damage state. Until "
                "then the run button is disabled."
            )
            damage_state = None

    st.markdown("**Analysis**")
    if "tier" not in st.session_state:
        st.session_state["tier"] = TIERS[0]
    if "loading" not in st.session_state:
        st.session_state["loading"] = "compression"
    tier = st.selectbox(
        "Tier",
        TIERS,
        key="tier",
        help=(
            "**empirical** — closed form: Soutis knockdown on the damage "
            "area (CAI), Whitney-Nuismer open-hole equivalent (TAI). "
            "Milliseconds.\n\n"
            "**semi_analytical** — CAI is the lower of the Soutis bound and "
            "the stress at which the weakest buckled sublaminate grows its "
            "delamination; TAI as empirical. Well under a second.\n\n"
            "Both closed-form tiers: CAI calibrated on measured NCAMP "
            "AS4/8552 and Lovejoy & Scotti IM7/8552 coupons, about 6-9% "
            "error when fitted on one dataset and checked on the other. TAI "
            "is not validated against measured data.\n\n"
            "**fe3d** — 3D hex FE, first-ply failure (LaRC05 for CAI, "
            "Tsai-Wu for TAI) with a softened damage zone. Qualitative: the "
            "damage-zone stiffness is fitted in-sample to three measured CAI "
            "coupons (Lovejoy & Scotti within 1.2%, NCAMP 21% low). In TAI "
            "the zone keeps its fibres and the broken fibres are the closed-form "
            "equivalent hole, so fe3d TAI falls with energy 2-8% below the "
            "closed-form TAI (3.6% off the one measured TAI dataset, partly "
            "in-sample). About 5 s "
            "for an 8-ply 150 x 100 mm panel at the 5 mm default mesh, longer "
            "for thick laminates "
            f"(solves damaged and undamaged panels); capped at {FE3D_MAX_DOF:,} "
            "DOF."
        ),
    )
    loading = st.radio(
        "Loading mode",
        ["compression", "tension"],
        horizontal=True,
        key="loading",
        help="CAI = compression-after-impact, TAI = tension-after-impact.",
    )

    mesh_params: MeshParams | None = None
    if tier == "fe3d":
        with st.expander("Mesh parameters (fe3d only)", expanded=expert_mode):
            elements_per_ply = st.number_input(
                "Elements per ply", min_value=1, max_value=10, value=1, step=1
            )
            in_plane_size_mm = st.number_input(
                "In-plane size [mm]",
                min_value=0.5,
                max_value=50.0,
                value=5.0,
                step=0.5,
                help=(
                    "Target in-plane element edge length. Smaller values = "
                    "finer mesh, quadratically more elements, memory and "
                    "time. The 5 mm default resolves damage areas above "
                    "about 1,000 mm² (at ~1,500 mm² it is within ~3-10% of a "
                    "2 mm mesh); a smaller damage zone is only a few "
                    "elements across and its knockdown can be ~25% high, so "
                    "use 2.5 mm or finer (an 8-ply 150 x 100 mm panel takes "
                    "~35 s at 2.5 mm). An 8-ply 150 x 100 mm panel "
                    "peaks at ~0.2 GB at 5 mm and ~0.5 GB at 2.5 mm (more if "
                    "the solver falls back to a direct factorisation). On "
                    "Streamlit Community Cloud (1 GB RAM) keep this >= 5 mm."
                ),
            )
            mesh_params = MeshParams(
                elements_per_ply=int(elements_per_ply),
                in_plane_size_mm=float(in_plane_size_mm),
            )

    run_clicked = st.button("Run analysis", type="primary", width="stretch")


# Sync the live sidebar config into the URL so the page link is shareable.
# Runs once per script execution, after all widgets have committed their
# values to ``st.session_state``.
_sync_query_params()


# ---------------------------------------------------------------------------
# Validate config and live-preview E_onset before running
# ---------------------------------------------------------------------------

panel = None
config_dict: dict | None = None
config_valid = bool(layup_deg) and Lx_mm > 0 and Ly_mm > 0
if input_mode == "C-scan inspection" and damage_state is None:
    config_valid = False
if input_mode == "Impact event" and impact_event is None:
    config_valid = False

if config_valid:
    panel = PanelGeometry(Lx_mm=float(Lx_mm), Ly_mm=float(Ly_mm), boundary=boundary)
    config_dict = {
        "material": material_choice,
        "layup_deg": layup_deg,
        "ply_thickness_mm": ply_thickness_mm,
        "panel": {
            "Lx_mm": panel.Lx_mm,
            "Ly_mm": panel.Ly_mm,
            "boundary": panel.boundary,
        },
        "loading": loading,
        "tier": tier,
        "impact": (
            None
            if impact_event is None
            else {
                "energy_J": impact_event.energy_J,
                "impactor": {
                    "diameter_mm": impact_event.impactor.diameter_mm,
                    "shape": impact_event.impactor.shape,
                },
                "mass_kg": impact_event.mass_kg,
                "location_xy_mm": (
                    None
                    if impact_event.location_xy_mm is None
                    else list(impact_event.location_xy_mm)
                ),
            }
        ),
        "damage": (
            None
            if damage_state is None
            else {
                "dent_depth_mm": damage_state.dent_depth_mm,
                "fiber_break_radius_mm": damage_state.fiber_break_radius_mm,
                "delaminations": [
                    {
                        "interface_index": e.interface_index,
                        "centroid_mm": list(e.centroid_mm),
                        "major_mm": e.major_mm,
                        "minor_mm": e.minor_mm,
                        "orientation_deg": e.orientation_deg,
                    }
                    for e in damage_state.delaminations
                ],
            }
        ),
        "mesh": (
            None
            if mesh_params is None
            else {
                "elements_per_ply": mesh_params.elements_per_ply,
                "in_plane_size_mm": mesh_params.in_plane_size_mm,
            }
        ),
    }


# Persistent result across reruns.
if "last_result" not in st.session_state:
    st.session_state["last_result"] = None


def _live_onset_preview() -> tuple[str | None, str | None]:
    """Olsson onset energy + DPA preview for the impact-driven path."""
    if input_mode != "Impact event" or panel is None or impact_event is None:
        return None, None
    return _onset_preview_cached(
        material_choice,
        tuple(layup_deg),
        float(ply_thickness_mm),
        float(panel.Lx_mm),
        float(panel.Ly_mm),
        panel.boundary,
        float(impact_event.energy_J),
        float(impact_event.impactor.diameter_mm),
        impact_event.impactor.shape,
        float(impact_event.mass_kg),
        (
            None
            if impact_event.location_xy_mm is None
            else (float(impact_event.location_xy_mm[0]), float(impact_event.location_xy_mm[1]))
        ),
    )


preview, preview_error = _live_onset_preview()
if preview:
    st.info(preview)
elif preview_error:
    st.caption(f"Onset preview unavailable — {preview_error}")


# ---------------------------------------------------------------------------
# Run analysis
# ---------------------------------------------------------------------------

if run_clicked:
    if not config_valid or config_dict is None:
        st.error(
            "Configuration is incomplete — fix the highlighted sidebar inputs " "before running."
        )
    else:
        # fe3d mesh-size guard: estimate DOF count before the long-running
        # solve so we can refuse impossibly large meshes early.
        if tier == "fe3d":
            try:
                cfg_check = _config_from_dict(config_dict)
                stats = estimate_fe_mesh_size(cfg_check)
            except Exception as exc:
                st.error(f"Failed to estimate mesh size: {exc}")
                stats = None
            if stats and stats["n_dof"] > FE3D_MAX_DOF:
                st.error(
                    f"fe3d mesh has {stats['n_elements']:,} elements "
                    f"({stats['n_dof']:,} DOFs), exceeding the safe cap "
                    f"of {FE3D_MAX_DOF:,}. Increase **In-plane size** or "
                    f"reduce **Elements per ply**."
                )
                st.stop()
        with st.spinner("Running BVID analysis…"):
            try:
                # Time the cached call so we can surface cache-hit vs
                # cache-miss as a subtle toast (#112). ``st.cache_data``
                # returns its memoised payload without entering the wrapped
                # function, so anything < ~50 ms is a hit in practice
                # (a real BvidAnalysis solve is at least sub-second on the
                # empirical tier and seconds-to-minutes on fe3d).
                _t0 = time.perf_counter()
                payload = _run_analysis_cached(config_dict)
                _dt = time.perf_counter() - _t0
                st.session_state["last_result"] = payload
                try:
                    if _dt < 0.05:
                        st.toast("Cache hit — instant", icon="⚡")
                    else:
                        st.toast(f"Solved in {_dt:.1f}s", icon="✓")
                except Exception:
                    # ``st.toast`` is fairly new; never let UX sugar
                    # block a successful analysis from rendering.
                    pass
            except Exception as exc:
                st.session_state["last_result"] = None
                st.error(f"Analysis failed: {exc}")


# ---------------------------------------------------------------------------
# Results tabs
# ---------------------------------------------------------------------------

payload = st.session_state.get("last_result")

if payload is not None and config_dict is not None and config_dict != payload["config_dict"]:
    st.warning(
        "Sidebar inputs changed since the last run — the results below are "
        "stale. Click **Run analysis** to refresh."
    )

tabs = st.tabs(
    [
        "Summary",
        "Damage Map",
        "Knockdown Curve",
        "3D Damage",
        "Buckling",
        "Damage Severity",
        "Sweep",
        "Export NCR",
    ]
)


# --- Summary tab ----------------------------------------------------------

with tabs[0]:
    st.subheader("Summary")
    if payload is None:
        st.info("Run an analysis to see its summary here.")
    else:
        result_dict = payload["result_dict"]
        col1, col2, col3 = st.columns(3)
        col1.metric(
            "Pristine strength",
            f"{result_dict['pristine_strength_MPa']:.1f} MPa",
        )
        col2.metric(
            "Residual strength",
            f"{result_dict['residual_strength_MPa']:.1f} MPa",
        )
        col3.metric("Knockdown", f"{result_dict['knockdown']:.3f}")
        st.text(payload["summary_text"])
        if payload["notes"]:
            st.warning("Runtime notes:\n\n" + "\n\n".join(payload["notes"]))

        st.download_button(
            "Download results JSON",
            data=json.dumps(result_dict, indent=2),
            file_name="bvidfe_results.json",
            mime="application/json",
        )


# --- Damage Map tab -------------------------------------------------------

with tabs[1]:
    st.subheader("Damage map")
    if payload is None:
        st.info("Run an analysis to see the damage map.")
    else:
        try:
            cfg = payload["config"]
            fig = plot_damage_map(payload["damage"], cfg.panel)
            st.pyplot(fig, width="content")
        except Exception as exc:
            st.error(f"Could not render the damage map: {exc}")


# --- Knockdown Curve tab --------------------------------------------------

with tabs[2]:
    st.subheader("Knockdown vs impact energy")
    if payload is None:
        st.info(
            "Run an impact-driven analysis to seed this view, or use the "
            "Sweep tab to compute a full curve."
        )
    elif payload["config"].impact is None:
        st.info(
            "Knockdown curves require an impact-driven configuration. "
            "Switch the sidebar to **Impact event** mode."
        )
    else:
        cfg = payload["config"]
        run_tier = payload["result_dict"]["tier_used"]
        # No key: the default follows the tier of the latest run.
        curve_tier = st.selectbox(
            "Tier for curve",
            TIERS,
            index=TIERS.index(run_tier) if run_tier in TIERS else 0,
            help=_SWEEP_TIER_HELP,
        )
        # 7-point sweep around the analysed energy, cached so it doesn't
        # rerun every script pass.
        sweep_E = _knockdown_curve_energies(cfg.impact.energy_J)
        curve_ok = True
        if curve_tier == "fe3d":
            st.warning(f"fe3d curve: {len(sweep_E)} fe3d analyses. {_SWEEP_TIER_HELP}")
            try:
                stats = estimate_fe_mesh_size(replace(cfg, tier="fe3d"))
            except Exception as exc:
                st.error(f"Failed to estimate fe3d mesh size: {exc}")
                curve_ok = False
            else:
                if stats["n_dof"] > FE3D_MAX_DOF:
                    st.error(
                        f"fe3d mesh has {stats['n_elements']:,} elements "
                        f"({stats['n_dof']:,} DOFs), exceeding the safe cap "
                        f"of {FE3D_MAX_DOF:,}. Increase **In-plane size** or "
                        f"reduce **Elements per ply** in the sidebar."
                    )
                    curve_ok = False
        if curve_ok:
            try:
                with st.spinner(f"Computing {curve_tier} knockdown sweep…"):
                    df = _sweep_cached(payload["config_dict"], curve_tier, sweep_E)
                ok = df.dropna(subset=["knockdown"])
                fig = plot_knockdown_curve(
                    ok["energy_J"].tolist(),
                    ok["knockdown"].tolist(),
                    tier_label=curve_tier,
                )
                st.pyplot(fig, width="content")
                _show_df(df)
            except Exception as exc:
                st.error(f"Could not compute the knockdown curve: {exc}")


# --- 3D Damage tab --------------------------------------------------------

with tabs[3]:
    st.subheader("3D damaged mesh")
    if payload is None:
        st.info("Run an analysis to see the 3D damaged mesh.")
    else:
        mesh = payload["mesh"]
        if mesh is None:
            st.warning("The display mesh could not be built for this configuration.")
        else:
            try:
                st.plotly_chart(
                    mesh_damage_figure(mesh, title="Out-of-plane stiffness loss (0 = pristine)"),
                    width="stretch",
                )
                st.caption(
                    f"Mesh: {mesh.n_elements:,} elements, {mesh.n_nodes:,} nodes. "
                    "Red elements touch a delaminated interface inside its "
                    "ellipse, or sit in the fiber-break core. The in-plane "
                    "softening of the whole damage zone (fe3d) is not shown."
                )
            except Exception as exc:
                st.error(f"Could not render the 3D mesh: {exc}")


# --- Buckling tab ---------------------------------------------------------

with tabs[4]:
    st.subheader("Buckling indicator")
    if payload is None:
        st.info("Run an analysis to see the buckling indicator.")
    else:
        eigs = payload["buckling_eigenvalues"]
        tier_used = payload["result_dict"]["tier_used"]
        if not eigs:
            st.info(
                f"Tier '{tier_used}' does not produce a buckling indicator. "
                "Switch to **semi_analytical** or **fe3d** to populate this "
                "tab."
            )
        else:
            # Units depend on tier:
            # - semi_analytical: N_cr (force per unit sublaminate width, N/mm)
            # - fe3d: sigma_crit (compressive buckling stress, MPa) via the
            #   Rayleigh-Ritz closed form (#129).
            if tier_used == "fe3d":
                col_label, unit_note = (
                    "σ_crit [MPa]",
                    "fe3d: the lower of whole-panel and sublaminate buckling "
                    "stress, from the Rayleigh-Ritz closed form (#129). "
                    "Buckling onset only; it does not set the strength.",
                )
            else:
                col_label, unit_note = (
                    "N_cr [N/mm]",
                    "semi_analytical: buckling force per unit width of the "
                    "sublaminate whose delamination grows first. Buckling "
                    "onset only; the growth stress sets the strength.",
                )
            st.caption(unit_note)
            df_eigs = pd.DataFrame({"mode": list(range(1, len(eigs) + 1)), col_label: eigs})
            st.bar_chart(df_eigs.set_index("mode"))
            st.dataframe(
                df_eigs.style.format({col_label: "{:.4g}"}),
                width="stretch",
            )


# --- Damage Severity tab --------------------------------------------------

with tabs[5]:
    st.subheader("Damage severity (through-thickness sum)")
    if payload is None:
        st.info("Run an analysis to see the damage-severity heatmap.")
    else:
        mesh = payload["mesh"]
        if mesh is None:
            st.warning("The display mesh could not be built for this configuration.")
        else:
            try:
                cfg_mesh = payload["mesh_config"]
                grid = element_grid_shape(cfg_mesh)
                severity = (1.0 - mesh.damage_factors).reshape(grid).sum(axis=0)
                fig, ax = plt.subplots(figsize=(6, 5))
                im = ax.imshow(
                    severity,
                    origin="lower",
                    extent=(0, cfg_mesh.panel.Lx_mm, 0, cfg_mesh.panel.Ly_mm),
                    cmap="hot_r",
                    vmin=0.0,
                    vmax=max(1.0, float(severity.max())),
                )
                ax.set_xlabel("x [mm]")
                ax.set_ylabel("y [mm]")
                ax.set_title(
                    f"Damage severity — tier={payload['result_dict']['tier_used']}  "
                    f"KD={payload['result_dict']['knockdown']:.3f}"
                )
                ax.set_aspect("equal")
                fig.colorbar(im, ax=ax, label="Stacked OOP loss")
                fig.tight_layout()
                st.pyplot(fig, width="content")
            except Exception as exc:
                st.error(f"Could not render the damage-severity heatmap: {exc}")


# --- Sweep tab ------------------------------------------------------------

with tabs[6]:
    st.subheader("Parametric energy sweep")
    if not config_valid or config_dict is None:
        st.info("Fix the sidebar inputs before running a sweep.")
    elif input_mode != "Impact event":
        st.info("Energy sweeps require **Impact event** mode in the sidebar.")
    else:
        default_E = ", ".join(f"{v:g}" for v in [5, 10, 20, 30, 40])
        energies_text = st.text_input("Energies [J] (comma-separated)", value=default_E)
        sweep_tier = st.selectbox(
            "Tier for sweep",
            TIERS,
            index=0,
            help=_SWEEP_TIER_HELP,
        )
        sweep_clicked = st.button("Run sweep", type="primary")
        if sweep_clicked:
            try:
                energies = [float(x.strip()) for x in energies_text.split(",") if x.strip()]
            except ValueError:
                st.error("Energies must be a comma-separated list of numbers.")
                energies = []
            if energies and not all(math.isfinite(e) and e > 0 for e in energies):
                st.error("Energies must be finite and positive.")
                energies = []
            if len(energies) > MAX_SWEEP_POINTS:
                st.error(
                    f"Sweep is capped at {MAX_SWEEP_POINTS} energy points "
                    f"(got {len(energies)}). Shorten the list."
                )
                energies = []
            if energies and sweep_tier == "fe3d":
                try:
                    stats = estimate_fe_mesh_size(
                        replace(_config_from_dict(config_dict), tier="fe3d")
                    )
                except Exception as exc:
                    st.error(f"Failed to estimate fe3d mesh size: {exc}")
                    energies = []
                else:
                    if stats["n_dof"] > FE3D_MAX_DOF:
                        st.error(
                            f"fe3d mesh has {stats['n_elements']:,} elements "
                            f"({stats['n_dof']:,} DOFs), exceeding the safe cap "
                            f"of {FE3D_MAX_DOF:,}. Increase **In-plane size** or "
                            f"reduce **Elements per ply** in the sidebar."
                        )
                        energies = []
            if energies:
                # Live per-point progress (#112). The bar is wired into a
                # callback passed to ``sweep_energies`` so each iteration
                # updates the UI as it lands; ``progress_callback`` is
                # ignored by ``st.cache_data``-style memoisation so we call
                # the sweep entrypoint directly here. A buggy callback
                # cannot abort the sweep itself — the library catches it.
                prog = st.progress(0.0, text=f"0/{len(energies)} starting…")

                def _on_point(i: int, total: int, inputs: dict) -> None:
                    e_val = inputs.get("energy_J", "?")
                    e_txt = f"{e_val:g}" if isinstance(e_val, (int, float)) else str(e_val)
                    prog.progress(i / total, text=f"{i}/{total}: E={e_txt} J")

                try:
                    with st.spinner(f"Running {len(energies)}-point {sweep_tier} sweep…"):
                        cfg_sweep = replace(_config_from_dict(config_dict), tier=sweep_tier)
                        df = sweep_energies(
                            cfg_sweep,
                            list(energies),
                            on_error="warn",
                            progress_callback=_on_point,
                        )
                except Exception as exc:
                    st.error(f"Sweep failed: {exc}")
                else:
                    prog.empty()
                    # Failed-point summary (#112). The fixed-schema sweep
                    # DataFrame (#38) always carries an ``error`` column;
                    # surface failed rows in a small table above the success
                    # plot so the user sees which points to investigate.
                    if "error" in df.columns:
                        failed_df = df[df["error"].notna() & (df["error"] != "")]
                    else:  # defensive — schema is fixed by parametric_sweep
                        failed_df = df.loc[df["knockdown"].isna()]
                    if df["knockdown"].isna().all():
                        st.error(
                            "Every sweep point failed — check the sidebar "
                            "configuration and the energy values."
                        )
                        if not failed_df.empty:
                            st.dataframe(
                                failed_df[["energy_J", "error"]],
                                width="stretch",
                            )
                    else:
                        if not failed_df.empty:
                            st.warning(f"⚠ {len(failed_df)} point(s) failed:")
                            cols_to_show = [
                                c for c in ("energy_J", "error") if c in failed_df.columns
                            ]
                            st.dataframe(
                                failed_df[cols_to_show],
                                width="stretch",
                            )
                        _show_df(df)
                        ok = df.dropna(subset=["knockdown"])
                        fig = plot_knockdown_curve(
                            ok["energy_J"].tolist(),
                            ok["knockdown"].tolist(),
                            tier_label=sweep_tier,
                        )
                        st.pyplot(fig, width="content")
                        csv_buf = io.StringIO()
                        df.to_csv(csv_buf, index=False)
                        st.download_button(
                            "Download sweep CSV",
                            data=csv_buf.getvalue(),
                            file_name=f"bvidfe_sweep_{sweep_tier}.csv",
                            mime="text/csv",
                        )


# --- Export NCR tab -------------------------------------------------------

with tabs[7]:
    st.subheader("Create Nonconformance Report (NCR)")
    st.caption(
        "For an engineer in the field: log the nonconformance, attach the "
        "BVID-FE residual-strength analysis, and export a structured NCR "
        "with a **proposed** MRB disposition. The recommendation is "
        "decision support only — a qualified Material Review Board must "
        "review, modify, and approve it."
    )

    # The NCR needs a material/layup/panel configuration. Prefer the last
    # analysed config so results attach; fall back to the live sidebar
    # config so the NC can be logged before analysis.
    ncr_payload = st.session_state.get("last_result")
    if ncr_payload is not None:
        ncr_cfg = ncr_payload["config"]
        ncr_result = ncr_payload["result_dict"]
    elif config_dict is not None:
        ncr_cfg = _config_from_dict(config_dict)
        ncr_result = None
    else:
        ncr_cfg = None
        ncr_result = None

    if ncr_cfg is None:
        st.info(
            "Complete the sidebar configuration (material, layup, panel) " "before creating an NCR."
        )
    else:
        if ncr_result is None:
            st.warning(
                "No BVID-FE analysis attached — the NCR will record the "
                "nonconformance with engineering analysis marked pending. "
                "Run an analysis to attach residual strength and a "
                "disposition recommendation."
            )

        today = date.today()
        c1, c2, c3 = st.columns(3)
        with c1:
            ncr_number = st.text_input(
                "NCR number",
                value=f"NCR-{today:%Y%m%d}-001",
                key="ncr_number",
            )
            revision = st.text_input("Revision", value="-", key="ncr_rev")
            originator = st.text_input("Originated by (engineer)", value="", key="ncr_originator")
            ncr_date = st.date_input("Date", value=today, key="ncr_date")
        with c2:
            program = st.text_input("Program / project", value="", key="ncr_program")
            part_number = st.text_input("Part number", value="", key="ncr_pn")
            part_name = st.text_input("Part name", value="", key="ncr_part_name")
            serial_lot = st.text_input("Serial / lot", value="", key="ncr_serial")
        with c3:
            work_order = st.text_input("Work order", value="", key="ncr_wo")
            quantity = st.text_input("Quantity affected", value="1", key="ncr_qty")
            station = st.text_input("Found at / station", value="", key="ncr_station")
            found_during = st.selectbox(
                "Found during", FOUND_DURING, index=3, key="ncr_found_during"
            )

        nc_type = st.selectbox("Nonconformance type", NC_TYPES, index=0, key="ncr_nc_type")
        dimensions = st.text_input(
            "As-found dimensions / location",
            value="",
            placeholder="e.g. 25 mm dia. dent, 0.4 mm deep, 120 mm from edge",
            key="ncr_dims",
        )
        description = st.text_area(
            "Nonconformance description",
            value="",
            height=120,
            placeholder=(
                "Describe what was found, how, and any inspection (visual, "
                "tap-test, ultrasonic C-scan) performed."
            ),
            key="ncr_desc",
        )
        dent_threshold_mm = st.number_input(
            "BVID detectability threshold [mm]",
            min_value=0.05,
            max_value=5.0,
            value=BVID_DENT_THRESHOLD_MM,
            step=0.05,
            help=(
                "Residual dent depth above which damage is treated as "
                "visible (VID). Override with your program's "
                "damage-tolerance value."
            ),
            key="ncr_dent_thr",
        )

        form = {
            "ncr_number": ncr_number.strip() or f"NCR-{today:%Y%m%d}-001",
            "revision": revision.strip() or "-",
            "originator": originator.strip() or "—",
            "date": ncr_date.isoformat(),
            "program": program.strip() or "—",
            "part_number": part_number.strip() or "—",
            "part_name": part_name.strip() or "—",
            "serial_lot": serial_lot.strip() or "—",
            "work_order": work_order.strip() or "—",
            "quantity": quantity.strip() or "—",
            "station": station.strip() or "—",
            "found_during": found_during,
            "nc_type": nc_type,
            "dimensions": dimensions.strip() or "—",
            "description": description.strip(),
        }

        rec = (
            None
            if ncr_result is None
            else recommend_disposition(ncr_result, float(dent_threshold_mm))
        )
        ncr_md, ncr_dict = build_ncr(form, ncr_cfg, ncr_result, rec)

        if rec is not None:
            st.success(
                f"Proposed disposition path: **{rec['recommended_path']}** " f"— {rec['status']}"
            )
            if rec["flags"]:
                st.warning("Flags:\n\n" + "\n\n".join(f"- {f}" for f in rec["flags"]))

        with st.expander("Preview NCR", expanded=True):
            st.markdown(ncr_md)

        safe_id = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in form["ncr_number"])
        dl1, dl2 = st.columns(2)
        with dl1:
            st.download_button(
                "Download NCR (Markdown)",
                data=ncr_md,
                file_name=f"{safe_id}.md",
                mime="text/markdown",
                width="stretch",
            )
        with dl2:
            st.download_button(
                "Download NCR (JSON)",
                data=json.dumps(ncr_dict, indent=2),
                file_name=f"{safe_id}.json",
                mime="application/json",
                width="stretch",
            )
