"""Parametric sweep utilities producing pandas DataFrames and CSV output.

All four sweep entry points (`sweep_energies`, `sweep_layups`,
`sweep_thicknesses`, `sweep_transverse_ratios`) accept two optional
parameters:

  * ``on_error`` — ``"raise"`` (default), ``"skip"``, or ``"warn"``.
    With ``"raise"`` a per-iteration failure aborts the whole sweep and
    no CSV is written (the legacy behaviour). With ``"skip"`` the failed
    row is filled with ``NaN`` numeric fields plus an ``error`` column
    containing the exception message, and the sweep continues; partial
    results still reach the DataFrame and the CSV. ``"warn"`` is identical
    to ``"skip"`` but additionally emits a ``UserWarning`` per failure.
    A point whose config fails validation (e.g. a per-ply thickness list
    that does not match a swept layup, or a negative energy) is a failed
    point like any other. An unknown ``on_error`` raises ``ValueError``
    before any point runs, even for an empty sweep.

  * ``progress_callback`` — optional callable invoked once per iteration
    to report progress. Two signatures are accepted (introspected at
    runtime, so existing callers don't need to change):

      - ``progress_callback(i_done, n_total)`` — legacy 2-arg form.
      - ``progress_callback(i_done, n_total, inputs)`` — 3-arg form where
        ``inputs`` is a small ``dict`` describing the iteration's swept
        inputs (``energy_J``, ``layup``, ``ply_thickness_mm`` or
        ``transverse_stress_ratio``). Used by
        the Streamlit UI to render a live ``E=…J`` progress label
        (#112).

    Callback exceptions are caught and logged at DEBUG so a misbehaving
    UI never aborts the sweep itself.

Each returns a ``pandas.DataFrame`` with the swept input first
(``energy_J``, ``layup`` as a ``"0/45/-45/90"`` string, or
``ply_thickness_mm``) followed by ``knockdown``, ``residual_MPa``,
``pristine_MPa``, ``dpa_mm2``, ``dent_mm``, ``n_delaminations``,
``tier_used`` and ``error`` (NaN on success). A failed point keeps
``n_delaminations = 0`` and ``tier_used = ""`` with NaN numbers. With
``csv_path`` the frame is also written as CSV (no index). Each impact
warning category (``SmallMassQuasiStaticWarning``, ``DPACapClipWarning``) is
shown at most once per sweep.

``sweep_transverse_ratios`` (biaxial loading, experimental and not
validated; issue #115) puts ``transverse_stress_ratio`` first and adds
``residual_sigma_y_MPa``, ``pristine_sigma_y_MPa`` and
``biaxial_interaction_factor`` before ``error``; the other sweeps' columns
are unchanged.

Together these make sweeps robust to a single bad input (a degenerate
mesh, a Tsai-Wu invalid combination, an out-of-range Olsson regime) and
preserve the partial CSV that previously evaporated when the sweep
aborted mid-run.
"""

from __future__ import annotations

import inspect
import logging
import math
import warnings
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable, Iterator, List, Optional, Sequence, Union

import pandas as pd

from bvidfe._types import _ON_ERROR_MODES, OnErrorMode
from bvidfe._validation import require_choice
from bvidfe.analysis import AnalysisConfig, AnalysisResults, BvidAnalysis
from bvidfe.impact.mapping import DPACapClipWarning, SmallMassQuasiStaticWarning

logger = logging.getLogger(__name__)

# Accept both 2-arg ``(i_done, n_total)`` and 3-arg
# ``(i_done, n_total, inputs)`` progress callbacks — the latter (#112) carries
# a small dict of the iteration's swept inputs so a UI can render a live
# per-point label without re-implementing the loop.
_ProgressCallback = Union[
    Callable[[int, int], None],
    Callable[[int, int, dict], None],
]


_DEDUP_WARNING_CATEGORIES = (SmallMassQuasiStaticWarning, DPACapClipWarning)


@contextmanager
def _dedupe_impact_warnings() -> Iterator[None]:
    """Emit each impact-mapping diagnostic at most once for the whole sweep.

    ``impact_to_damage`` raises ``SmallMassQuasiStaticWarning`` /
    ``DPACapClipWarning`` once per sweep point. Python's ``"once"`` /
    ``"default"`` filter actions key on the *message text*, which here
    embeds per-point values (the predicted DPA, the mass ratio), so the
    filter never recognises them as duplicates and the diagnostic floods
    stderr once per iteration. Dedupe by *category* explicitly instead.
    Non-impact warnings (e.g. the ``on_error="warn"`` per-failure
    ``UserWarning``) pass through unchanged.
    """
    seen: set[type] = set()
    with warnings.catch_warnings():
        outer_show = warnings.showwarning
        # See every warning so we can decide per category; the original
        # showwarning is still used for the first of each / pass-through.
        warnings.simplefilter("always")

        def _show(message, category, filename, lineno, file=None, line=None):
            if issubclass(category, _DEDUP_WARNING_CATEGORIES):
                if category in seen:
                    return
                seen.add(category)
            outer_show(message, category, filename, lineno, file, line)

        warnings.showwarning = _show
        yield


def _run_one(cfg: AnalysisConfig) -> dict:
    """Run a single analysis and return a dict of key fields."""
    return _result_row(BvidAnalysis(cfg).run())


def _result_row(result: AnalysisResults) -> dict:
    """The per-run columns of one result (``_RESULT_COLS`` but ``error``)."""
    return {
        "knockdown": result.knockdown,
        "residual_MPa": result.residual_strength_MPa,
        "pristine_MPa": result.pristine_strength_MPa,
        "dpa_mm2": result.dpa_mm2,
        "dent_mm": result.damage.dent_depth_mm,
        "n_delaminations": len(result.damage.delaminations),
        "tier_used": result.tier_used,
    }


def _nan_row() -> dict:
    """Row schema matching _run_one with NaN numerics — used for failed
    iterations under on_error in {'skip', 'warn'}."""
    return {
        "knockdown": math.nan,
        "residual_MPa": math.nan,
        "pristine_MPa": math.nan,
        "dpa_mm2": math.nan,
        "dent_mm": math.nan,
        "n_delaminations": 0,
        "tier_used": "",
    }


def _check_on_error(on_error: str) -> None:
    """Reject an unknown ``on_error`` before the sweep starts, so a typo is
    caught even when the sweep list is empty or every point succeeds."""
    require_choice(on_error, _ON_ERROR_MODES, "on_error")


def _try_run_one(
    build_cfg: Callable[[], AnalysisConfig],
    on_error: OnErrorMode,
    label: str,
    run: Optional[Callable[[AnalysisConfig], dict]] = None,
    nan_row: Callable[[], dict] = _nan_row,
) -> dict:
    """Build the point's config with ``build_cfg()`` and run it, translating
    any exception from either step to a NaN-filled row according to
    ``on_error``. Building inside the handler means a point that fails
    ``AnalysisConfig`` / ``ImpactEvent`` validation is recorded like a
    run-time failure. ``label`` is used in the warning text. ``run`` maps
    the config to the row (``_run_one`` when ``None``, looked up at call
    time) and ``nan_row`` gives a failed row's columns."""
    try:
        return (run or _run_one)(build_cfg())
    except Exception as exc:  # noqa: BLE001
        if on_error == "raise":
            raise
        row = nan_row()
        row["error"] = f"{type(exc).__name__}: {exc}"
        if on_error == "warn":
            warnings.warn(
                f"sweep iteration {label} failed: {row['error']}",
                stacklevel=3,
            )
        return row


# Fixed result schema. The swept independent variable goes first; `error`
# is always present (NaN/empty on success) so downstream `df["error"]`
# never KeyErrors and column order is identical across runs (#38).
_RESULT_COLS = (
    "knockdown",
    "residual_MPa",
    "pristine_MPa",
    "dpa_mm2",
    "dent_mm",
    "n_delaminations",
    "tier_used",
    "error",
)


# Extra columns of ``sweep_transverse_ratios``, inserted before ``error``.
_BIAXIAL_COLS = (
    "residual_sigma_y_MPa",
    "pristine_sigma_y_MPa",
    "biaxial_interaction_factor",
)


def _finalize(
    rows: List[dict], swept_col: str, result_cols: Sequence[str] = _RESULT_COLS
) -> pd.DataFrame:
    """Build the DataFrame with a deterministic column order regardless of
    how many iterations failed: ``[swept_col, *result_cols]`` with
    ``error`` always present (#38)."""
    return pd.DataFrame(rows).reindex(columns=[swept_col, *result_cols])


def _write_csv(df: pd.DataFrame, csv_path: Optional[Path]) -> None:
    if csv_path is not None:
        df.to_csv(Path(csv_path), index=False)


def _callback_accepts_inputs(cb: Callable[..., Any]) -> bool:
    """Inspect ``cb`` to decide if it accepts the 3-arg ``(i, n, inputs)``
    form. Falls back to the legacy 2-arg form on any introspection failure
    (built-ins, C extensions, partials with unusual signatures, ...)."""
    try:
        params = inspect.signature(cb).parameters
    except (TypeError, ValueError):
        return False
    # Count positional-capable params; *args is treated as accepting all.
    positional = 0
    for p in params.values():
        if p.kind is inspect.Parameter.VAR_POSITIONAL:
            return True
        if p.kind in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        ):
            positional += 1
    return positional >= 3


def _emit_progress(
    cb: Optional[_ProgressCallback],
    i: int,
    n: int,
    inputs: Optional[dict] = None,
) -> None:
    """Invoke the progress callback, tolerating both 2-arg and 3-arg forms.

    A buggy or slow callback (e.g. a Streamlit ``st.progress`` rendering
    against a torn-down session) must never abort the sweep — exceptions
    are swallowed at DEBUG.
    """
    if cb is None:
        return
    try:
        if inputs is not None and _callback_accepts_inputs(cb):
            cb(i, n, inputs)  # type: ignore[call-arg]
        else:
            cb(i, n)  # type: ignore[call-arg]
    except Exception:  # noqa: BLE001 — callback errors must not sink the sweep
        logger.debug("progress_callback raised; ignoring", exc_info=True)


def sweep_energies(
    base_cfg: AnalysisConfig,
    energies_J: Sequence[float],
    csv_path: Optional[Path | str] = None,
    *,
    on_error: OnErrorMode = "raise",
    progress_callback: Optional[_ProgressCallback] = None,
) -> pd.DataFrame:
    """Sweep impact energies; base_cfg must have `impact` set.

    Each point is ``base_cfg`` with ``impact.energy_J`` replaced. Raises
    ``ValueError`` when ``base_cfg.impact`` is ``None``. See the module
    docstring for the columns, ``csv_path``, ``on_error`` and
    ``progress_callback``.

    Examples
    --------
    >>> from bvidfe import AnalysisConfig, ImpactEvent, PanelGeometry
    >>> cfg = AnalysisConfig(material="IM7/8552", layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ...                      ply_thickness_mm=0.152, panel=PanelGeometry(150, 100),
    ...                      impact=ImpactEvent(energy_J=10.0))
    >>> df = sweep_energies(cfg, [5, 10, 20])
    >>> df["knockdown"].round(2).tolist()
    [0.55, 0.45, 0.36]
    """
    _check_on_error(on_error)
    if base_cfg.impact is None:
        raise ValueError("sweep_energies requires base_cfg.impact to be set")
    impact = base_cfg.impact
    rows: List[dict] = []
    n = len(energies_J)
    # Only the impact energy (→ damage) varies here; the mesh-defining inputs
    # are constant, so build_fe_mesh reuses the cached mesh skeleton across
    # iterations and only re-derives the per-element damage factors (#23).
    with _dedupe_impact_warnings():
        for i, E in enumerate(energies_J):
            E = float(E)

            def _build(E: float = E) -> AnalysisConfig:
                return replace(base_cfg, impact=replace(impact, energy_J=E))

            row = _try_run_one(_build, on_error, f"energy_J={E:g}")
            row["energy_J"] = E
            rows.append(row)
            _emit_progress(progress_callback, i + 1, n, {"energy_J": E})
    df = _finalize(rows, "energy_J")
    _write_csv(df, Path(csv_path) if csv_path else None)
    return df


def sweep_layups(
    base_cfg: AnalysisConfig,
    layups: Sequence[Sequence[float]],
    csv_path: Optional[Path | str] = None,
    *,
    on_error: OnErrorMode = "raise",
    progress_callback: Optional[_ProgressCallback] = None,
) -> pd.DataFrame:
    """Sweep layup sequences.

    Each point is ``base_cfg`` with ``layup_deg`` replaced; the ``layup``
    column holds it as a ``"0/45/-45/90"`` string. A per-ply
    ``ply_thickness_mm`` list in ``base_cfg`` must match every layup's
    length (as must a given damage state's interface indices); a layup
    that does not is a failed point under ``on_error``. See the module docstring for the columns, ``csv_path``,
    ``on_error`` and ``progress_callback``.
    """
    _check_on_error(on_error)
    rows: List[dict] = []
    n = len(layups)
    with _dedupe_impact_warnings():
        for i, layup_in in enumerate(layups):
            layup = list(layup_in)
            layup_str = "/".join(f"{a:g}" for a in layup)

            def _build(layup: List[float] = layup) -> AnalysisConfig:
                return replace(base_cfg, layup_deg=layup)

            row = _try_run_one(_build, on_error, f"layup={layup_str}")
            row["layup"] = layup_str
            rows.append(row)
            _emit_progress(progress_callback, i + 1, n, {"layup": layup_str})
    df = _finalize(rows, "layup")
    _write_csv(df, Path(csv_path) if csv_path else None)
    return df


def sweep_thicknesses(
    base_cfg: AnalysisConfig,
    ply_thicknesses_mm: Sequence[float],
    csv_path: Optional[Path | str] = None,
    *,
    on_error: OnErrorMode = "raise",
    progress_callback: Optional[_ProgressCallback] = None,
) -> pd.DataFrame:
    """Sweep ply thickness values.

    Each point is ``base_cfg`` with a uniform ``ply_thickness_mm`` (mm)
    replacing the base value, a per-ply list included. See the module
    docstring for the columns, ``csv_path``, ``on_error`` and
    ``progress_callback``.
    """
    _check_on_error(on_error)
    rows: List[dict] = []
    n = len(ply_thicknesses_mm)
    with _dedupe_impact_warnings():
        for i, t in enumerate(ply_thicknesses_mm):
            t = float(t)

            def _build(t: float = t) -> AnalysisConfig:
                return replace(base_cfg, ply_thickness_mm=t)

            row = _try_run_one(_build, on_error, f"ply_thickness_mm={t:g}")
            row["ply_thickness_mm"] = t
            rows.append(row)
            _emit_progress(progress_callback, i + 1, n, {"ply_thickness_mm": t})
    df = _finalize(rows, "ply_thickness_mm")
    _write_csv(df, Path(csv_path) if csv_path else None)
    return df


def _run_one_ratio(cfg: AnalysisConfig) -> dict:
    """``_run_one``'s columns plus the biaxial ones, from a single run.

    A uniaxial point (``transverse_stress_ratio == 0``) has ``None`` biaxial
    result fields; the sweep writes their values by definition instead
    (sigma_y = s * 0 * lambda = 0 at pristine and damaged failure, and an
    interaction factor of 1), so the ratio-0 point lies on the envelope."""
    result = BvidAnalysis(cfg).run()
    row = _result_row(result)
    uniaxial = cfg.transverse_stress_ratio == 0
    row["residual_sigma_y_MPa"] = 0.0 if uniaxial else result.residual_transverse_stress_MPa
    row["pristine_sigma_y_MPa"] = 0.0 if uniaxial else result.pristine_transverse_stress_MPa
    row["biaxial_interaction_factor"] = 1.0 if uniaxial else result.biaxial_interaction_factor
    return row


def _nan_row_ratio() -> dict:
    """``_nan_row`` plus NaN biaxial columns."""
    return {**_nan_row(), **{c: math.nan for c in _BIAXIAL_COLS}}


def _ratio_point_runner(base_cfg: AnalysisConfig) -> Callable[[AnalysisConfig], dict]:
    """The per-point runner of ``sweep_transverse_ratios``.

    Every point is a full ``BvidAnalysis`` run today, in every tier. This is
    the seam for a tier that can share work across ratios: fe3d's design
    (issue #115) solves two basis load cases (unit sigma_x, unit sigma_y)
    per mesh, damage and loading and superposes them for each ratio, so a
    runner that caches those solves by (mesh signature, damage, loading)
    goes here, returning the same row as ``_run_one_ratio``."""
    return _run_one_ratio


def sweep_transverse_ratios(
    base_cfg: AnalysisConfig,
    ratios: Sequence[float],
    csv_path: Optional[Path | str] = None,
    *,
    on_error: OnErrorMode = "raise",
    progress_callback: Optional[_ProgressCallback] = None,
) -> pd.DataFrame:
    """Sweep the biaxial ratio ``transverse_stress_ratio`` (sigma_y / sigma_x).

    **Experimental, not validated** (issue #115): no biaxial
    residual-strength data exist to check these numbers, so do not read
    the envelope as an allowable. Each point is ``base_cfg`` with
    ``transverse_stress_ratio`` replaced (each in ``[-1, 1]``; a ratio
    outside it is a failed point under ``on_error``), at ``base_cfg``'s
    ``loading`` and ``tier``. Sweep each ``loading`` separately for both
    halves of the envelope (``bvidfe.viz.plots_2d.biaxial_envelope_plot``).

    The columns are ``transverse_stress_ratio``, then the other sweeps'
    (``knockdown``, ``residual_MPa`` and ``pristine_MPa``, which are
    ``|sigma_x|`` at failure along the load ray, ``dpa_mm2``, ``dent_mm``,
    ``n_delaminations``, ``tier_used``), then ``residual_sigma_y_MPa`` and
    ``pristine_sigma_y_MPa`` (the signed sigma_y at damaged and pristine
    failure, ``AnalysisResults.residual_transverse_stress_MPa`` and
    ``pristine_transverse_stress_MPa``), ``biaxial_interaction_factor`` and
    ``error``. A ratio-0 point has sigma_y = 0 and an interaction factor of
    1 by definition (its result fields are ``None``). A tier without a
    biaxial model raises ``NotImplementedError`` at a nonzero ratio: the
    sweep re-raises it with ``on_error="raise"`` and records it in ``error``
    otherwise. See the module docstring for ``csv_path``, ``on_error`` and
    ``progress_callback`` (whose ``inputs`` is
    ``{"transverse_stress_ratio": beta}``).

    Examples
    --------
    >>> from bvidfe import AnalysisConfig, ImpactEvent, PanelGeometry
    >>> cfg = AnalysisConfig(material="IM7/8552", layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ...                      ply_thickness_mm=0.152, panel=PanelGeometry(150, 100),
    ...                      impact=ImpactEvent(energy_J=10.0))
    >>> df = sweep_transverse_ratios(cfg, [0.0, 0.5])
    >>> df["pristine_sigma_y_MPa"].round(1).tolist()
    [0.0, -351.9]
    """
    _check_on_error(on_error)
    run = _ratio_point_runner(base_cfg)
    rows: List[dict] = []
    n = len(ratios)
    with _dedupe_impact_warnings():
        for i, beta_in in enumerate(ratios):
            beta = float(beta_in)

            def _build(beta: float = beta) -> AnalysisConfig:
                return replace(base_cfg, transverse_stress_ratio=beta)

            row = _try_run_one(
                _build,
                on_error,
                f"transverse_stress_ratio={beta:g}",
                run=run,
                nan_row=_nan_row_ratio,
            )
            row["transverse_stress_ratio"] = beta
            rows.append(row)
            _emit_progress(progress_callback, i + 1, n, {"transverse_stress_ratio": beta})
    df = _finalize(rows, "transverse_stress_ratio", (*_RESULT_COLS[:-1], *_BIAXIAL_COLS, "error"))
    _write_csv(df, Path(csv_path) if csv_path else None)
    return df
