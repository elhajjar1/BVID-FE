"""Streamlit ``app.py`` importability smoke test.

A typo or broken import in ``app.py`` is invisible to the rest of the
test suite — nothing in ``tests/`` imports the file — so the failure
only surfaces when someone visits the deployed Streamlit URL. This
module imports ``app`` from the repository root and fails fast if that
fails.

The test sets ``STREAMLIT_BROWSER_GATHER_USAGE_STATS=false`` so the
import does not attempt any usage-stat network handshake. ``app.py``
runs Streamlit API calls at module top level (``st.set_page_config``,
``st.title``, sidebar widgets); these are tolerated by Streamlit in
"bare mode" (just a noisy warning) and do not raise.

If Streamlit is somehow unavailable in the runtime environment the test
is skipped rather than errored — the package depends on ``streamlit>=1.51``
but defensive skipping keeps this smoke test useful in stripped-down
environments.
"""

from __future__ import annotations

import importlib
import json
import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


streamlit = pytest.importorskip("streamlit", reason="streamlit not installed")


def test_app_module_imports(monkeypatch: pytest.MonkeyPatch) -> None:
    """``import app`` must succeed without raising.

    Catches typos in import statements, syntax errors, and breakage in
    any module ``app.py`` pulls in at import time (``bvidfe.analysis``,
    ``bvidfe.viz``, ``bvidfe.sweep``, ...).
    """
    monkeypatch.setenv("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
    # Telemetry off too — belt-and-braces for offline CI environments.
    monkeypatch.setenv("STREAMLIT_BROWSER_SERVER_ADDRESS", "localhost")

    # Make the repo root importable so ``import app`` resolves to
    # ``<repo>/app.py``. The repo is src-layout, so the root is not
    # automatically on sys.path the way ``src/`` is.
    monkeypatch.syspath_prepend(str(REPO_ROOT))

    # Force a fresh import — otherwise a previously cached module would
    # mask an import-time regression introduced after pytest started.
    sys.modules.pop("app", None)

    app_module = importlib.import_module("app")

    # Spot-check a couple of public names that the Streamlit page relies
    # on; if app.py is imported but its globals never got built, those
    # assertions catch a half-initialised module (e.g. early ``return``).
    assert hasattr(app_module, "MATERIAL_NAMES"), "app.py did not finish initialising"
    assert app_module.MATERIAL_NAMES, "MATERIAL_NAMES should not be empty"
    assert hasattr(app_module, "_parse_layup"), "expected _parse_layup helper in app.py"


def _import_app(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
    monkeypatch.syspath_prepend(str(REPO_ROOT))
    return importlib.import_module("app")


@pytest.mark.parametrize("thickness", [0.152, [0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1]])
def test_config_from_dict_round_trips_ply_thickness(monkeypatch, thickness) -> None:
    """``_config_from_dict`` used to call ``float()`` on ``ply_thickness_mm``
    and so rejected per-ply thickness lists."""
    app_module = _import_app(monkeypatch)
    d = {
        "material": "IM7/8552",
        "layup_deg": [0, 45, -45, 90, 90, -45, 45, 0],
        "ply_thickness_mm": thickness,
        "panel": {"Lx_mm": 150.0, "Ly_mm": 100.0, "boundary": "simply_supported"},
        "loading": "compression",
        "tier": "empirical",
        "impact": {
            "energy_J": 20.0,
            "impactor": {"diameter_mm": 16.0, "shape": "hemispherical"},
            "mass_kg": 5.5,
            "location_xy_mm": None,
        },
        "damage": None,
    }
    cfg = app_module._config_from_dict(d)
    assert cfg.ply_thickness_mm == thickness
    assert cfg.impact.location_xy_mm is None
    assert "mm" in app_module._format_ply_thickness(cfg.ply_thickness_mm)


def test_app_does_not_use_deprecated_use_container_width() -> None:
    """Streamlit deprecated ``use_container_width`` (removal after 2025-12-31);
    ``app.py`` must use ``width="stretch"`` / ``width="content"`` instead.
    A static check, because ``st.button`` / ``st.download_button`` accept the
    old argument without logging a warning, so the run test below misses them."""
    source = (REPO_ROOT / "app.py").read_text(encoding="utf-8")
    hits = [i for i, line in enumerate(source.splitlines(), 1) if "use_container_width" in line]
    assert not hits, f"app.py uses deprecated use_container_width on line(s) {hits}"


class _RecordHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def test_app_run_emits_no_streamlit_deprecation_warnings(monkeypatch) -> None:
    """Run the page headless on its defaults (empirical tier), click
    "Run analysis", and check that nothing raised and that Streamlit logged
    no deprecation warning. Streamlit routes those through the
    ``streamlit.deprecation_util`` logger, which does not propagate."""
    from streamlit.testing.v1 import AppTest

    monkeypatch.setenv("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
    handler = _RecordHandler()
    logger = logging.getLogger("streamlit.deprecation_util")
    logger.addHandler(handler)
    try:
        at = AppTest.from_file(str(REPO_ROOT / "app.py"), default_timeout=120).run()
        assert not at.exception, [e.value for e in at.exception]
        (run_button,) = [b for b in at.button if b.label == "Run analysis"]
        run_button.click().run()
        assert not at.exception, [e.value for e in at.exception]
    finally:
        logger.removeHandler(handler)
    deprecations = [r.getMessage() for r in handler.records]
    assert not deprecations, deprecations
    # The run actually produced results (not just an idle page).
    assert any(m.label == "Knockdown" for m in at.metric)


# ---------------------------------------------------------------------------
# Share-link values and Expert-mode state (#171)
# ---------------------------------------------------------------------------


def _app_test(query: dict | None = None):
    """A fresh headless run of the page, with ``query`` as the URL parameters."""
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(str(REPO_ROOT / "app.py"), default_timeout=120)
    for k, v in (query or {}).items():
        at.query_params[k] = v
    return at.run()


def _no_exception(at) -> None:
    assert not at.exception, [e.value for e in at.exception]


def _run_analysis(at):
    (run_button,) = [b for b in at.button if b.label == "Run analysis"]
    return run_button.click().run()


@pytest.fixture
def _no_stats(monkeypatch):
    monkeypatch.setenv("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")


def test_url_nan_energy_renders_with_warning_not_traceback(_no_stats) -> None:
    """``?E=nan`` used to turn the whole page into a ``ValueError``
    traceback from ``ImpactEvent``; the value is dropped and reported."""
    at = _app_test({"E": "nan"})
    _no_exception(at)
    assert any("E=nan" in w.value and "not finite" in w.value for w in at.warning)
    assert at.number_input(key="energy_J").value == 30.0


@pytest.mark.parametrize(
    "query", [{"E": "500"}, {"E": "-3"}, {"E": "abc"}, {"lx": "5"}, {"lx": "inf"}]
)
def test_url_out_of_range_value_warns_and_keeps_default(_no_stats, query) -> None:
    """An out-of-bounds value used to be silently replaced by the widget's
    minimum (0.1 J, 10 mm); it is now dropped with a warning naming it."""
    at = _app_test(query)
    _no_exception(at)
    ((key, value),) = query.items()
    assert any(f"{key}={value}" in w.value for w in at.warning), [w.value for w in at.warning]
    assert at.number_input(key="energy_J").value == 30.0
    assert at.number_input(key="Lx_mm").value == 150.0


def test_url_location_outside_panel_is_dropped(_no_stats) -> None:
    at = _app_test({"x0": "500", "y0": "30"})
    _no_exception(at)
    assert any("x0=500" in w.value for w in at.warning)
    # The valid y0 still moves the impact off the centre.
    assert at.session_state["loc_centre"] is False
    assert at.session_state["loc_y_mm"] == 30.0


def test_invalid_impact_event_shows_error_not_exception(_no_stats, monkeypatch) -> None:
    """The sidebar's ``ImpactEvent`` construction is guarded: a rejection
    is an ``st.error`` and the Run button reports an incomplete config."""
    import bvidfe.impact.mapping as mapping

    def _reject(*args, **kwargs):
        raise ValueError("ImpactEvent.energy_J must be finite (got nan)")

    monkeypatch.setattr(mapping, "ImpactEvent", _reject)
    at = _app_test()
    _no_exception(at)
    assert any("Invalid impact event" in e.value for e in at.error)
    at = _run_analysis(at)
    _no_exception(at)
    assert any("Configuration is incomplete" in e.value for e in at.error)


def test_state_bounds_match_widget_bounds(monkeypatch) -> None:
    """``_STATE_BOUNDS`` must reject exactly what each widget would."""
    app_module = _import_app(monkeypatch)
    at = _app_test()
    at.toggle(key="expert_mode").set_value(True).run()
    _no_exception(at)
    for key, (lo, hi) in app_module._STATE_BOUNDS.items():
        widget = at.number_input(key=key)
        assert (widget.min, widget.max) == (lo, hi), key


def test_expert_values_survive_expert_mode_toggle(_no_stats) -> None:
    """The Expert-mode widgets had no ``key``: switching Expert mode off and
    on reset them to their literals. They now keep their values, and apply
    while Expert mode is off."""
    at = _app_test()
    at.toggle(key="expert_mode").set_value(True).run()
    at.number_input(key="ply_thickness_mm").set_value(0.188).run()
    at.selectbox(key="boundary").set_value("clamped").run()
    at.selectbox(key="impactor_shape").set_value("flat").run()
    at.number_input(key="mass_kg").set_value(2.0).run()
    at.checkbox(key="loc_centre").uncheck().run()
    at.number_input(key="loc_x_mm").set_value(40.0).run()
    at.number_input(key="loc_y_mm").set_value(30.0).run()
    _no_exception(at)

    at.toggle(key="expert_mode").set_value(False).run()
    _no_exception(at)
    assert any("Expert-mode settings in effect" in c.value for c in at.caption)
    at = _run_analysis(at)
    _no_exception(at)
    cfg = at.session_state["last_result"]["config"]
    assert cfg.ply_thickness_mm == 0.188
    assert cfg.panel.boundary == "clamped"
    assert cfg.impact.impactor.shape == "flat"
    assert cfg.impact.mass_kg == 2.0
    assert cfg.impact.location_xy_mm == (40.0, 30.0)

    at.toggle(key="expert_mode").set_value(True).run()
    _no_exception(at)
    assert at.number_input(key="ply_thickness_mm").value == 0.188
    assert at.selectbox(key="boundary").value == "clamped"
    assert at.selectbox(key="impactor_shape").value == "flat"
    assert at.number_input(key="mass_kg").value == 2.0
    assert at.checkbox(key="loc_centre").value is False
    assert at.number_input(key="loc_x_mm").value == 40.0
    assert at.number_input(key="loc_y_mm").value == 30.0


def test_url_round_trip_includes_expert_keys(_no_stats) -> None:
    """A link carrying the Expert-mode keys hydrates them, and the page
    writes them back so the link it shows reproduces the same analysis."""
    query = {
        "E": "12",
        "tp": "0.188",
        "bc": "clamped",
        "sh": "conical",
        "m": "2.5",
        "x0": "60",
        "y0": "40",
    }
    at = _app_test(query)
    _no_exception(at)
    assert not [w.value for w in at.warning if "shared link" in w.value]
    for k, v in query.items():
        assert at.query_params[k] == v, (k, at.query_params[k])
    assert at.session_state["ply_thickness_mm"] == 0.188
    assert at.session_state["boundary"] == "clamped"
    assert at.session_state["impactor_shape"] == "conical"
    assert at.session_state["mass_kg"] == 2.5
    assert at.session_state["loc_centre"] is False
    assert (at.session_state["loc_x_mm"], at.session_state["loc_y_mm"]) == (60.0, 40.0)


def test_url_centre_impact_omits_location_keys(_no_stats) -> None:
    at = _app_test()
    _no_exception(at)
    for k in ("tp", "bc", "sh", "m"):
        assert k in at.query_params
    assert "x0" not in at.query_params and "y0" not in at.query_params


def test_centre_impact_passes_location_none(_no_stats) -> None:
    """*Centre of panel* runs ``location_xy_mm=None`` (the library's panel
    centre, #171), not explicit centre coordinates."""
    at = _run_analysis(_app_test())
    _no_exception(at)
    payload = at.session_state["last_result"]
    assert payload["config"].impact.location_xy_mm is None
    assert payload["config_dict"]["impact"]["location_xy_mm"] is None


def test_off_centre_corner_point_is_an_error_not_a_centre_run(_no_stats) -> None:
    """An off-centre point at the corner ``(0, 0)`` is not inside the panel.
    The library still reads ``(0, 0)`` as the deprecated centre sentinel, so
    the page must not pass it on: it is an error and the run is blocked."""
    at = _app_test({"x0": "0", "y0": "0"})
    _no_exception(at)
    assert at.session_state["loc_centre"] is False
    assert any("corner" in e.value and "inside the panel" in e.value for e in at.error), [
        e.value for e in at.error
    ]
    at = _run_analysis(at)
    _no_exception(at)
    assert at.session_state["last_result"] is None
    assert any("Configuration is incomplete" in e.value for e in at.error)


def test_presets_reset_expert_values(monkeypatch) -> None:
    app_module = _import_app(monkeypatch)
    for name, preset in app_module.PRESETS.items():
        if name == app_module.PRESET_CUSTOM:
            continue
        for k, v in app_module._EXPERT_DEFAULTS.items():
            assert preset[k] == v, (name, k)


def test_onset_preview_error_is_a_caption(_no_stats) -> None:
    """An impact point on the panel edge cannot be analysed; the preview
    used to vanish silently and now says why."""
    at = _app_test({"x0": "75", "y0": "0"})
    _no_exception(at)
    assert any(
        c.value.startswith("Onset preview unavailable") and "inside the panel" in c.value
        for c in at.caption
    ), [c.value for c in at.caption]


@pytest.mark.parametrize("base_E", [0.1, 1.0, 1.25, 20.0, 200.0])
def test_knockdown_curve_energies_ascend(monkeypatch, base_E) -> None:
    app_module = _import_app(monkeypatch)
    grid = app_module._knockdown_curve_energies(base_E)
    assert len(grid) == 7
    assert all(b > a for a, b in zip(grid, grid[1:]))
    assert grid[0] > 0
    assert grid[0] <= base_E <= grid[-1]


def test_knockdown_curve_tab_low_energy_and_run_tier(_no_stats) -> None:
    """At E = 0.1 J the curve grid used to run 1 J -> 0.2 J (descending),
    and the curve was always empirical. It now defaults to the run's tier."""
    at = _app_test({"E": "0.1", "t": "semi_analytical"})
    at = _run_analysis(at)
    _no_exception(at)
    (curve_tier,) = [s for s in at.selectbox if s.label == "Tier for curve"]
    assert curve_tier.value == "semi_analytical"
    curves = [d.value for d in at.dataframe if "energy_J" in d.value.columns]
    assert curves, "no knockdown-curve table rendered"
    energies = curves[0]["energy_J"].tolist()
    assert len(energies) == 7
    assert all(b > a for a, b in zip(energies, energies[1:]))
    assert set(curves[0]["tier_used"]) == {"semi_analytical"}


@pytest.mark.parametrize("tier", ["empirical", "fe3d"])
def test_summary_shows_the_stiffness(_no_stats, monkeypatch, tier) -> None:
    """Issue #114: the Summary tab shows E11 / E22 / G12, fe3d's E11 as a
    residual with its knockdown and the rest as the pristine CLT value.
    fe3d's solves are replaced by fakes (E11 knockdown 0.8): the app also
    runs a 7-point knockdown curve at the run's tier."""
    from bvidfe.analysis import bvid

    monkeypatch.setattr(
        bvid,
        "fe3d_cai_channels_and_modulus",
        lambda cfg, damage, lam, cap: (None, 100.0, 0.8 if damage.delaminations else 1.0),
    )
    at = _run_analysis(_app_test({"t": tier, "ld": "compression"}))
    _no_exception(at)
    labels = {m.label for m in at.metric}
    result_dict = at.session_state["last_result"]["result_dict"]
    if tier == "fe3d":
        assert result_dict["E11_knockdown"] == pytest.approx(0.8)
        assert "E11 residual" in labels
        assert {"E22 (pristine, CLT)", "G12 (pristine, CLT)"} <= labels
        (e11,) = [m for m in at.metric if m.label == "E11 residual"]
        assert e11.value == f"{result_dict['residual_E11_MPa'] / 1000:.1f} GPa"
    else:
        assert {"E11 (pristine, CLT)", "E22 (pristine, CLT)", "G12 (pristine, CLT)"} <= labels


# ---------------------------------------------------------------------------
# Layup range (#169) and the config / results JSON (#170)
# ---------------------------------------------------------------------------


def test_parse_layup_accepts_any_finite_angle_like_the_cli(monkeypatch) -> None:
    """The app rejected angles outside [-90, 90] while the CLI and
    ``AnalysisConfig`` take any finite angle (used modulo 180 deg)."""
    import bvidfe.cli as cli

    app_module = _import_app(monkeypatch)
    text = "0, 135, -135, 270"
    assert app_module._parse_layup(text) == [0.0, 135.0, -135.0, 270.0]
    assert app_module._parse_layup(text) == cli._parse_layup(text.replace(" ", ""))
    for bad in ("0, nan", "0, inf", "0, -inf"):
        with pytest.raises(ValueError, match="finite"):
            app_module._parse_layup(bad)
    with pytest.raises(ValueError, match="cap"):
        app_module._parse_layup(", ".join(["0"] * (app_module.MAX_PLIES + 1)))


def test_layup_outside_90_runs_in_the_app(_no_stats) -> None:
    layup = [0.0, 135.0, -135.0, 90.0, 90.0, -135.0, 135.0, 0.0]
    at = _app_test({"lay": ",".join(f"{a:g}" for a in layup)})
    _no_exception(at)
    assert not [e.value for e in at.error if "layup" in e.value.lower()]
    at = _run_analysis(at)
    _no_exception(at)
    assert at.session_state["last_result"]["config"].layup_deg == layup


_QI8 = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]


def test_config_dict_damage_block_is_the_cscan_layout(monkeypatch) -> None:
    """The app's config dict (the run's cache key and the source of its
    config) wrote the damage block by hand, without ``schema_version``; it
    is ``damage_state_to_dict``, which ``damage_state_from_dict`` reads
    back, and the whole dict loads with ``AnalysisConfig.from_dict``."""
    from bvidfe.analysis import AnalysisConfig
    from bvidfe.core.geometry import PanelGeometry
    from bvidfe.damage.io import SCHEMA_VERSION, damage_state_from_dict
    from bvidfe.damage.state import DamageState, DelaminationEllipse

    app_module = _import_app(monkeypatch)
    damage = DamageState(
        [
            DelaminationEllipse(1, (75.0, 50.0), 20.0, 12.0, 45.0),
            DelaminationEllipse(3, (70.0, 48.0), 25.0, 10.0, -45.0),
        ],
        dent_depth_mm=0.3,
        fiber_break_radius_mm=1.5,
    )
    panel = PanelGeometry(150.0, 100.0)
    d = app_module._config_dict(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=panel,
        loading="compression",
        tier="empirical",
        impact_event=None,
        damage_state=damage,
        mesh_params=None,
    )
    d = json.loads(json.dumps(d))
    assert d["damage"]["schema_version"] == SCHEMA_VERSION
    assert damage_state_from_dict(d["damage"]) == damage
    cfg = AnalysisConfig.from_dict(d)
    assert (cfg.damage, cfg.panel, cfg.impact) == (damage, panel, None)
    assert app_module._config_from_dict(d) == cfg


def test_config_dict_impact_block_loads(monkeypatch) -> None:
    from bvidfe.analysis import AnalysisConfig, MeshParams
    from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
    from bvidfe.impact.mapping import ImpactEvent

    app_module = _import_app(monkeypatch)
    impact = ImpactEvent(20.0, ImpactorGeometry(12.7, "conical"), 2.0, (60.0, 40.0))
    kw: dict = dict(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=[0.152] * 8,
        panel=PanelGeometry(150.0, 100.0, "clamped"),
        loading="tension",
        tier="fe3d",
        impact_event=impact,
        damage_state=None,
        mesh_params=MeshParams(2, 2.5),
    )
    cfg = AnalysisConfig.from_dict(json.loads(json.dumps(app_module._config_dict(**kw))))
    assert cfg == AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=[0.152] * 8,
        panel=kw["panel"],
        loading="tension",
        tier="fe3d",
        impact=impact,
        mesh=MeshParams(2, 2.5),
    )
    centre = app_module._config_dict(**{**kw, "impact_event": ImpactEvent(20.0)})
    assert centre["impact"]["location_xy_mm"] is None
    assert AnalysisConfig.from_dict(json.loads(json.dumps(centre))).impact == ImpactEvent(20.0)


def test_downloaded_results_json_reloads(_no_stats) -> None:
    """The Summary tab's *Download results JSON* (``json.dumps`` of the
    run's ``result_dict``): its ``damage`` reads back with
    ``damage_state_from_dict`` and its ``config_snapshot`` with
    ``AnalysisConfig.from_dict``, giving the analysed damage and config."""
    from bvidfe.analysis import AnalysisConfig
    from bvidfe.damage.io import damage_state_from_dict

    at = _run_analysis(_app_test())
    _no_exception(at)
    payload = at.session_state["last_result"]
    downloaded = json.loads(json.dumps(payload["result_dict"], indent=2))
    assert damage_state_from_dict(downloaded["damage"]) == payload["damage"]
    assert AnalysisConfig.from_dict(downloaded["config_snapshot"]) == payload["config"]
    assert payload["config"] == AnalysisConfig.from_dict(
        json.loads(json.dumps(payload["config_dict"]))
    )


def test_unused_ply_thickness_helper_is_gone(monkeypatch) -> None:
    """``_ply_thickness_from_json`` was unused (``AnalysisConfig.from_dict``
    reads a uniform or per-ply thickness)."""
    assert not hasattr(_import_app(monkeypatch), "_ply_thickness_from_json")


# ---------------------------------------------------------------------------
# fe3d stage log in an st.status block (#110)
# ---------------------------------------------------------------------------

# A small, fast fe3d run: 768 elements, ~3k DOFs at the 5 mm default mesh.
_SMALL_FE3D = {"t": "fe3d", "lx": "60", "ly": "40", "E": "2"}


def _stage_handlers() -> list:
    return [
        h
        for h in logging.getLogger("bvidfe.fe3d").handlers
        if type(h).__name__ == "_StageStatusHandler"
    ]


def test_fe3d_run_shows_its_stage_log_in_a_status_block(_no_stats) -> None:
    streamlit.cache_data.clear()
    at = _run_analysis(_app_test(_SMALL_FE3D))
    _no_exception(at)
    assert not at.error, [e.value for e in at.error]
    assert any(m.label == "Knockdown" for m in at.metric)
    (status,) = at.status
    assert status.label.startswith("fe3d stage log (") and status.label.endswith(" lines)")
    (log,) = status.text
    for line in [
        "fe3d FPF (damaged, ",
        "fe3d FPF (undamaged): mesh build start",
        "element build done",
        "assembly done",
        "BCs done",
        "solve done: direct LU",
        "FPF evaluation done",
        "fe3d buckling",
    ]:
        assert line in log.value, line
    assert not _stage_handlers()

    # Same config again: a cache hit, nothing solved and nothing logged.
    at = _run_analysis(at)
    _no_exception(at)
    (status,) = at.status
    assert "cached result" in status.label
    assert not _stage_handlers()


def test_non_fe3d_run_has_no_status_block(_no_stats) -> None:
    at = _run_analysis(_app_test({"t": "semi_analytical"}))
    _no_exception(at)
    assert not at.status


def test_fe3d_stage_handler_is_removed_when_the_run_fails(_no_stats, monkeypatch) -> None:
    import bvidfe.analysis.bvid as bvid

    def _fail(*args, **kwargs):
        logging.getLogger("bvidfe.fe3d").info("about to fail")
        raise RuntimeError("solver exploded")

    monkeypatch.setattr(bvid, "fe3d_cai_channels_and_modulus", _fail)
    streamlit.cache_data.clear()
    at = _run_analysis(_app_test({**_SMALL_FE3D, "E": "1.5"}))
    _no_exception(at)
    assert any("solver exploded" in e.value for e in at.error)
    (status,) = at.status
    assert status.proto.state == status.proto.State.ERROR
    assert "about to fail" in status.text[0].value
    assert not _stage_handlers()


def test_stage_handler_updates_the_label_live_for_its_own_thread_only(monkeypatch) -> None:
    """Each INFO line becomes the status label as it is logged; records
    from other threads (other sessions' runs) are ignored."""
    import threading

    app_module = _import_app(monkeypatch)

    class _Status:
        def __init__(self):
            self.labels: list[str] = []

        def update(self, *, label):
            self.labels.append(label)

    status = _Status()
    handler = app_module._StageStatusHandler(status)
    log = logging.getLogger("bvidfe.fe3d.test_stage_handler")
    monkeypatch.setattr(log, "propagate", False)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)
    try:
        log.info("mesh build start")
        log.debug("not a stage line")
        other = threading.Thread(target=log.info, args=("another session",))
        other.start()
        other.join()
        log.info("solve done: direct LU (0.05s)")
    finally:
        log.removeHandler(handler)
        log.setLevel(logging.NOTSET)
    assert status.labels == ["fe3d: mesh build start", "fe3d: solve done: direct LU (0.05s)"]
    assert [line.split(" ", 1)[1] for line in handler.lines] == [
        "mesh build start",
        "solve done: direct LU (0.05s)",
    ]


# ---------------------------------------------------------------------------
# Biaxial loading: the sigma_y / sigma_x input, URL key ``tr``, Summary,
# ratio sweep and NCR (issue #115)
# ---------------------------------------------------------------------------


def _fake_unsupported_tier(monkeypatch, tier: str) -> None:
    """Make ``tier`` raise ``NotImplementedError`` at a nonzero ratio, as a
    tier without a biaxial model does, whichever tiers support it today."""
    from bvidfe.analysis import BvidAnalysis

    real = BvidAnalysis.run

    def _run(self):
        if self.config.tier == tier and self.config.transverse_stress_ratio != 0:
            raise NotImplementedError(f"transverse_stress_ratio != 0 is not supported by {tier}")
        return real(self)

    monkeypatch.setattr(BvidAnalysis, "run", _run)
    streamlit.cache_data.clear()


def _sweep_tab_select_ratio(at):
    at.radio(key="sweep_kind").set_value("σy/σx ratio (experimental)").run()
    _no_exception(at)
    return at


def _run_ratio_sweep(at):
    (button,) = [b for b in at.button if b.label == "Run ratio sweep"]
    return button.click().run()


def test_default_ratio_is_zero_and_the_app_is_unchanged(_no_stats) -> None:
    """With the default 0: no ``tr`` in the URL, a uniaxial run, no biaxial
    metrics or caveat, and an NCR without a biaxial line."""
    at = _app_test()
    _no_exception(at)
    widget = at.number_input(key="transverse_stress_ratio")
    assert widget.label == "σy/σx (experimental)"
    assert (widget.value, widget.min, widget.max) == (0.0, -1.0, 1.0)
    assert "tr" not in at.query_params
    at = _run_analysis(at)
    _no_exception(at)
    payload = at.session_state["last_result"]
    assert payload["config"].transverse_stress_ratio == 0.0
    assert payload["result_dict"]["biaxial_interaction_factor"] is None
    labels = {m.label for m in at.metric}
    assert "Biaxial interaction factor" not in labels
    assert not [w.value for w in at.warning if "Biaxial" in w.value]
    ncr_md = [md.value for md in at.markdown if "Nonconformance Report" in md.value]
    assert ncr_md and "Biaxial" not in ncr_md[0]


def test_url_tr_runs_a_biaxial_analysis_and_round_trips(_no_stats) -> None:
    from bvidfe.analysis import AnalysisConfig, BvidAnalysis

    at = _app_test({"tr": "0.5", "ld": "tension"})
    _no_exception(at)
    assert not [w.value for w in at.warning if "shared link" in w.value]
    assert at.number_input(key="transverse_stress_ratio").value == 0.5
    assert at.query_params["tr"] == "0.5"
    assert any("experimental, not validated" in c.value for c in at.caption)
    at = _run_analysis(at)
    _no_exception(at)
    payload = at.session_state["last_result"]
    cfg = payload["config"]
    assert (cfg.transverse_stress_ratio, cfg.loading) == (0.5, "tension")
    assert payload["config_dict"]["transverse_stress_ratio"] == 0.5
    expected = BvidAnalysis(AnalysisConfig.from_dict(payload["config_dict"])).run()
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["σy at pristine failure"] == (
        f"{expected.pristine_transverse_stress_MPa:.1f} MPa"
    )
    assert metrics["σy at residual failure"] == (
        f"{expected.residual_transverse_stress_MPa:.1f} MPa"
    )
    assert metrics["Biaxial interaction factor"] == f"{expected.biaxial_interaction_factor:.3f}"
    assert any("σy/σx = +0.50" in w.value and "not validated" in w.value for w in at.warning)
    # The NCR export records the ratio and flags it.
    ncr_md = [md.value for md in at.markdown if "Nonconformance Report" in md.value]
    assert ncr_md and "σy/σx = +0.50 (tension)" in ncr_md[0]
    assert "NOT VALIDATED" in ncr_md[0]


@pytest.mark.parametrize("value", ["1.5", "-2", "nan", "abc"])
def test_url_tr_out_of_range_is_dropped(_no_stats, value) -> None:
    at = _app_test({"tr": value})
    _no_exception(at)
    assert any(f"tr={value}" in w.value for w in at.warning), [w.value for w in at.warning]
    assert at.number_input(key="transverse_stress_ratio").value == 0.0
    assert "tr" not in at.query_params


def test_preset_resets_the_ratio(_no_stats) -> None:
    at = _app_test({"tr": "-0.5"})
    at.selectbox(key="preset_choice").set_value("Empirical quick (IM7/8552 quasi-iso, 20 J)").run()
    _no_exception(at)
    assert at.number_input(key="transverse_stress_ratio").value == 0.0


def test_unsupported_tier_shows_a_message_not_a_traceback(_no_stats, monkeypatch) -> None:
    _fake_unsupported_tier(monkeypatch, "semi_analytical")
    at = _run_analysis(_app_test({"tr": "0.5", "t": "semi_analytical"}))
    _no_exception(at)
    assert at.session_state["last_result"] is None
    (err,) = [e.value for e in at.error if "cannot run this configuration" in e.value]
    assert "semi_analytical" in err and "set σy/σx to 0" in err


def test_ratio_sweep_tab_draws_the_envelope(_no_stats) -> None:
    at = _sweep_tab_select_ratio(_app_test())
    assert any("not validated" in w.value for w in at.warning)
    at.text_input(key="ratio_sweep_text").set_value("-1, 0, 0.5").run()
    at = _run_ratio_sweep(at)
    _no_exception(at)
    assert not at.error, [e.value for e in at.error]
    tables = [d.value for d in at.dataframe if "transverse_stress_ratio" in d.value.columns]
    assert tables, "no ratio-sweep table rendered"
    df = tables[0]
    assert df["loading"].tolist() == ["compression"] * 3 + ["tension"] * 3
    assert df["transverse_stress_ratio"].tolist() == [-1.0, 0.0, 0.5] * 2
    assert df["error"].isna().all()
    assert {"residual_sigma_y_MPa", "pristine_sigma_y_MPa", "biaxial_interaction_factor"} <= set(
        df.columns
    )


def test_ratio_sweep_rejects_bad_ratios(_no_stats) -> None:
    at = _sweep_tab_select_ratio(_app_test())
    at.text_input(key="ratio_sweep_text").set_value("0, 1.5").run()
    at = _run_ratio_sweep(at)
    _no_exception(at)
    assert any("Invalid ratios" in e.value and "between -1 and 1" in e.value for e in at.error)


def test_ratio_sweep_on_an_unsupported_tier_says_so(_no_stats, monkeypatch) -> None:
    _fake_unsupported_tier(monkeypatch, "semi_analytical")
    at = _sweep_tab_select_ratio(_app_test())
    at.selectbox(key="ratio_sweep_tier").set_value("semi_analytical").run()
    at.checkbox(key="ratio_sweep_both").uncheck().run()
    at = _run_ratio_sweep(at)
    _no_exception(at)
    assert any(
        "semi_analytical" in w.value and "cannot run 4 of these point(s)" in w.value
        for w in at.warning
    ), [w.value for w in at.warning]
    tables = [d.value for d in at.dataframe if "transverse_stress_ratio" in d.value.columns]
    df = tables[0]
    assert df.loc[df["transverse_stress_ratio"] == 0.0, "error"].isna().all()
    assert df.loc[df["transverse_stress_ratio"] != 0.0, "error"].notna().all()


def test_ratio_sweep_does_not_need_impact_mode(_no_stats) -> None:
    """The energy sweep needs Impact event mode, the ratio sweep only a
    complete config: with C-scan mode and no file it asks for the inputs."""
    at = _app_test()
    at.radio(key="input_mode").set_value("C-scan inspection").run()
    _no_exception(at)
    at = _sweep_tab_select_ratio(at)
    assert any("Fix the sidebar inputs" in i.value for i in at.info)
    assert not any("require **Impact event** mode" in i.value for i in at.info)
