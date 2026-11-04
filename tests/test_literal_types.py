"""The ``Literal`` options of :mod:`bvidfe._types` (#82).

Each runtime frozenset is derived from its ``Literal``, and every
construction point rejects a typo'd value with a ``ValueError`` that names
the field, shows the bad value and lists the valid ones.
"""

from __future__ import annotations

import typing

import pytest

from bvidfe import _types
from bvidfe._types import (
    _BOUNDARY_KINDS,
    _CRITERION_NAMES,
    _IMPACTOR_SHAPES,
    _LOADING_MODES,
    _ON_ERROR_MODES,
    _TIER_NAMES,
)
from bvidfe.analysis.config import AnalysisConfig
from bvidfe.analysis.semi_analytical import (
    panel_buckling_load,
    semi_analytical_cai,
    sublaminate_buckling_load,
)
from bvidfe.core import geometry
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure import evaluator
from bvidfe.failure.evaluator import FailureEvaluator
from bvidfe.impact.dent_model import dent_depth_mm
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.sweep.parametric_sweep import sweep_energies, sweep_layups, sweep_thicknesses

_MAT = MATERIAL_LIBRARY["IM7/8552"]
_LAYUP = [0, 45, -45, 90, 90, -45, 45, 0]


def _assert_choice_error(exc: pytest.ExceptionInfo, field: str, bad: object, valid) -> None:
    msg = str(exc.value)
    assert field in msg
    assert repr(bad) in msg
    for v in valid:
        assert repr(v) in msg, (v, msg)


def _cfg(**kw) -> AnalysisConfig:
    base = dict(
        material="IM7/8552",
        layup_deg=_LAYUP,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        impact=ImpactEvent(energy_J=10.0),
    )
    base.update(kw)
    return AnalysisConfig(**base)


# ---------------------------------------------------------------------------
# The frozensets follow the Literals
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "alias,values",
    [
        (_types.TierName, _TIER_NAMES),
        (_types.LoadingMode, _LOADING_MODES),
        (_types.CriterionName, _CRITERION_NAMES),
        (_types.BoundaryKind, _BOUNDARY_KINDS),
        (_types.ImpactorShape, _IMPACTOR_SHAPES),
        (_types.OnErrorMode, _ON_ERROR_MODES),
    ],
)
def test_frozenset_matches_literal(alias, values) -> None:
    assert isinstance(values, frozenset)
    assert values == frozenset(typing.get_args(alias))


def test_geometry_reexports_the_shared_aliases() -> None:
    assert geometry.BoundaryKind is _types.BoundaryKind
    assert geometry.ImpactorShape is _types.ImpactorShape


def test_criterion_registries_cover_exactly_the_literal() -> None:
    assert set(evaluator._CRITERION_REGISTRY) == _CRITERION_NAMES
    assert set(evaluator._CRITERION_SCALAR_REGISTRY) == _CRITERION_NAMES


# ---------------------------------------------------------------------------
# A typo at each construction point
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad", ["fe3D", "semi-analytical", "", None, 3])
def test_analysis_config_tier_typo(bad) -> None:
    with pytest.raises(ValueError) as exc:
        _cfg(tier=bad)
    _assert_choice_error(exc, "tier", bad, _TIER_NAMES)


@pytest.mark.parametrize("bad", ["compresion", "Tension", None])
def test_analysis_config_loading_typo(bad) -> None:
    with pytest.raises(ValueError) as exc:
        _cfg(loading=bad)
    _assert_choice_error(exc, "loading", bad, _LOADING_MODES)


def test_analysis_config_from_dict_tier_typo() -> None:
    d = {
        "material": "IM7/8552",
        "layup_deg": _LAYUP,
        "ply_thickness_mm": 0.152,
        "panel": {"Lx_mm": 150, "Ly_mm": 100},
        "impact": {"energy_J": 10.0},
        "tier": "fe_3d",
    }
    with pytest.raises(ValueError) as exc:
        AnalysisConfig.from_dict(d)
    _assert_choice_error(exc, "tier", "fe_3d", _TIER_NAMES)


@pytest.mark.parametrize("bad", ["simply-supported", "Clamped", "", None])
def test_panel_geometry_boundary_typo(bad) -> None:
    with pytest.raises(ValueError) as exc:
        PanelGeometry(150, 100, boundary=bad)
    _assert_choice_error(exc, "PanelGeometry.boundary", bad, _BOUNDARY_KINDS)


def test_panel_geometry_from_dict_boundary_typo() -> None:
    with pytest.raises(ValueError) as exc:
        PanelGeometry.from_dict({"Lx_mm": 150, "Ly_mm": 100, "boundary": "clampd"})
    _assert_choice_error(exc, "PanelGeometry.boundary", "clampd", _BOUNDARY_KINDS)


@pytest.mark.parametrize("bad", ["hemisphere", "Flat", "cone", None])
def test_impactor_geometry_shape_typo(bad) -> None:
    with pytest.raises(ValueError) as exc:
        ImpactorGeometry(shape=bad)
    _assert_choice_error(exc, "ImpactorGeometry.shape", bad, _IMPACTOR_SHAPES)


def test_impactor_geometry_from_dict_shape_typo() -> None:
    with pytest.raises(ValueError) as exc:
        ImpactorGeometry.from_dict({"shape": "conic"})
    _assert_choice_error(exc, "ImpactorGeometry.shape", "conic", _IMPACTOR_SHAPES)


@pytest.mark.parametrize("bad", ["tsaiwu", "LaRC05", "hashin", None])
def test_failure_evaluator_criterion_typo(bad) -> None:
    with pytest.raises(ValueError) as exc:
        FailureEvaluator(_MAT, criterion=bad)
    _assert_choice_error(exc, "criterion", bad, _CRITERION_NAMES)


@pytest.mark.parametrize(
    "sweep,values",
    [
        (sweep_energies, [10.0]),
        (sweep_layups, [_LAYUP]),
        (sweep_thicknesses, [0.152]),
    ],
)
def test_sweep_on_error_typo(sweep, values) -> None:
    with pytest.raises(ValueError) as exc:
        sweep(_cfg(), values, on_error="ignore")
    _assert_choice_error(exc, "on_error", "ignore", _ON_ERROR_MODES)


# ---------------------------------------------------------------------------
# Public helpers that take the same options no longer fall back silently
# ---------------------------------------------------------------------------


def test_buckling_helpers_reject_boundary_typo() -> None:
    lam = Laminate(_MAT, _LAYUP, 0.152)
    ell = DelaminationEllipse(
        interface_index=0,
        centroid_mm=(75.0, 50.0),
        major_mm=10.0,
        minor_mm=8.0,
        orientation_deg=0.0,
    )
    with pytest.raises(ValueError) as exc:
        panel_buckling_load(lam, 150.0, 100.0, "clampd")  # type: ignore[arg-type]
    _assert_choice_error(exc, "boundary", "clampd", _BOUNDARY_KINDS)
    with pytest.raises(ValueError):
        sublaminate_buckling_load(lam, ell, boundary="SSSS")  # type: ignore[arg-type]
    damage = DamageState(delaminations=[ell], dent_depth_mm=0.2)
    with pytest.raises(ValueError):
        semi_analytical_cai(lam, damage, 1000.0, boundary="fixed")  # type: ignore[arg-type]
    # Every valid value still works.
    for b in _BOUNDARY_KINDS:
        assert panel_buckling_load(lam, 150.0, 100.0, b) > 0  # type: ignore[arg-type]


def test_dent_depth_rejects_shape_typo() -> None:
    with pytest.raises(ValueError) as exc:
        dent_depth_mm(_MAT, 20.0, 5.0, 1.2, shape="cone")  # type: ignore[arg-type]
    _assert_choice_error(exc, "shape", "cone", _IMPACTOR_SHAPES)
    for s in _IMPACTOR_SHAPES:
        assert dent_depth_mm(_MAT, 20.0, 5.0, 1.2, shape=s) >= 0  # type: ignore[arg-type]
