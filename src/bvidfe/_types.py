"""Shared Literal type aliases and their runtime-validation sets.

This module is the single source of truth for the string-enum-style options
that flow through the public API (``AnalysisConfig``, ``PanelGeometry``,
``ImpactorGeometry``, ``FailureEvaluator``, the sweep entry points and the
CLI). Each option is a ``typing.Literal`` alias, so mypy / Pyright and IDE
completion see the valid values in the signatures, and each has a matching
frozenset derived from the Literal with ``typing.get_args``, so the runtime
check cannot drift from the static type. Adding a value means editing the
Literal only.

The values stay plain strings (not ``Enum`` members) so JSON configs, the
CLI and ``AnalysisResults.to_dict()`` keep the same spelling.

Construction-time validation (each raises ``ValueError`` listing the valid
values):

- ``AnalysisConfig.tier`` / ``.loading`` against ``_TIER_NAMES`` /
  ``_LOADING_MODES`` (``bvidfe.analysis.config``);
- ``PanelGeometry.boundary`` / ``ImpactorGeometry.shape`` against
  ``_BOUNDARY_KINDS`` / ``_IMPACTOR_SHAPES`` (``bvidfe.core.geometry``);
- ``FailureEvaluator(criterion=...)`` and the fe3d first-ply-failure solve
  against ``_CRITERION_NAMES`` (``bvidfe.failure.evaluator``,
  ``bvidfe.analysis.fe_tier``);
- the sweeps' ``on_error`` against ``_ON_ERROR_MODES``
  (``bvidfe.sweep.parametric_sweep``).
"""

from __future__ import annotations

from typing import Literal, get_args

# ---------------------------------------------------------------------------
# Literal aliases
# ---------------------------------------------------------------------------

#: Identifier for the analysis tier selected on ``AnalysisConfig.tier``.
TierName = Literal["empirical", "semi_analytical", "fe3d"]

#: In-plane loading mode selected on ``AnalysisConfig.loading``.
LoadingMode = Literal["compression", "tension"]

#: Failure-criterion name. Extend this Literal in lockstep with the registry
#: in :mod:`bvidfe.failure.evaluator` (and the analytic-scaling branches of
#: :mod:`bvidfe.analysis.fe_tier`) when adding a new criterion.
CriterionName = Literal["tsai_wu", "larc05", "puck"]

#: Panel edge condition on ``PanelGeometry.boundary``. Re-exported by
#: :mod:`bvidfe.core.geometry`.
BoundaryKind = Literal["clamped", "simply_supported", "free"]

#: Impactor nose shape on ``ImpactorGeometry.shape``. Re-exported by
#: :mod:`bvidfe.core.geometry`.
ImpactorShape = Literal["hemispherical", "flat", "conical"]

#: What a parametric sweep does with a point that raises: re-raise it, or
#: record a NaN row (``"skip"`` silently, ``"warn"`` with a ``UserWarning``).
OnErrorMode = Literal["raise", "skip", "warn"]


# ---------------------------------------------------------------------------
# Runtime-validation sets, derived from the Literals so they cannot drift
# ---------------------------------------------------------------------------

_TIER_NAMES: frozenset[str] = frozenset(get_args(TierName))
_LOADING_MODES: frozenset[str] = frozenset(get_args(LoadingMode))
_CRITERION_NAMES: frozenset[str] = frozenset(get_args(CriterionName))
_BOUNDARY_KINDS: frozenset[str] = frozenset(get_args(BoundaryKind))
_IMPACTOR_SHAPES: frozenset[str] = frozenset(get_args(ImpactorShape))
_ON_ERROR_MODES: frozenset[str] = frozenset(get_args(OnErrorMode))


__all__ = [
    "TierName",
    "LoadingMode",
    "CriterionName",
    "BoundaryKind",
    "ImpactorShape",
    "OnErrorMode",
    "_TIER_NAMES",
    "_LOADING_MODES",
    "_CRITERION_NAMES",
    "_BOUNDARY_KINDS",
    "_IMPACTOR_SHAPES",
    "_ON_ERROR_MODES",
]
