"""Analysis results dataclasses."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from bvidfe.damage.state import DamageState


@dataclass
class FieldResults:
    """Per-mesh field outputs reserved for the 3D FE tier.

    .. note::
       **Currently always ``None``.** ``BvidAnalysis.run()`` does not yet
       populate ``AnalysisResults.field_results`` for any tier (the fe3d
       tier returns scalar residual strength + the scalar buckling stress
       stored in ``buckling_eigenvalues``).
       This dataclass documents the *intended* contract for a future
       field-output release; downstream code should treat
       ``result.field_results`` as optional and handle ``None``.

    Attributes
    ----------
    displacement : np.ndarray
        Nodal displacement, shape ``(n_nodes, 3)``, float64, units mm,
        in the global (x, y, z) coordinate system, original (undeformed)
        node ordering.
    stress_global : np.ndarray
        Per-element stress in the global frame, shape
        ``(n_elements, 6)``, units MPa, Voigt order
        ``[σxx, σyy, σzz, σyz, σxz, σxy]``.
    strain_global : np.ndarray
        Per-element strain in the global frame, shape
        ``(n_elements, 6)``, dimensionless, same Voigt order as
        ``stress_global``.
    stress_local : np.ndarray
        Per-element stress rotated into each ply's material (fibre) axes,
        shape ``(n_elements, 6)``, units MPa, Voigt order
        ``[σ11, σ22, σ33, σ23, σ13, σ12]``.
    failure_index : np.ndarray
        Per-element maximum failure index, shape ``(n_elements,)``,
        dimensionless (≥ 1.0 indicates predicted first-ply failure).
    buckling_mode_shape : np.ndarray, optional
        First buckling eigenvector as nodal displacement, shape
        ``(n_nodes, 3)``, normalised (dimensionless). ``None`` when no
        buckling solve ran or it did not converge.
    cohesive_damage : np.ndarray, optional
        Per-interface-element scalar cohesive damage variable in
        ``[0, 1]`` (0 = intact, 1 = fully separated). ``None`` until
        true cohesive surfaces land (see README "Limitations").
    """

    displacement: np.ndarray
    stress_global: np.ndarray
    strain_global: np.ndarray
    stress_local: np.ndarray
    failure_index: np.ndarray
    buckling_mode_shape: Optional[np.ndarray] = None
    cohesive_damage: Optional[np.ndarray] = None


@dataclass
class AnalysisResults:
    """Outcome of a single ``BvidAnalysis.run()`` invocation.

    The strength fields are defined identically across all three tiers:

    - ``pristine_strength_MPa`` is the unnotched laminate strength:
      ply-discount last-ply failure on classical lamination theory (see
      ``bvidfe.failure.laminate_strength``). Same denominator for every
      tier.
    - ``residual_strength_MPa`` is the tier-specific damaged strength. For
      fe3d it is ``pristine_strength_MPa`` times fe3d's damaged strength
      over fe3d's own undamaged strength, so an undamaged fe3d run gives
      a knockdown of 1.0; the raw fe3d stresses are recorded in ``notes``.
    - ``knockdown = residual_strength_MPa / pristine_strength_MPa``,
      assigned in ``BvidAnalysis.run()``.

    Knockdowns from different tiers are on the same scale (same baseline),
    but capture different failure mechanisms — see the README "Knockdown
    definition and cross-tier comparability" section before comparing
    values across tiers.

    ``notes`` carries free-form runtime diagnostics emitted by the analysis
    backends that affect the interpretation of ``knockdown``. fe3d runs
    always add one line with the raw damaged and undamaged fe3d stresses
    behind the knockdown; fe3d compression runs also record the governing
    delamination growth stress and the whole-panel buckling stress that is
    left out of the knockdown. Otherwise empty when the run produced no
    diagnostic-worthy events.

    ``buckling_eigenvalues`` is buckling *onset*, not strength: the
    sublaminate buckling load (N/mm) at ``critical_sublaminate`` for
    semi_analytical, and the lower of whole-panel and sublaminate buckling
    stress (MPa) for fe3d. ``critical_sublaminate`` is the interface whose
    delamination grows at the lowest stress (CAI only).

    Attributes
    ----------
    residual_strength_MPa : float
        Damaged strength (MPa) from the selected tier.
    pristine_strength_MPa : float
        Unnotched laminate strength (MPa), identical for every tier.
    knockdown : float
        ``residual_strength_MPa / pristine_strength_MPa`` (``0.0`` if the
        pristine strength is not positive).
    damage : DamageState
        The damage state analysed: the given one, or the one predicted from
        the impact (placed on ``panel``).
    dpa_mm2 : float
        ``damage.projected_damage_area_mm2``, the plan-view union area (mm²).
    tier_used : str
        The tier that ran.
    config_snapshot : dict
        ``dataclasses.asdict`` of the ``AnalysisConfig`` (a custom material
        expanded to a dict), for provenance.
    buckling_eigenvalues : list of float, optional
        One value: N/mm for semi_analytical CAI, MPa for fe3d CAI. ``None``
        for the empirical tier, for tension, and when there is no
        delamination (semi_analytical) or no positive buckling stress (fe3d).
    critical_sublaminate : int, optional
        Interface index from
        :func:`~bvidfe.analysis.semi_analytical.weakest_sublaminate_growth`
        for semi_analytical and fe3d CAI; ``None`` otherwise.
    field_results : FieldResults, optional
        Always ``None`` at present.
    notes : list of str
        Human-readable diagnostics (see above).
    warnings : list of str
        Machine-readable tags; always empty at present (see below).
    """

    residual_strength_MPa: float
    pristine_strength_MPa: float
    knockdown: float
    damage: DamageState
    dpa_mm2: float
    tier_used: str
    config_snapshot: Dict[str, Any]
    buckling_eigenvalues: Optional[List[float]] = None
    critical_sublaminate: Optional[int] = None
    field_results: Optional[FieldResults] = None
    notes: List[str] = field(default_factory=list)
    #: Every ``UserWarning`` raised during ``BvidAnalysis.run`` (e.g.
    #: ``DPACapClipWarning``, ``SmallMassQuasiStaticWarning``,
    #: ``DamageOutsidePanelWarning``), as ``"<WarningClass>: <message>"``,
    #: de-duplicated, in the order raised. The warnings are still emitted
    #: through Python's ``warnings`` machinery as well.
    warnings: List[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"AnalysisResults(tier={self.tier_used!r}, "
            f"kd={self.knockdown:.3f}, "
            f"sigma={self.residual_strength_MPa:.1f}/"
            f"{self.pristine_strength_MPa:.1f}MPa, "
            f"dpa={self.dpa_mm2:.1f}mm2, "
            f"notes={len(self.notes)}, warnings={len(self.warnings)})"
        )

    def summary(self) -> str:
        """Multi-line, human-readable report: tier, strengths, knockdown,
        damage area, dent depth, number of delaminations, then any
        ``notes`` and ``warnings``."""
        lines = [
            "BVID Analysis Results",
            f"  tier_used              : {self.tier_used}",
            f"  pristine_strength_MPa  : {self.pristine_strength_MPa:.1f}",
            f"  residual_strength_MPa  : {self.residual_strength_MPa:.1f}",
            f"  knockdown              : {self.knockdown:.3f}",
            f"  dpa_mm2                : {self.dpa_mm2:.1f}",
            f"  dent_depth_mm          : {self.damage.dent_depth_mm:.3f}",
            f"  n_delaminations        : {len(self.damage.delaminations)}",
        ]
        if self.notes:
            lines.append("  notes:")
            for note in self.notes:
                lines.append(f"    - {note}")
        if self.warnings:
            lines.append("  warnings:")
            for w in self.warnings:
                lines.append(f"    - {w}")
        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serialisable dict of the result (the CLI's JSON output).

        Keys: ``residual_strength_MPa``, ``pristine_strength_MPa``,
        ``knockdown``, ``dpa_mm2``, ``tier_used``, ``damage`` (in the
        C-scan JSON layout without ``schema_version``),
        ``buckling_eigenvalues``, ``critical_sublaminate``,
        ``config_snapshot``, ``notes`` and ``warnings``. ``field_results``
        is left out. ``config_snapshot`` holds the config's dataclasses as
        dicts, with tuples kept as tuples (``json.dumps`` writes them as
        lists).
        """
        out = {
            "residual_strength_MPa": self.residual_strength_MPa,
            "pristine_strength_MPa": self.pristine_strength_MPa,
            "knockdown": self.knockdown,
            "dpa_mm2": self.dpa_mm2,
            "tier_used": self.tier_used,
            "damage": {
                "dent_depth_mm": self.damage.dent_depth_mm,
                "fiber_break_radius_mm": self.damage.fiber_break_radius_mm,
                "delaminations": [
                    {
                        "interface_index": d.interface_index,
                        "centroid_mm": list(d.centroid_mm),
                        "major_mm": d.major_mm,
                        "minor_mm": d.minor_mm,
                        "orientation_deg": d.orientation_deg,
                    }
                    for d in self.damage.delaminations
                ],
            },
            "buckling_eigenvalues": self.buckling_eigenvalues,
            "critical_sublaminate": self.critical_sublaminate,
            "config_snapshot": self.config_snapshot,
            "notes": list(self.notes),
            "warnings": list(self.warnings),
        }
        return out
