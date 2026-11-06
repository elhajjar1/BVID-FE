"""Analysis results dataclasses."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from bvidfe._types import TierName
from bvidfe.damage.io import damage_state_to_dict
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
        true cohesive surfaces land (see "Limitations" in docs/physics_models.md).
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
      fe3d it is ``pristine_strength_MPa`` times fe3d's damaged over its
      own undamaged first-ply-failure stress, so an undamaged fe3d run
      gives a knockdown of 1.0; in compression the delamination growth
      stress instead where it is lower (the closed-form stress
      semi_analytical uses, not rescaled). The raw fe3d stresses are
      recorded in ``notes``.
    - ``knockdown = residual_strength_MPa / pristine_strength_MPa``,
      assigned in ``BvidAnalysis.run()``.

    Knockdowns from different tiers are on the same scale (same baseline),
    but capture different failure mechanisms — see "Knockdown definition and
    cross-tier comparability" in docs/physics_models.md before comparing
    values across tiers.

    ``notes`` carries free-form runtime diagnostics emitted by the analysis
    backends that affect the interpretation of ``knockdown``. fe3d runs
    always add one line with the raw damaged and undamaged fe3d stresses
    behind the knockdown; fe3d compression runs also record the weakest
    delamination's growth stress, which of it and first-ply failure
    governs, and the whole-panel buckling stress that is left out of the
    knockdown. Otherwise empty when the run produced no
    diagnostic-worthy events.

    The stiffness fields (issue #114) are the laminate's effective membrane
    moduli in the panel axes: ``E11`` along x (the loading direction and
    the layup's 0 deg direction), ``E22`` along y and ``G12`` the in-plane
    shear modulus, i.e. the CLT ``Ex``, ``Ey`` and ``Gxy`` of
    :meth:`~bvidfe.core.laminate.Laminate.effective_engineering_constants`
    (not the ply's ``E11`` / ``E22`` / ``G12``). ``pristine_*`` is that CLT
    value for every tier. Only fe3d computes a residual, and only
    ``residual_E11_MPa``: ``pristine_E11_MPa`` times fe3d's damaged over
    undamaged secant axial modulus of the whole panel (the loaded-edge
    reaction forces of the two first-ply failure solves, at the same end
    strain), so in the run's loading sense (compression for CAI, tension
    for TAI). The closed-form tiers model the damage as a strength
    equivalent (a hole, a notch, buckled sublaminates) with no stiffness,
    so their residual fields are ``None``, as are fe3d's E22 and G12
    (they would need further solves). See "Residual stiffness" in
    docs/physics_models.md; no residual stiffness measurement validates
    them.

    ``buckling_eigenvalues`` is buckling *onset*, not strength: the
    sublaminate buckling load (N/mm) at ``critical_sublaminate`` for
    semi_analytical, and the lower of whole-panel and sublaminate buckling
    stress (MPa) for fe3d. ``critical_sublaminate`` is the interface whose
    delamination grows at the lowest stress (CAI, and semi_analytical TAI
    under a lateral compression, ``transverse_stress_ratio < 0``).

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
    tier_used : {"empirical", "semi_analytical", "fe3d"}
        The tier that ran (``AnalysisConfig.tier``).
    config_snapshot : dict
        ``dataclasses.asdict`` of the ``AnalysisConfig`` (a custom material
        expanded to a dict), for provenance.
    buckling_eigenvalues : list of float, optional
        One value: N/mm for semi_analytical CAI, MPa for fe3d CAI. ``None``
        for the empirical tier, for tension (except semi_analytical under a
        lateral compression, where a film can buckle and grow), and when
        there is no delamination (semi_analytical) or no positive buckling
        stress (fe3d). On a biaxial semi_analytical run the force is along
        the panel axis the film is more compressed along.
    critical_sublaminate : int, optional
        Interface index from
        :func:`~bvidfe.analysis.semi_analytical.weakest_sublaminate_growth`
        for semi_analytical and fe3d CAI, and semi_analytical and fe3d TAI
        when a laterally compressed film can grow; ``None`` otherwise.
    field_results : FieldResults, optional
        Always ``None`` at present.
    notes : list of str
        Human-readable diagnostics (see above).
    warnings : list of str
        Every ``UserWarning`` raised during ``BvidAnalysis.run``, as
        ``"<WarningClass>: <message>"`` (see below).
    pristine_E11_MPa, pristine_E22_MPa, pristine_G12_MPa : float, optional
        CLT effective membrane moduli (MPa) of the intact laminate, set by
        ``BvidAnalysis.run`` for every tier.
    residual_E11_MPa, residual_E22_MPa, residual_G12_MPa : float, optional
        Residual moduli (MPa); only ``residual_E11_MPa`` for fe3d, ``None``
        otherwise.
    E11_knockdown, E22_knockdown, G12_knockdown : float, optional
        ``residual / pristine`` for each modulus, ``None`` with the residual.
    pristine_transverse_stress_MPa, residual_transverse_stress_MPa : float, optional
        Biaxial runs only (``AnalysisConfig.transverse_stress_ratio = beta
        != 0``, experimental): the signed sigma_y (MPa) at pristine and at
        damaged failure, ``s * beta * pristine_strength_MPa`` and ``s * beta
        * residual_strength_MPa`` (``s = -1`` in compression, ``+1`` in
        tension). ``None`` at ``beta == 0``.
    biaxial_interaction_factor : float, optional
        Biaxial runs only: the damaged axial strength along the ray over
        that of the same damaged panel under the uniaxial load,
        ``residual_strength_MPa(beta) / residual_strength_MPa(0)`` (1 = the
        lateral stress has no effect, below 1 = it costs strength).
        ``None`` at ``beta == 0``.

    With ``beta != 0`` the strength fields are ``|sigma_x|`` at failure
    along the proportional load ray ``(sigma_x, sigma_y) = lambda * (s, s *
    beta)``: ``pristine_strength_MPa`` the ply-discount envelope along the
    ray, ``residual_strength_MPa`` the tier's damaged strength along it, and
    ``knockdown`` their ratio as before. No biaxial residual-strength data
    validate these; see "Biaxial loading" in docs/physics_models.md.
    """

    residual_strength_MPa: float
    pristine_strength_MPa: float
    knockdown: float
    damage: DamageState
    dpa_mm2: float
    tier_used: TierName
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
    #: CLT effective membrane modulus along x (MPa) of the intact laminate,
    #: ``Laminate.effective_engineering_constants()[0]``; every tier.
    pristine_E11_MPa: Optional[float] = None
    #: Damaged panel's axial modulus along x (MPa): ``pristine_E11_MPa *
    #: E11_knockdown``. fe3d only (``None`` for the closed-form tiers).
    residual_E11_MPa: Optional[float] = None
    #: fe3d's damaged over undamaged secant axial modulus of the panel, in
    #: the run's loading sense; ``None`` for the closed-form tiers.
    E11_knockdown: Optional[float] = None
    #: CLT effective membrane modulus along y (MPa) of the intact laminate.
    pristine_E22_MPa: Optional[float] = None
    #: Not computed by any tier yet: always ``None``.
    residual_E22_MPa: Optional[float] = None
    #: Not computed by any tier yet: always ``None``.
    E22_knockdown: Optional[float] = None
    #: CLT effective in-plane shear modulus (MPa) of the intact laminate.
    pristine_G12_MPa: Optional[float] = None
    #: Not computed by any tier yet: always ``None``.
    residual_G12_MPa: Optional[float] = None
    #: Not computed by any tier yet: always ``None``.
    G12_knockdown: Optional[float] = None
    #: Signed sigma_y (MPa) at pristine failure along the biaxial load ray;
    #: ``None`` unless ``transverse_stress_ratio != 0``.
    pristine_transverse_stress_MPa: Optional[float] = None
    #: Signed sigma_y (MPa) at damaged failure along the biaxial load ray;
    #: ``None`` unless ``transverse_stress_ratio != 0``.
    residual_transverse_stress_MPa: Optional[float] = None
    #: ``residual_strength_MPa`` over the same run's uniaxial residual
    #: strength; ``None`` unless ``transverse_stress_ratio != 0``.
    biaxial_interaction_factor: Optional[float] = None

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
        for name in ("E11", "E22", "G12"):
            pristine = getattr(self, f"pristine_{name}_MPa")
            if pristine is None:
                continue
            residual = getattr(self, f"residual_{name}_MPa")
            kd = getattr(self, f"{name}_knockdown")
            line = f"  {name + '_MPa':23s}: {pristine:.0f} pristine"
            if residual is not None and kd is not None:
                line += f", {residual:.0f} residual (knockdown {kd:.3f})"
            else:
                line += ", residual not modelled by this tier"
            lines.append(line)
        if self.biaxial_interaction_factor is not None:
            ratio = self.config_snapshot.get("transverse_stress_ratio")
            label = "" if ratio is None else f" (sigma_y/sigma_x = {ratio:+.2f})"
            lines.append(f"  biaxial{label}, experimental:")
            for key in ("pristine_transverse_stress_MPa", "residual_transverse_stress_MPa"):
                value = getattr(self, key)
                if value is not None:
                    lines.append(f"    {key:31s}: {value:.1f}")
            lines.append(
                f"    {'biaxial_interaction_factor':31s}: {self.biaxial_interaction_factor:.3f}"
            )
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
        ``knockdown``, ``dpa_mm2``, ``tier_used``, ``damage``,
        ``buckling_eigenvalues``, ``critical_sublaminate``,
        ``config_snapshot``, ``notes``, ``warnings``, the nine stiffness
        fields (``pristine_E11_MPa``, ``residual_E11_MPa``,
        ``E11_knockdown`` and the same for ``E22`` and ``G12``) and the three
        biaxial fields (``pristine_transverse_stress_MPa``,
        ``residual_transverse_stress_MPa``, ``biaxial_interaction_factor``),
        ``None`` as ``null``. ``field_results`` is left out. ``damage`` is
        :func:`~bvidfe.damage.io.damage_state_to_dict`, the schema-v1.0
        C-scan layout with ``schema_version``, so it can be saved and read
        back with ``load_cscan_json`` (or the CLI's ``--cscan``) for an
        inspection-driven rerun. ``config_snapshot`` holds the config's
        dataclasses as dicts, with tuples kept as tuples (``json.dumps``
        writes them as lists); ``AnalysisConfig.from_dict`` reads it back.
        """
        out = {
            "residual_strength_MPa": self.residual_strength_MPa,
            "pristine_strength_MPa": self.pristine_strength_MPa,
            "knockdown": self.knockdown,
            "dpa_mm2": self.dpa_mm2,
            "tier_used": self.tier_used,
            "damage": damage_state_to_dict(self.damage),
            "buckling_eigenvalues": self.buckling_eigenvalues,
            "critical_sublaminate": self.critical_sublaminate,
            "config_snapshot": self.config_snapshot,
            "notes": list(self.notes),
            "warnings": list(self.warnings),
        }
        for name in ("E11", "E22", "G12"):
            for key in (f"pristine_{name}_MPa", f"residual_{name}_MPa", f"{name}_knockdown"):
                out[key] = getattr(self, key)
        for key in (
            "pristine_transverse_stress_MPa",
            "residual_transverse_stress_MPa",
            "biaxial_interaction_factor",
        ):
            out[key] = getattr(self, key)
        return out
