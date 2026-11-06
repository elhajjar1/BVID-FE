"""BvidAnalysis: orchestrates impact-to-damage and tier dispatch."""

from __future__ import annotations

import math
import warnings
from contextlib import nullcontext
from copy import deepcopy
from dataclasses import asdict, replace
from typing import ContextManager, Union

from bvidfe._types import LoadingMode
from bvidfe.analysis.config import AnalysisConfig
from bvidfe.analysis.fe_tier import (
    fe3d_basis_cache,
    fe3d_cai_buckling,
    fe3d_cai_channels_and_modulus,
    fe3d_tai_and_modulus,
    tai_fibre_damage_hole,
)
from bvidfe.analysis.results import AnalysisResults
from bvidfe.analysis.semi_analytical import (
    SemiAnalyticalResult,
    SublaminateGrowth,
    _semi_analytical_tai_result,
    panel_buckling_load,
    semi_analytical_cai,
    semi_analytical_tai,
    weakest_sublaminate_growth,
)
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY, OrthotropicMaterial
from bvidfe.damage.state import DamageState
from bvidfe.failure.laminate_strength import unnotched_strength
from bvidfe.failure.lekhnitskii import biaxial_point_stress_factor
from bvidfe.failure.soutis_openhole import (
    SOUTIS_REFERENCE_AREA_MM2,
    fibre_damage_tai,
    lekhnitskii_kt_infinity,
    soutis_cai,
    tai_hole_radius_mm,
    transverse_damage_width_mm,
)
from bvidfe.impact.mapping import impact_to_damage


class DamageOutsidePanelWarning(UserWarning):
    """A delamination's centroid lies outside the panel. The closed-form
    tiers still count its area in the DPA, but the fe3d mesh (which spans
    the panel only) may not see it at all."""


# Registry for re-emitting the warnings ``BvidAnalysis.run`` records, so the
# caller's "default" (once per location) action deduplicates repeats within a
# run (e.g. the damaged and undamaged fe3d solves). Python clears it whenever
# the filters change, which ``catch_warnings`` does, so each run shows its own.
_REEMIT_REGISTRY: dict = {}


def _resolve_material(m: Union[str, OrthotropicMaterial]) -> OrthotropicMaterial:
    if isinstance(m, str):
        try:
            return MATERIAL_LIBRARY[m]
        except KeyError:
            raise ValueError(
                f"material {m!r} is not in MATERIAL_LIBRARY; choose one of "
                f"{sorted(MATERIAL_LIBRARY)} (case-sensitive)"
            ) from None
    return m


def _warn_damage_outside_panel(damage: DamageState, Lx_mm: float, Ly_mm: float) -> None:
    """Warn when a delamination centroid is outside the ``[0, Lx] x [0, Ly]``
    panel frame (e.g. a C-scan in a different coordinate origin)."""
    outside = [
        (i, e.centroid_mm)
        for i, e in enumerate(damage.delaminations)
        if not (0.0 <= e.centroid_mm[0] <= Lx_mm and 0.0 <= e.centroid_mm[1] <= Ly_mm)
    ]
    if not outside:
        return
    listed = ", ".join(f"#{i} at ({x:g}, {y:g})" for i, (x, y) in outside[:5])
    more = f" and {len(outside) - 5} more" if len(outside) > 5 else ""
    warnings.warn(
        f"{len(outside)} delamination(s) have a centroid outside the "
        f"{Lx_mm:g} x {Ly_mm:g} mm panel ({listed}{more}); centroids are in the panel "
        f"frame with the origin at a corner. The empirical and semi_analytical tiers "
        f"still count their area in the DPA, but the fe3d mesh may not see them.",
        DamageOutsidePanelWarning,
        stacklevel=4,  # the caller of BvidAnalysis.run
    )


def _pristine_strength(
    lam: Laminate, loading: LoadingMode, transverse_stress_ratio: float = 0.0
) -> float:
    """Unnotched laminate strength in the loading direction (MPa): ply-discount
    last-ply failure on classical lamination theory, see
    ``bvidfe.failure.laminate_strength``; ``|sigma_x|`` along the biaxial
    load ray for ``transverse_stress_ratio != 0``."""
    return unnotched_strength(lam, loading, transverse_stress_ratio)


def _tai_note(
    lam: Laminate,
    damage: DamageState,
    panel_width_mm: float,
    transverse_stress_ratio: float = 0.0,
) -> str:
    """``AnalysisResults.notes`` entry for a closed-form TAI run: the
    equivalent hole and how little the TAI model is validated; on a biaxial
    run also the point-stress factor ``F(beta)``."""
    mat = lam.material
    width = transverse_damage_width_mm(damage, panel_width_mm)
    R = tai_hole_radius_mm(mat, damage, panel_width_mm)
    note = (
        f"TAI: equivalent hole {2 * R:.2f} mm across (tai_notch_fraction "
        f"{mat.tai_notch_fraction:g} x {width:.1f} mm footprint width across the load), "
        f"Whitney-Nuismer point stress at d0 = {mat.wn_d0_mm:g} mm, Kt_inf = "
        f"{lekhnitskii_kt_infinity(lam):.2f}, finite-width corrected for a "
        f"{panel_width_mm:g} mm wide panel. The TAI model is not validated: "
        f"tai_notch_fraction is fitted in-sample to one measured condition "
        f"(Korbelin et al. 2020, M21/T800S quasi-isotropic, 25 J)."
    )
    if transverse_stress_ratio != 0 and R > 0:
        F = biaxial_point_stress_factor(lam, R, mat.wn_d0_mm, transverse_stress_ratio)
        note += (
            f" Biaxial: the point-stress knockdown is multiplied by F = {F:.3f}, the "
            f"Lekhnitskii far-field sigma_y term of the axial stress at the "
            f"characteristic point (0, R + d0); the finite-width correction stays "
            f"uniaxial."
        )
    return note


def _biaxial_note(
    tier: str,
    loading: LoadingMode,
    transverse_stress_ratio: float,
    residual_uniaxial_MPa: float,
    interaction_factor: float,
    semi_detail: str = "",
) -> str:
    """``AnalysisResults.notes`` entry of every biaxial run (D12 of the
    #115 design): experimental and unvalidated, the load ray, and the tier's
    assumption (for ``semi_analytical`` followed by ``semi_detail``, see
    :func:`_semi_biaxial_detail`)."""
    beta = transverse_stress_ratio
    s = 1 if loading == "tension" else -1
    note = (
        f"Biaxial sigma_y/sigma_x = {beta:+.2f}: experimental, not validated (no "
        f"biaxial residual-strength data); the load is the ray (sigma_x, sigma_y) = "
        f"lambda * ({s:+d}, {s * beta:+.2f}), the strength fields are |sigma_x| at "
        f"failure along it and pristine_strength_MPa is the ply-discount envelope "
        f"along it"
    )
    if loading == "tension" and beta > 0:
        note += (
            " (under biaxial tension it may end when every ply has cracked, "
            "conservative where the fibres could carry on; see Limitations)"
        )
    if tier == "empirical" and loading == "compression":
        note += (
            "; empirical CAI assumes the Soutis knockdown does not depend on the "
            "ratio (homothetic: residual = knockdown x pristine along the ray)"
        )
    elif tier == "empirical":
        note += (
            "; empirical TAI adds the far-field sigma_y term at the hole's "
            "characteristic point (see the TAI note)"
        )
    elif tier == "semi_analytical":
        bound = (
            "the homothetic Soutis bound of the empirical tier"
            if loading == "compression"
            else "the fibre-damage hole with the far-field sigma_y term (see the TAI note)"
        )
        note += (
            "; semi_analytical puts the film's lateral stress (the applied sigma_y "
            "through the sublaminate's A; the Poisson-mismatch part is left out, as "
            "in the uniaxial model) into the clamped-ellipse buckling and the 2D "
            f"blister, and takes the lower of {bound} and the growth stress" + semi_detail
        )
    elif tier == "fe3d":
        channels = (
            "takes the lower of that and semi_analytical's biaxial growth stress"
            if loading == "compression"
            else "notches the fibre damage by 1 / (knockdown x F) (see the TAI note), and "
            "takes the lower of that and semi_analytical's biaxial growth stress when a "
            "film is compressed along some axis"
        )
        note += (
            "; fe3d superposes two linear solves of each mesh, x platens with the y "
            "edges free and y platens with the x faces held, to the ratio of the "
            "face-averaged reaction stresses (near the damage the local ratio differs), "
            "scales first-ply failure along the ray, and " + channels + "; the damage "
            "zone is softened for the sense of sigma_x alone, as in the uniaxial model, "
            "and E11_knockdown is the uniaxial (x solve) one"
        )
    return note + (
        f"; biaxial_interaction_factor {interaction_factor:.3f} = residual over the "
        f"uniaxial residual strength {residual_uniaxial_MPa:.1f} MPa."
    )


def _semi_biaxial_detail(
    sa_result: SemiAnalyticalResult,
    lam: Laminate,
    cfg: AnalysisConfig,
) -> str:
    """The growing interface and the whole-panel buckling load of a biaxial
    ``semi_analytical`` run, for :func:`_biaxial_note`."""
    growth = sa_result.growth
    if growth is None:
        detail = (
            ": no delamination grows (every film is in tension along both axes)"
            if cfg.loading == "tension"
            else ": no delamination grows"
        )
    else:
        other = "y" if growth.film_stress_axis == "x" else "x"
        detail = (
            f": interface {growth.interface_index} grows at |sigma_x| = "
            f"{growth.growth_stress_MPa:.1f} MPa along its {growth.growth_axis} axis; its "
            f"film is more compressed along {growth.film_stress_axis}, with "
            f"{growth.film_stress_ratio:+.3f} times that compression along {other}, and "
            f"buckles at a film stress of {growth.film_buckling_stress_MPa:.1f} MPa"
        )
    panel = cfg.panel
    N_panel = panel_buckling_load(
        lam,
        panel.Lx_mm,
        panel.Ly_mm,
        panel.boundary,
        loading=cfg.loading,
        transverse_stress_ratio=cfg.transverse_stress_ratio,
    )
    if math.isfinite(N_panel):
        factor = (
            ""
            if panel.boundary == "simply_supported"
            else f" times the {panel.boundary} panel factor, a uniaxial ratio"
        )
        detail += (
            f"; the whole panel buckles at |sigma_x| = {N_panel / lam.thickness_mm:.1f} MPa "
            f"along the ray (simply supported Rayleigh-Ritz{factor}; reported only, "
            f"not in the residual)"
        )
    return detail


def _fe3d_stiffness(Ex_damaged: float, Ex_undamaged: float, sense: str) -> tuple[float, str]:
    """fe3d's E11 knockdown, the damaged over the undamaged panel's secant
    axial modulus (capped at 1 against solver round-off: softening can only
    lower the reaction force at a prescribed end strain), and its
    ``AnalysisResults.notes`` line."""
    knockdown = min(1.0, Ex_damaged / Ex_undamaged)
    note = (
        f"fe3d E11 knockdown {knockdown:.3f} = damaged {Ex_damaged:.0f} MPa / undamaged "
        f"{Ex_undamaged:.0f} MPa (the whole panel's secant axial modulus in {sense}, "
        f"from the first-ply failure solves), applied to pristine_E11_MPa; fe3d does "
        f"not compute E22 or G12, and no stiffness measurement validates it"
    )
    return knockdown, note


def _growth_stress(growth: SublaminateGrowth | None) -> float:
    """The growth stress of ``growth`` (MPa), ``inf`` for ``None``."""
    return growth.growth_stress_MPa if growth is not None else math.inf


def _fe3d_growth_note(growth: SublaminateGrowth, transverse_stress_ratio: float) -> str:
    """``AnalysisResults.notes`` entry of fe3d's delamination growth
    channel; on a biaxial run also the film axis."""
    note = (
        f"fe3d delamination growth: interface {growth.interface_index} grows at "
        f"{growth.growth_stress_MPa:.1f} MPa (its sublaminate buckles at a film "
        f"stress of {growth.film_buckling_stress_MPa:.1f} MPa)"
    )
    if transverse_stress_ratio != 0:
        note += (
            f"; biaxial: |sigma_x| along the load ray, the film more compressed along "
            f"{growth.film_stress_axis} (semi_analytical's biaxial growth stress)"
        )
    return note


def _config_snapshot_dict(cfg: AnalysisConfig) -> dict:
    """Serialise an AnalysisConfig for provenance, handling the material union."""
    d = {}
    for name, value in asdict(cfg).items():
        d[name] = value
    if isinstance(cfg.material, OrthotropicMaterial):
        d["material"] = asdict(cfg.material)
    return deepcopy(d)


class BvidAnalysis:
    """Runs one BVID residual-strength analysis for an :class:`AnalysisConfig`.

    ``run()`` resolves the material, builds the laminate, obtains the damage
    state (the config's ``damage``, or :func:`~bvidfe.impact.mapping.impact_to_damage`
    on ``impact_panel`` or ``panel`` for an ``impact``), computes the
    pristine (unnotched) strength and dispatches to the selected tier.

    Examples
    --------
    >>> from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, PanelGeometry
    >>> cfg = AnalysisConfig(
    ...     material="IM7/8552",
    ...     layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ...     ply_thickness_mm=0.152,
    ...     panel=PanelGeometry(150, 100),
    ...     impact=ImpactEvent(energy_J=10.0),
    ... )
    >>> result = BvidAnalysis(cfg).run()
    >>> round(result.knockdown, 2)
    0.45
    """

    def __init__(self, config: AnalysisConfig) -> None:
        self.config = config

    def run(self) -> AnalysisResults:
        """Run the analysis and return an :class:`AnalysisResults`.

        The config is not modified, so ``run()`` can be called again.

        Every ``UserWarning`` raised during the run (e.g.
        ``DPACapClipWarning``, ``DamageOutsidePanelWarning``) is recorded as
        ``"<Category>: <message>"`` in ``AnalysisResults.warnings`` and then
        re-emitted with its original category and location, so the
        caller's warning filters (and ``pytest.warns``) still see it.

        Raises
        ------
        bvidfe.analysis.fe_tier.FE3DSizeError
            For ``tier="fe3d"`` when the mesh exceeds ``FE3D_MAX_DOF``.

        Warns
        -----
        SmallMassQuasiStaticWarning, DPACapClipWarning, UserWarning
            From the impact mapping (impact-driven configs), see
            :func:`~bvidfe.impact.mapping.impact_to_damage`.
        DamageOutsidePanelWarning
            When a delamination centroid lies outside the panel.
        """
        with warnings.catch_warnings(record=True) as caught:
            # Record everything; the caller's filters apply on re-emission.
            warnings.simplefilter("always")
            result = self._run()
        recorded: list[str] = []
        for w in caught:
            if issubclass(w.category, UserWarning):
                text = f"{w.category.__name__}: {w.message}"
                if text not in recorded:
                    recorded.append(text)
        result.warnings.extend(recorded)
        for w in caught:
            warnings.warn_explicit(
                w.message,
                w.category,
                w.filename,
                w.lineno,
                registry=_REEMIT_REGISTRY,
                source=w.source,
            )
        return result

    def _run(self) -> AnalysisResults:
        beta = self.config.transverse_stress_ratio
        mat = _resolve_material(self.config.material)
        lam = Laminate(
            material=mat,
            layup_deg=self.config.layup_deg,
            ply_thickness_mm=self.config.ply_thickness_mm,
        )
        damage = self._resolve_damage(lam)
        _warn_damage_outside_panel(damage, self.config.panel.Lx_mm, self.config.panel.Ly_mm)
        sigma_0 = _pristine_strength(lam, self.config.loading, beta)
        notes: list[str] = []
        E11_knockdown: float | None = None
        fe3d_uniaxial: float | None = None

        if self.config.tier == "empirical":
            sigma = self._empirical(lam, damage, sigma_0, beta)
            buckling_eigs = None
            critical_interface = None
            field_results = None
        elif self.config.tier == "semi_analytical":
            sa_result = self._semi_analytical(lam, damage, sigma_0, beta)
            sigma = sa_result.residual_strength_MPa
            critical_interface = sa_result.critical_interface_index
            N_cr = sa_result.critical_buckling_load_N
            buckling_eigs = [N_cr] if N_cr is not None else None
            field_results = None
        elif self.config.tier == "fe3d":
            (
                sigma,
                buckling_eigs,
                critical_interface,
                fe3d_notes,
                E11_knockdown,
                fe3d_uniaxial,
            ) = self._fe3d(lam, damage, sigma_0)
            notes.extend(fe3d_notes)
            field_results = None
        else:
            raise NotImplementedError(f"tier '{self.config.tier}' is not recognized")
        sigma_y_0: float | None = None
        sigma_y: float | None = None
        interaction_factor: float | None = None
        if beta != 0:
            # The uniaxial residual on the same damage: for the closed-form
            # tiers a second evaluation (for semi_analytical a second growth
            # solve), for fe3d from the basis X solves of the run itself.
            s = 1.0 if self.config.loading == "tension" else -1.0
            sigma_0_uniaxial = _pristine_strength(lam, self.config.loading)
            semi_detail = ""
            if fe3d_uniaxial is not None:
                sigma_uniaxial = fe3d_uniaxial
            elif self.config.tier == "semi_analytical":
                sigma_uniaxial = self._semi_analytical(
                    lam, damage, sigma_0_uniaxial
                ).residual_strength_MPa
                semi_detail = _semi_biaxial_detail(sa_result, lam, self.config)
            else:
                sigma_uniaxial = self._empirical(lam, damage, sigma_0_uniaxial, 0.0)
            interaction_factor = float(sigma / sigma_uniaxial if sigma_uniaxial > 0 else 0.0)
            sigma_y_0 = float(s * beta * sigma_0)
            sigma_y = float(s * beta * sigma)
            notes.insert(
                0,
                _biaxial_note(
                    self.config.tier,
                    self.config.loading,
                    beta,
                    sigma_uniaxial,
                    interaction_factor,
                    semi_detail,
                ),
            )
        if self.config.loading == "tension" and self.config.tier != "fe3d":
            notes.append(_tai_note(lam, damage, self.config.panel.Ly_mm, beta))
        # Residual stiffness (issue #114): CLT moduli for every tier, the
        # E11 residual from fe3d's reaction-force ratio only.
        E11_0, E22_0, G12_0, _ = lam.effective_engineering_constants()

        return AnalysisResults(
            residual_strength_MPa=sigma,
            pristine_strength_MPa=sigma_0,
            knockdown=sigma / sigma_0 if sigma_0 > 0 else 0.0,
            damage=damage,
            dpa_mm2=damage.projected_damage_area_mm2,
            tier_used=self.config.tier,
            config_snapshot=_config_snapshot_dict(self.config),
            buckling_eigenvalues=buckling_eigs,
            critical_sublaminate=critical_interface,
            field_results=field_results,
            notes=notes,
            pristine_E11_MPa=float(E11_0),
            residual_E11_MPa=None if E11_knockdown is None else float(E11_0 * E11_knockdown),
            E11_knockdown=E11_knockdown,
            pristine_E22_MPa=float(E22_0),
            pristine_G12_MPa=float(G12_0),
            pristine_transverse_stress_MPa=sigma_y_0,
            residual_transverse_stress_MPa=sigma_y,
            biaxial_interaction_factor=interaction_factor,
        )

    def _fe3d(
        self, lam: Laminate, damage: DamageState, sigma_0: float
    ) -> tuple[float, list[float] | None, int | None, list[str], float, float | None]:
        """fe3d residual strength on the shared pristine scale.

        fe3d's first-ply-failure stress (on the 3D mesh) is not the same
        quantity as ``_pristine_strength``, so dividing one by the other
        reported a knockdown below 1 for an undamaged panel. It is instead
        taken as fe3d's damaged over its own undamaged first-ply failure,
        applied to ``sigma_0``, so that ``pristine_strength_MPa`` stays the
        same for every tier and ``knockdown = residual / pristine`` still
        holds. This costs a second fe3d solve on the undamaged mesh.

        In compression the delamination growth stress is a closed-form
        far-field stress, the one ``semi_analytical`` uses, so it is compared
        with the normalised first-ply failure as it is, not rescaled:
        ``residual = min(sigma_growth, sigma_0 * FPF_damaged / FPF_undamaged)``
        (issue #166; the growth stress used to be rescaled by
        ``sigma_0 / FPF_undamaged`` too, which put fe3d above
        ``semi_analytical`` wherever growth governs). Buckling onset (whole
        panel or sublaminate) is reported but does not enter the residual;
        see ``fe3d_cai``.

        The same two solves give the E11 knockdown: the damaged over the
        undamaged panel's secant axial modulus (loaded-edge reaction force
        at the same end strain), recorded in ``notes`` with the raw moduli.

        Biaxial (``transverse_stress_ratio != 0``, experimental, issue
        #115): each mesh also gets the transverse basis solve, and the
        failure stresses are ``|sigma_x|`` at first-ply failure along the
        load ray (``fe_tier._ratio_displacement``); the growth stress is
        ``semi_analytical``'s biaxial one. The uniaxial residual, for
        ``biaxial_interaction_factor``, comes from the same meshes' x
        solves (an ``fe3d_basis_cache``), so the run solves four systems,
        not six.

        Returns ``(residual_MPa, buckling_eigenvalues, critical_interface,
        notes, E11_knockdown, uniaxial_residual_MPa)``, the last ``None``
        at ``transverse_stress_ratio == 0``.
        """
        undamaged = DamageState(delaminations=[], dent_depth_mm=0.0)
        if self.config.transverse_stress_ratio == 0:
            cache: ContextManager = nullcontext()
        else:
            # One entry: each mesh's ratio and uniaxial evaluations are
            # consecutive, so the damaged mesh is released before the
            # undamaged one is built (an enclosing sweep's cache is reused).
            cache = fe3d_basis_cache(max_entries=1)
        with cache:
            if self.config.loading != "compression":
                return self._fe3d_tai(lam, damage, undamaged, sigma_0)
            return self._fe3d_cai(lam, damage, undamaged, sigma_0)

    def _fe3d_cai(
        self, lam: Laminate, damage: DamageState, undamaged: DamageState, sigma_0: float
    ) -> tuple[float, list[float] | None, int | None, list[str], float, float | None]:
        """fe3d compression after impact, see ``_fe3d``."""
        cfg = self.config
        beta = cfg.transverse_stress_ratio
        cfg_uniaxial = replace(cfg, transverse_stress_ratio=0.0)
        growth, fpf_damaged, Ex_damaged = fe3d_cai_channels_and_modulus(cfg, damage, lam, math.inf)
        if beta != 0:
            growth_0, fpf_damaged_0, _ = fe3d_cai_channels_and_modulus(
                cfg_uniaxial, damage, lam, math.inf
            )
        _, fpf_undamaged, Ex_undamaged = fe3d_cai_channels_and_modulus(
            cfg, undamaged, lam, math.inf
        )
        uniaxial: float | None = None
        if beta != 0:
            _, fpf_undamaged_0, _ = fe3d_cai_channels_and_modulus(
                cfg_uniaxial, undamaged, lam, math.inf
            )
            uniaxial = min(
                _growth_stress(growth_0),
                _pristine_strength(lam, cfg.loading) * min(1.0, fpf_damaged_0 / fpf_undamaged_0),
            )
        sigma_fpf = sigma_0 * min(1.0, fpf_damaged / fpf_undamaged)
        notes = [
            f"fe3d first-ply failure = damaged {fpf_damaged:.1f} MPa / undamaged "
            f"{fpf_undamaged:.1f} MPa (fe3d's own failure stresses), applied to "
            f"pristine_strength_MPa: {sigma_fpf:.1f} MPa"
        ]
        E11_knockdown, stiffness_note = _fe3d_stiffness(Ex_damaged, Ex_undamaged, "compression")
        notes.append(stiffness_note)
        critical_interface = None
        sigma_growth = _growth_stress(growth)
        if growth is not None and math.isfinite(sigma_growth):
            critical_interface = growth.interface_index
            notes.append(_fe3d_growth_note(growth, beta))
        if sigma_growth < sigma_fpf:
            notes.append(
                f"fe3d CAI: delamination growth governs ({sigma_growth:.1f} MPa, below "
                f"first-ply failure at {sigma_fpf:.1f} MPa)"
            )
        elif math.isfinite(sigma_growth):
            notes.append(
                f"fe3d CAI: first-ply failure governs ({sigma_fpf:.1f} MPa, below "
                f"delamination growth at {sigma_growth:.1f} MPa)"
            )
        else:
            notes.append(
                f"fe3d CAI: first-ply failure governs ({sigma_fpf:.1f} MPa; no "
                f"delamination grows)"
            )
        _, lambda_crit, buckling_notes = fe3d_cai_buckling(cfg, damage, lam, math.inf)
        notes.extend(buckling_notes)
        buckling_eigs = [lambda_crit] if lambda_crit > 0 else None
        return (
            float(min(sigma_growth, sigma_fpf)),
            buckling_eigs,
            critical_interface,
            notes,
            E11_knockdown,
            None if uniaxial is None else float(uniaxial),
        )

    def _fe3d_tai(
        self, lam: Laminate, damage: DamageState, undamaged: DamageState, sigma_0: float
    ) -> tuple[float, list[float] | None, int | None, list[str], float, float | None]:
        """fe3d tension after impact: fe3d's damaged over undamaged
        first-ply failure (``fe3d_tai``, uncapped so that two results above
        ``_pristine_strength`` still keep their true ratio), applied to
        ``sigma_0``; the E11 knockdown in tension as in ``_fe3d``.

        Biaxial (``transverse_stress_ratio != 0``, see ``_fe3d``): the
        fibre-damage notch takes the empirical tier's ``F(beta)``, and
        when a film is compressed along some axis (``beta < 0``) the
        residual is the lower of that and ``semi_analytical``'s biaxial
        growth stress, as in compression."""
        cfg = self.config
        beta = cfg.transverse_stress_ratio
        cfg_uniaxial = replace(cfg, transverse_stress_ratio=0.0)
        sigma_damaged, Ex_damaged = fe3d_tai_and_modulus(cfg, damage, lam, math.inf)
        if beta != 0:
            sigma_damaged_0, _ = fe3d_tai_and_modulus(cfg_uniaxial, damage, lam, math.inf)
        sigma_undamaged, Ex_undamaged = fe3d_tai_and_modulus(cfg, undamaged, lam, math.inf)
        uniaxial: float | None = None
        if beta != 0:
            sigma_undamaged_0, _ = fe3d_tai_and_modulus(cfg_uniaxial, undamaged, lam, math.inf)
            uniaxial = _pristine_strength(lam, cfg.loading) * min(
                1.0, sigma_damaged_0 / sigma_undamaged_0
            )
        knockdown = min(1.0, sigma_damaged / sigma_undamaged)
        notes = [
            f"fe3d knockdown = damaged {sigma_damaged:.1f} MPa / undamaged "
            f"{sigma_undamaged:.1f} MPa (fe3d's own failure stresses), applied to "
            f"pristine_strength_MPa"
        ]
        E11_knockdown, stiffness_note = _fe3d_stiffness(Ex_damaged, Ex_undamaged, "tension")
        notes.append(stiffness_note)
        hole = tai_fibre_damage_hole(cfg, damage, lam)
        if hole is not None:
            size = (
                "the damage state's fibre-break core, not meshed"
                if damage.fiber_break_radius_mm >= hole[0]
                else "below the mesh"
            )
            biaxial = "" if beta == 0 else " (1 / (knockdown x F), F as in the TAI note)"
            notes.append(
                f"fe3d TAI: fibre damage as the closed-form equivalent hole "
                f"({2 * hole[0]:.2f} mm across, {size}), a Whitney-Nuismer "
                f"notch x{1 / hole[1]:.3f}{biaxial} on the stress of the elements under it"
            )
            notes.append(_tai_note(lam, damage, cfg.panel.Ly_mm, beta))
        sigma = sigma_0 * knockdown
        critical_interface = None
        if beta != 0:
            # The same call as semi_analytical's (_semi_analytical_tai_result).
            growth = weakest_sublaminate_growth(
                lam,
                damage,
                boundary=cfg.panel.boundary,
                loading="tension",
                transverse_stress_ratio=beta,
            )
            sigma_growth = _growth_stress(growth)
            if growth is not None and math.isfinite(sigma_growth):
                critical_interface = growth.interface_index
                notes.append(_fe3d_growth_note(growth, beta))
                governs = (
                    f"delamination growth governs ({sigma_growth:.1f} MPa, below first-ply "
                    f"failure at {sigma:.1f} MPa)"
                    if sigma_growth < sigma
                    else f"first-ply failure governs ({sigma:.1f} MPa, below delamination "
                    f"growth at {sigma_growth:.1f} MPa)"
                )
                notes.append(f"fe3d TAI: {governs}")
                sigma = min(sigma, sigma_growth)
            elif damage.delaminations:
                notes.append(
                    "fe3d TAI: no delamination grows (every film is in tension along both axes)"
                )
        return (
            float(sigma),
            None,
            critical_interface,
            notes,
            E11_knockdown,
            None if uniaxial is None else float(uniaxial),
        )

    def _resolve_damage(self, lam: Laminate) -> DamageState:
        if self.config.damage is not None:
            return self.config.damage
        assert self.config.impact is not None  # asserted in __post_init__
        panel = self.config.panel
        impact_panel = self.config.impact_panel or panel
        damage = impact_to_damage(self.config.impact, lam, impact_panel)
        if impact_panel is panel:
            return damage
        # Move the damage from the impact panel's frame to the concentric
        # test panel's.
        dx = (panel.Lx_mm - impact_panel.Lx_mm) / 2
        dy = (panel.Ly_mm - impact_panel.Ly_mm) / 2
        return replace(
            damage,
            delaminations=[
                replace(e, centroid_mm=(e.centroid_mm[0] + dx, e.centroid_mm[1] + dy))
                for e in damage.delaminations
            ],
        )

    def _semi_analytical(
        self,
        lam: Laminate,
        damage: DamageState,
        sigma_0: float,
        transverse_stress_ratio: float = 0.0,
    ) -> SemiAnalyticalResult:
        """semi_analytical residual strength for the pristine strength
        ``sigma_0`` along the load ray. Biaxial (``transverse_stress_ratio
        != 0``, experimental): the growth stress of the biaxial film load,
        against the homothetic Soutis bound in compression and the hole
        model with ``F(beta)`` in tension (where a laterally compressed
        film, ``beta < 0``, can grow too)."""
        boundary = self.config.panel.boundary
        if transverse_stress_ratio != 0:
            if self.config.loading == "compression":
                return semi_analytical_cai(
                    lam,
                    damage,
                    sigma_0,
                    boundary=boundary,
                    transverse_stress_ratio=transverse_stress_ratio,
                )
            return _semi_analytical_tai_result(
                lam,
                damage,
                sigma_0,
                self.config.panel.Ly_mm,
                boundary=boundary,
                transverse_stress_ratio=transverse_stress_ratio,
            )
        if self.config.loading == "compression":
            return semi_analytical_cai(lam, damage, sigma_0, boundary=boundary)
        # tension: no sublaminate buckling channel — wrap the scalar result
        # so callers consume a uniform ``SemiAnalyticalResult`` shape.
        sigma = semi_analytical_tai(lam, damage, sigma_0, self.config.panel.Ly_mm)
        return SemiAnalyticalResult(
            residual_strength_MPa=sigma,
            critical_interface_index=None,
            critical_buckling_load_N=None,
        )

    def _empirical(
        self,
        lam: Laminate,
        damage: DamageState,
        sigma_0: float,
        transverse_stress_ratio: float = 0.0,
    ) -> float:
        """Empirical residual strength (MPa) for the pristine strength
        ``sigma_0`` along the load ray. Biaxial (``transverse_stress_ratio
        != 0``, experimental): CAI scales the Soutis knockdown, which does
        not depend on the ratio, onto ``sigma_0`` (homothetic, D6 of the
        #115 design); TAI multiplies the hole's point-stress knockdown by
        the Lekhnitskii factor ``F(beta)`` (D7)."""
        dpa = damage.projected_damage_area_mm2
        mat = lam.material
        if self.config.loading == "compression":
            return soutis_cai(mat, dpa, SOUTIS_REFERENCE_AREA_MM2, sigma_0)
        Kt_inf = lekhnitskii_kt_infinity(lam)
        if transverse_stress_ratio == 0:
            return fibre_damage_tai(mat, damage, sigma_0, Kt_inf, self.config.panel.Ly_mm)
        return fibre_damage_tai(
            mat,
            damage,
            sigma_0,
            Kt_inf,
            self.config.panel.Ly_mm,
            transverse_stress_ratio=transverse_stress_ratio,
            laminate=lam,
        )
