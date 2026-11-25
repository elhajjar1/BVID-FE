"""fe3d fields for the CAI load-ramp animation (private, not part of the API).

The fe3d tier (:mod:`bvidfe.analysis.fe_tier`) is linear elastic up to
first-ply failure. It solves each mesh once at a reference end strain of
5% and scales that state to the load at which the failure index first
reaches 1. :func:`fe3d_field_history` returns what that scaling rests on,
for the damaged panel and for the undamaged panel that ``BvidAnalysis``
normalises it by, so that an animation can draw every load step from one
reference state per mesh.

For each mesh (:class:`Fe3dPanelFields`) it returns:

- the mesh and its boundary quads, for rendering;
- the nodal displacement ``u_ref`` (mm) at the reference load and the
  far-field stress ``sigma_ref_MPa`` that load carries (the reaction force
  on the loaded edge over the gross section);
- per element, the multiplier on that reference state at which the
  element fails, from the lowest of its Gauss points (``c_pointwise``) and
  disk-averaged (``c_averaged``, the check that sets fe3d's first-ply
  failure).

The solve is linear, so at an applied far-field stress ``sigma`` (MPa,
``|sigma_x|``):

- the displacement is ``(sigma / sigma_ref_MPa) * u_ref``;
- each element's load ratio to its own failure is
  ``(sigma / sigma_ref_MPa) / c_e``;
- the first element to reach a ratio of 1 under ``c_averaged`` is
  ``fpf_element``, at ``fpf_stress_MPa``.

For the panel pair (:class:`Fe3dFields`) it returns first-ply failure,
the delamination growth channel, which of the two governs, and the
residual strength and knockdowns. These equal what
``BvidAnalysis(cfg).run()`` reports for the same config (with
``tier="fe3d"``).

**Not called by ``run()``.** The values come from the ``fe_tier`` helpers
that ``run()`` itself calls (``_fe3d_basis``, ``_ratio_displacement``,
``_averaged_critical_multipliers``, ``_fpf_stress_and_modulus``, the
growth channels and the TAI notch). ``run()`` does not import this module,
and no ``AnalysisResults`` field changes: ``field_results`` stays ``None``.

**Cost.** The same two solves as ``run()``, or four on a biaxial config.
Each element's Gauss-point stresses are evaluated twice: once for the
disk-averaged check and once for the pointwise one. The result then holds
per-element and nodal arrays only: no Gauss-point values and no element
objects. Measured on an 8-ply QI 150 x 100 mm panel (CAI, 10 J,
``OMP_NUM_THREADS=1``), with peak RSS for the whole process:

| Mesh | Elements | DOFs | Time | Peak RSS | Returned |
|---|---|---|---|---|---|
| 5 mm | 4,800 | 17.6k | 4.2 s (``run()``: 4.0 s) | 0.19 GB (``run()``: 0.21 GB) | 2.8 MB |
| 2.5 mm | 19,200 | 67.5k | 19.9 s (``run()``: 19.0 s) | 0.48 GB (``run()``: 0.49 GB) | 10.9 MB |

At 2.5 mm, each panel's ``u_ref`` is 0.54 MB and its two multiplier
arrays are 0.15 MB each. Most of the returned size is the mesh, chiefly
its per-element DOF maps (5.8 MB). The damaged and undamaged fields share
it (the ``fe_mesh`` skeleton cache), and they also share the boundary
faces (6,400 quads, 0.26 MB).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np

from bvidfe._types import CriterionName, LoadingMode
from bvidfe.analysis.bvid import (
    BvidAnalysis,
    _fe3d_stiffness,
    _growth_stress,
    _pristine_strength,
    _resolve_material,
    _warn_damage_outside_panel,
)
from bvidfe.analysis.config import AnalysisConfig
from bvidfe.analysis.fe_mesh import FeMesh, element_grid_shape
from bvidfe.analysis.fe_tier import (
    FPF_AVERAGING_RADIUS_MM,
    _averaged_critical_multipliers,
    _fe3d_basis,
    _fe3d_cai_growth,
    _fe3d_tai_growth,
    _fpf_stress_and_modulus,
    _pointwise_critical_multipliers,
    _ratio_displacement,
    _tai_fibre_damage_notch,
)
from bvidfe.analysis.semi_analytical import SublaminateGrowth
from bvidfe.core.laminate import Laminate
from bvidfe.damage.state import DamageState

#: Which channel sets the residual strength: first-ply failure on the
#: damaged mesh, or the weakest delamination's growth stress.
GoverningChannel = Literal["first_ply_failure", "delamination_growth"]


def _frozen(a: np.ndarray) -> np.ndarray:
    """``a`` marked read-only (the fields are a snapshot)."""
    a.setflags(write=False)
    return a


@dataclass(frozen=True, eq=False)
class Fe3dPanelFields:
    """The linear fe3d first-ply-failure solve of one mesh.

    Stresses are far-field ``|sigma_x|`` applied to this mesh (MPa, fe3d's
    own scale; see :class:`Fe3dFields` for the scale ``BvidAnalysis``
    reports). A multiplier ``c`` is relative to the reference load: the
    state at multiplier ``c`` is ``c * u_ref``, at the stress
    ``c * sigma_ref_MPa``.

    Attributes
    ----------
    mesh : FeMesh
        The structured hex mesh: ``node_coords`` (``(n_nodes, 3)``, mm),
        ``element_connectivity`` (``(n_elements, 8)``), and the per-element
        ``ply_indices``, ``ply_angles_deg``, ``in_plane_damage_factors`` and
        ``damage_factors`` (out-of-plane), each in ``(0, 1]``. Element
        ``e = k * ny * nx + j * nx + i`` (``grid_shape``). Its arrays are
        shared with ``fe_mesh``'s skeleton cache: do not modify them.
    boundary_faces : numpy.ndarray
        ``(n_faces, 5)`` int: the 4 node indices of each boundary quad,
        then its element (``viz.plotly_3d._boundary_faces``). Use
        ``boundary_faces[:, :4]`` for the quads and ``boundary_faces[:, 4]``
        to look up a per-element value.
    grid_shape : tuple of int
        ``(nz, ny, nx)``: a per-element array reshapes to it.
    u_ref : numpy.ndarray
        ``(n_dof,)`` nodal displacement (mm) at the reference load,
        node-major: DOF ``3 n + d`` is node ``n``'s ``u_x``, ``u_y``,
        ``u_z`` for ``d = 0, 1, 2`` (``displacement`` reshapes it).
        Biaxial: ``u_X + b u_Y`` of the two basis solves, ``b`` such that
        the face-averaged stresses are in the ratio ``sigma_y = beta
        sigma_x`` (``fe_tier._ratio_displacement``).
    reference_strain : float
        Signed end strain of the x platens at ``u_ref``: ``-0.05`` in
        compression, ``+0.05`` in tension.
    sigma_ref_MPa : float
        ``|sigma_x|`` at ``u_ref``: ``reaction_force_ref_N /
        gross_section_mm2``.
    sigma_y_ref_MPa : float
        Signed face-averaged ``sigma_y`` at ``u_ref`` (MPa; 0 on a
        uniaxial config, where the ``y`` faces are free).
    reaction_force_ref_N : float
        ``|F_x|`` on the loaded ``x_max`` face at ``u_ref`` (N).
    end_shortening_ref_mm : float
        Displacement of the ``x_max`` platen at ``u_ref``,
        ``|reference_strain| * Lx`` (mm; an extension in tension). The
        ``x_min`` platen stays at ``u_x = 0``, and on a biaxial config the
        transverse basis holds both x faces at ``u_x = 0``.
    gross_section_mm2 : float
        ``Ly * h``, the section ``sigma_x`` is averaged over.
    c_pointwise : numpy.ndarray
        ``(n_elements,)``: per element, the smallest multiplier at which
        the failure index at one of its 8 Gauss points reaches 1 (``inf``
        where none is loaded towards failure), the pointwise check of
        ``fe_tier._pointwise_critical_multipliers``. With the TAI notch,
        divided by ``notch`` (the Gauss-point check of the notched stress).
    c_averaged : numpy.ndarray
        ``(n_elements,)``: per element, the multiplier at which its failure
        index, averaged over its Gauss points and then over the
        ``FPF_AVERAGING_RADIUS_MM`` (2.5 mm) disk around it, reaches 1
        (``inf`` where it never does), the check that sets fe3d's
        first-ply failure (``fe_tier._averaged_critical_multipliers``, with
        the TAI notch). Its minimum, capped at 1, is ``fpf_multiplier``.
    notch : numpy.ndarray or None
        ``(n_elements,)``: the TAI fibre-damage notch on each element's
        stress (``1 / knockdown`` under the hole, 1 elsewhere;
        ``fe_tier._tai_fibre_damage_notch``), ``None`` in compression or
        without a hole.
    fpf_multiplier : float
        ``min(1, c_averaged.min())``: the multiplier at first-ply failure,
        1 (the reference strain) when nothing fails by then
        (``fpf_capped``).
    fpf_element : int
        ``argmin(c_averaged)``: the element whose averaged index reaches 1
        first, which sets ``fpf_stress_MPa`` (unless ``fpf_capped``).
    fpf_stress_MPa : float
        ``|sigma_x|`` at first-ply failure, ``fpf_multiplier *
        sigma_ref_MPa`` to round-off: the damaged or undamaged
        first-ply-failure stress of ``BvidAnalysis``'s fe3d note.
    secant_modulus_MPa : float
        The panel's secant axial modulus ``E_x`` from basis X (the
        uniaxial solve, also on a biaxial config), whose damaged over
        undamaged ratio is ``E11_knockdown``.

    Notes
    -----
    Relation between the two multipliers (``K`` the notch, 1 without):
    for every element ``e``, ``c_averaged[e] * K[e] >= min(1, min_j
    c_pointwise[j] * K[j])`` over the elements ``j`` of ``e``'s disk (the
    same layer, overlap weight > 0). Below that multiplier every Gauss
    point of the disk has an index below 1, so their weighted mean does
    too. The 1 covers the Tsai-Wu Gauss points the pointwise check skips:
    their index at the reference load is <= 0, and its quadratic term is
    positive semi-definite. At meshes of 5 mm and coarser, the disk is the
    element alone, so ``c_averaged >= min(1 / K, c_pointwise)`` element by
    element. At 2.5 mm the disk takes in the neighbours, and an element
    next to a stress peak can have ``c_averaged`` below its own
    ``c_pointwise``. Without a notch, the averaged first-ply failure is
    never below the pointwise one: on the 10 mm and 2.5 mm test panels the
    pointwise check fails 0.1-20% earlier (``min c_pointwise / min
    c_averaged`` from 0.80 to 0.999).
    """

    mesh: FeMesh
    boundary_faces: np.ndarray
    grid_shape: tuple[int, int, int]
    u_ref: np.ndarray
    reference_strain: float
    sigma_ref_MPa: float
    sigma_y_ref_MPa: float
    reaction_force_ref_N: float
    end_shortening_ref_mm: float
    gross_section_mm2: float
    c_pointwise: np.ndarray
    c_averaged: np.ndarray
    notch: Optional[np.ndarray]
    fpf_multiplier: float
    fpf_element: int
    fpf_stress_MPa: float
    secant_modulus_MPa: float

    @property
    def fpf_capped(self) -> bool:
        """True when no element fails by the reference strain, so the
        first-ply-failure stress is the reference stress."""
        return bool(self.c_averaged.min() > self.fpf_multiplier)

    def multiplier(self, stress_MPa: float) -> float:
        """The multiplier on the reference load at the applied far-field
        stress ``stress_MPa`` (``|sigma_x|``, MPa)."""
        return stress_MPa / self.sigma_ref_MPa

    def displacement(self, stress_MPa: float) -> np.ndarray:
        """Nodal displacement ``(n_nodes, 3)`` (mm) at the applied
        far-field stress ``stress_MPa``: ``(stress_MPa / sigma_ref_MPa) *
        u_ref``."""
        return self.multiplier(stress_MPa) * self.u_ref.reshape(-1, 3)

    def load_ratio(self, stress_MPa: float, *, pointwise: bool = False) -> np.ndarray:
        """Per element, the applied load over the element's own failure
        load at ``stress_MPa``: ``(stress_MPa / sigma_ref_MPa) / c_e``, with
        ``c_averaged`` (``c_pointwise`` with ``pointwise=True``); 0 where
        ``c_e`` is infinite. With ``c_averaged``, at ``fpf_stress_MPa`` it
        reaches 1 (to round-off) at ``fpf_element`` and is at most 1
        elsewhere."""
        c = self.c_pointwise if pointwise else self.c_averaged
        return self.multiplier(stress_MPa) / c

    def end_shortening_mm(self, stress_MPa: float) -> float:
        """Platen displacement (mm) at ``stress_MPa``. The load-displacement
        curve is the line through the origin and
        ``(end_shortening_ref_mm, reaction_force_ref_N)``."""
        return self.multiplier(stress_MPa) * self.end_shortening_ref_mm

    def reaction_force_N(self, stress_MPa: float) -> float:
        """``|F_x|`` (N) on the loaded face at ``stress_MPa``:
        ``stress_MPa * gross_section_mm2``."""
        return stress_MPa * self.gross_section_mm2


@dataclass(frozen=True, eq=False)
class Fe3dFields:
    """fe3d fields of a config's damaged panel and its undamaged reference.

    Two stress scales appear. fe3d's own (``damaged`` and ``undamaged``,
    and what their methods take) is the far-field stress applied to the
    mesh. ``BvidAnalysis`` reports first-ply failure on the pristine scale
    instead: ``fpf_residual_MPa = pristine_strength_MPa * min(1,
    damaged.fpf_stress_MPa / undamaged.fpf_stress_MPa)``. So a reported
    stress is the applied one times ``pristine_strength_MPa /
    undamaged.fpf_stress_MPa`` (``reported_stress`` and
    ``applied_stress`` convert). The growth stress is a closed-form stress
    on the reported scale, as is the residual.

    Attributes
    ----------
    config : AnalysisConfig
        The config (its ``tier`` is not read).
    damage : DamageState
        The damage state of the damaged panel.
    criterion : str
        The first-ply-failure criterion: ``"larc05"`` in compression,
        ``"tsai_wu"`` in tension.
    damaged, undamaged : Fe3dPanelFields
        The damaged mesh and the undamaged mesh.
    pristine_strength_MPa : float
        The unnotched strength along the load ray, ``run()``'s
        ``pristine_strength_MPa``.
    fpf_knockdown : float
        ``min(1, damaged.fpf_stress_MPa / undamaged.fpf_stress_MPa)``.
    fpf_residual_MPa : float
        First-ply failure on the reported scale, ``pristine_strength_MPa *
        fpf_knockdown``.
    growth : SublaminateGrowth or None
        The delamination growth channel: the weakest delamination along
        the compression ray in CAI, and in TAI on a biaxial ray only
        (``None`` at ``transverse_stress_ratio == 0`` or without
        delaminations).
    governing : str
        ``"delamination_growth"`` when ``growth_stress_MPa <
        fpf_residual_MPa``, else ``"first_ply_failure"`` (the
        ``"... governs"`` notes of ``run()``).
    residual_strength_MPa : float
        ``min(growth_stress_MPa, fpf_residual_MPa)``, ``run()``'s
        ``residual_strength_MPa``.
    knockdown : float
        ``residual_strength_MPa / pristine_strength_MPa``, ``run()``'s
        ``knockdown``.
    E11_knockdown : float
        ``min(1, damaged / undamaged secant_modulus_MPa)``, ``run()``'s
        ``E11_knockdown``.
    """

    config: AnalysisConfig
    damage: DamageState
    criterion: CriterionName
    damaged: Fe3dPanelFields
    undamaged: Fe3dPanelFields
    pristine_strength_MPa: float
    fpf_knockdown: float
    fpf_residual_MPa: float
    growth: Optional[SublaminateGrowth]
    governing: GoverningChannel
    residual_strength_MPa: float
    knockdown: float
    E11_knockdown: float

    @property
    def loading(self) -> LoadingMode:
        return self.config.loading

    @property
    def transverse_stress_ratio(self) -> float:
        return self.config.transverse_stress_ratio

    @property
    def growth_stress_MPa(self) -> float:
        """The growth stress (MPa, reported scale), ``inf`` without one."""
        return _growth_stress(self.growth)

    @property
    def critical_sublaminate(self) -> Optional[int]:
        """The growing interface when the growth stress is finite, as
        ``run()``'s ``critical_sublaminate``."""
        if self.growth is not None and math.isfinite(self.growth_stress_MPa):
            return self.growth.interface_index
        return None

    def reported_stress(self, applied_MPa: float) -> float:
        """An applied (fe3d) stress on the reported scale."""
        return applied_MPa * self.pristine_strength_MPa / self.undamaged.fpf_stress_MPa

    def applied_stress(self, reported_MPa: float) -> float:
        """A reported stress as the stress applied to the fe3d meshes."""
        return reported_MPa * self.undamaged.fpf_stress_MPa / self.pristine_strength_MPa

    @property
    def residual_applied_stress_MPa(self) -> float:
        """The applied stress at which the damaged panel reaches its
        residual strength: where a load ramp ends. When first-ply failure
        governs it is ``min(damaged, undamaged fpf_stress_MPa)``, so
        ``damaged.fpf_stress_MPa`` itself when ``fpf_knockdown < 1``;
        when growth governs it is ``applied_stress(growth_stress_MPa)``,
        and no element has failed by then."""
        if self.governing == "delamination_growth":
            return self.applied_stress(self.growth_stress_MPa)
        return min(self.damaged.fpf_stress_MPa, self.undamaged.fpf_stress_MPa)


def fe3d_field_history(cfg: AnalysisConfig, damage: Optional[DamageState] = None) -> Fe3dFields:
    """The fe3d fields of ``cfg`` (see the module docstring).

    ``damage`` defaults to the damage state ``BvidAnalysis(cfg).run()``
    resolves: ``cfg.damage``, or the impact mapping of ``cfg.impact``. With
    it, the strengths and knockdowns equal ``run()``'s for that config
    with ``tier="fe3d"``; ``cfg.tier`` itself is not read. Biaxial
    (``cfg.transverse_stress_ratio != 0``, experimental): the fields of
    the load ray, as ``run()``'s strengths are.

    Raises ``FE3DSizeError`` above the fe3d size cap, as ``run()`` does.
    Warnings (``DamageOutsidePanelWarning``, the mesh's
    ``FiberBreakCoreUnmeshedWarning``, the impact mapping's) are emitted
    as they arise, not recorded.
    """
    material = _resolve_material(cfg.material)
    lam = Laminate(
        material=material, layup_deg=cfg.layup_deg, ply_thickness_mm=cfg.ply_thickness_mm
    )
    if damage is None:
        damage = BvidAnalysis(cfg)._resolve_damage(lam)
    _warn_damage_outside_panel(damage, cfg.panel.Lx_mm, cfg.panel.Ly_mm)
    sigma_0 = _pristine_strength(lam, cfg.loading, cfg.transverse_stress_ratio)
    tension = cfg.loading != "compression"
    criterion: CriterionName = "tsai_wu" if tension else "larc05"

    damaged = _panel_fields(cfg, damage, lam, tension=tension, criterion=criterion)
    undamaged = _panel_fields(
        cfg,
        DamageState(delaminations=[], dent_depth_mm=0.0),
        lam,
        tension=tension,
        criterion=criterion,
        reuse_faces=damaged,
    )

    # The arithmetic of BvidAnalysis._fe3d_cai and _fe3d_tai.
    fpf_knockdown = min(1.0, damaged.fpf_stress_MPa / undamaged.fpf_stress_MPa)
    fpf_residual = sigma_0 * fpf_knockdown
    growth = _fe3d_tai_growth(cfg, damage, lam) if tension else _fe3d_cai_growth(cfg, damage, lam)
    sigma_growth = _growth_stress(growth)
    residual = float(min(sigma_growth, fpf_residual))
    governing: GoverningChannel = (
        "delamination_growth" if sigma_growth < fpf_residual else "first_ply_failure"
    )
    E11_knockdown, _ = _fe3d_stiffness(
        damaged.secant_modulus_MPa,
        undamaged.secant_modulus_MPa,
        "tension" if tension else "compression",
    )
    return Fe3dFields(
        config=cfg,
        damage=damage,
        criterion=criterion,
        damaged=damaged,
        undamaged=undamaged,
        pristine_strength_MPa=sigma_0,
        fpf_knockdown=fpf_knockdown,
        fpf_residual_MPa=fpf_residual,
        growth=growth,
        governing=governing,
        residual_strength_MPa=residual,
        knockdown=residual / sigma_0 if sigma_0 > 0 else 0.0,
        E11_knockdown=E11_knockdown,
    )


def _panel_fields(
    cfg: AnalysisConfig,
    damage: DamageState,
    lam: Laminate,
    *,
    tension: bool,
    criterion: CriterionName,
    reuse_faces: Optional[Fe3dPanelFields] = None,
) -> Fe3dPanelFields:
    """:class:`Fe3dPanelFields` of the mesh of ``damage``, from the solves
    and checks of ``fe_tier.fe3d_tai_and_modulus`` (``tension``) or
    ``_fe3d_cai_first_ply_failure_and_modulus``. ``reuse_faces``: fields
    whose boundary faces to share when the connectivity is the same."""
    basis, _ = _fe3d_basis(
        cfg,
        damage,
        lam,
        label="fe3d fields (TAI)" if tension else "fe3d fields (CAI)",
        tension=tension,
        strain_sign=1 if tension else -1,
    )
    mesh, elements = basis.mesh, basis.elements
    notch = _tai_fibre_damage_notch(cfg, damage, lam, mesh) if tension else None
    u_ref, force_x, sigma_y = _ratio_displacement(cfg, basis, lam, cfg.transverse_stress_ratio)
    c_averaged = _averaged_critical_multipliers(
        cfg, mesh, elements, u_ref, criterion, FPF_AVERAGING_RADIUS_MM, notch
    )
    c_pointwise = _pointwise_critical_multipliers(cfg, mesh, elements, u_ref, criterion)
    if notch is not None:
        c_pointwise = c_pointwise / notch
    c_min = float(c_averaged.min())
    fpf_stress, modulus = _fpf_stress_and_modulus(cfg, basis, lam, c_min, force_x)
    area = cfg.panel.Ly_mm * lam.thickness_mm

    if reuse_faces is not None and np.array_equal(
        reuse_faces.mesh.element_connectivity, mesh.element_connectivity
    ):
        faces = reuse_faces.boundary_faces
    else:
        # Lazy: the viz package imports plotly.
        from bvidfe.viz.plotly_3d import _boundary_faces

        faces = _frozen(_boundary_faces(mesh.element_connectivity))
    return Fe3dPanelFields(
        mesh=mesh,
        boundary_faces=faces,
        grid_shape=element_grid_shape(cfg),
        u_ref=_frozen(np.array(u_ref, dtype=float)),
        reference_strain=basis.strain_sign * basis.strain_cap,
        sigma_ref_MPa=abs(force_x) / area,
        sigma_y_ref_MPa=float(sigma_y),
        reaction_force_ref_N=abs(force_x),
        end_shortening_ref_mm=basis.strain_cap * cfg.panel.Lx_mm,
        gross_section_mm2=area,
        c_pointwise=_frozen(c_pointwise),
        c_averaged=_frozen(c_averaged),
        notch=None if notch is None else _frozen(np.asarray(notch, dtype=float)),
        fpf_multiplier=min(c_min, 1.0),
        fpf_element=int(np.argmin(c_averaged)),
        fpf_stress_MPa=fpf_stress,
        secant_modulus_MPa=modulus,
    )
