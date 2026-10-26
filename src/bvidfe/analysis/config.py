"""Analysis configuration dataclasses."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Union

from bvidfe._types import _LOADING_MODES, _TIER_NAMES, LoadingMode, TierName
from bvidfe._validation import require_finite, require_positive
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.material import MATERIAL_LIBRARY, OrthotropicMaterial
from bvidfe.damage.state import DamageState
from bvidfe.impact.mapping import ImpactEvent


def _validate_choice(value: object, allowed: Iterable[str], field: str) -> None:
    """Raise ``ValueError`` if *value* is not a member of *allowed*.

    The error message lists the allowed values (sorted for deterministic
    output) and the offending value, prefixed with the *field* name, e.g.
    ``'tier must be one of [...]; got "fe3D"'``.
    """
    allowed_sorted = sorted(allowed)
    if value not in allowed_sorted:
        raise ValueError(f"{field} must be one of {allowed_sorted}; got {value!r}")


@dataclass
class MeshParams:
    """Mesh resolution parameters for the ``tier="fe3d"`` solve.

    Only the fe3d tier reads it (``AnalysisConfig`` fills in a default
    ``MeshParams()`` when ``tier="fe3d"`` and ``mesh`` is ``None``); the
    closed-form tiers ignore it. The Streamlit app's "3D Damage" view
    builds its display mesh from it (from ``MeshParams()`` when unset).

    Parameters
    ----------
    elements_per_ply : int
        Number of hex elements stacked through the thickness of each ply
        (a positive ``int``; a float such as ``2.0`` raises ``ValueError``).
        Total through-thickness elements ``nz = n_plies * elements_per_ply``.
        Default ``1`` (one element per ply) is adequate for first-ply-failure
        screening; increase for smoother through-thickness stress gradients
        at proportionally higher cost.
    in_plane_size_mm : float
        Target in-plane element edge length in millimetres (``> 0``). The
        mesh uses ``nx = max(1, ceil(Lx_mm / in_plane_size_mm))`` and
        ``ny = max(1, ceil(Ly_mm / in_plane_size_mm))`` equal elements in x
        and y, so the actual edge length is at most the target. Smaller
        values resolve the damage-zone edge more finely at quadratically
        higher element count and memory. Default ``5.0`` mm. The first-ply
        failure index is averaged over a ``FPF_AVERAGING_RADIUS_MM``
        (2.5 mm) disk (:mod:`bvidfe.analysis.fe_tier`), so at the default
        mesh the average covers only the element itself. Problems above
        ``FE3D_MAX_DOF`` (150,000 DOFs unless ``BVIDFE_FE3D_MAX_DOF`` is
        set) raise ``FE3DSizeError``; ``fe_mesh.estimate_fe_mesh_size``
        gives the size beforehand. The Streamlit app's help text advises
        ``>= 5`` mm on its 1 GB-RAM cloud deployment.
    cohesive_zone_factor : float
        Reserved (``> 0``, default ``1.0``). Intended as a multiplier on
        the cohesive-zone characteristic length once true cohesive
        surfaces land (see README "Limitations"); **validated but not read
        by any code path** — the present fe3d tier uses a component-wise
        stiffness-reduction model, not cohesive elements, so changing it
        does not change any result. Leave at the default.

    Raises
    ------
    ValueError
        If ``elements_per_ply`` is not a positive ``int``, or
        ``in_plane_size_mm`` or ``cohesive_zone_factor`` is not ``> 0``.

    Examples
    --------
    >>> MeshParams(in_plane_size_mm=2.5)  # finer than the 5 mm default
    MeshParams(elements_per_ply=1, in_plane_size_mm=2.5, cohesive_zone_factor=1.0)
    """

    elements_per_ply: int = 1
    in_plane_size_mm: float = 5.0
    cohesive_zone_factor: float = 1.0

    def __post_init__(self) -> None:
        if (
            isinstance(self.elements_per_ply, bool)
            or not isinstance(self.elements_per_ply, int)
            or self.elements_per_ply <= 0
        ):
            raise ValueError(
                f"MeshParams.elements_per_ply must be a positive int "
                f"(got {self.elements_per_ply!r})"
            )
        require_positive(self.in_plane_size_mm, "MeshParams.in_plane_size_mm")
        require_positive(self.cohesive_zone_factor, "MeshParams.cohesive_zone_factor")


@dataclass
class AnalysisConfig:
    """BVID analysis configuration. Provide exactly ONE of `impact` or `damage`.

    ``ply_thickness_mm`` may be either a single positive ``float`` (uniform
    laminate) or a list/tuple of positive floats with one entry per ply
    (matching ``len(layup_deg)``). Per-ply thicknesses let users model
    laminates that mix plies of different fabric weights or prepreg gauges.

    ``impact_panel`` is the panel the impact is applied on, when it is not
    the panel being tested: e.g. a plate impacted on a 60 x 60 mm clamped
    window and then trimmed to an 80 x 50 mm compression coupon. The impact
    is mapped to damage on ``impact_panel`` and the damage is placed on
    ``panel``, the two panels taken as concentric. Only valid with
    ``impact``; defaults to ``panel``.

    Parameters
    ----------
    material : str or OrthotropicMaterial
        A ``MATERIAL_LIBRARY`` preset name (e.g. ``"IM7/8552"``) or a
        custom :class:`~bvidfe.core.material.OrthotropicMaterial`. An
        unknown preset name raises ``ValueError`` listing the presets.
    layup_deg : list of float
        Ply angles in degrees, in stacking order from ply 0 (the ``z = 0``
        face of the fe3d mesh); ``0`` is the x (load) direction.
    ply_thickness_mm : float or sequence of float
        Ply thickness(es) in mm, see above.
    panel : PanelGeometry
        The panel that is loaded; ``Lx_mm`` is along the load direction.
    loading : {"compression", "tension"}
        ``"compression"`` (CAI, default) or ``"tension"`` (TAI).
    tier : {"empirical", "semi_analytical", "fe3d"}
        Residual-strength model (default ``"empirical"``).
    impact : ImpactEvent, optional
        Impact-driven workflow: the damage is predicted from the impact.
    damage : DamageState, optional
        Inspection-driven workflow: a measured (e.g. C-scan) damage state in
        the frame of ``panel``.
    mesh : MeshParams, optional
        fe3d mesh resolution; set to ``MeshParams()`` for ``tier="fe3d"``
        when omitted, ignored by the other tiers.
    impact_panel : PanelGeometry, optional
        Panel the impact is applied on, see above. Its ``boundary`` is the
        one the impact model uses.

    Raises
    ------
    ValueError
        For an unknown ``tier``, ``loading`` or material preset; unless
        exactly one of ``impact`` and ``damage`` is set; for
        ``impact_panel`` without ``impact``; for an empty or string layup,
        a non-finite ply angle, a non-finite or non-positive ply thickness
        or a thickness list whose length is not ``len(layup_deg)``; for a
        ``panel``, ``impact``, ``damage`` or ``mesh`` of the wrong type; and
        for a delamination ``interface_index`` past the last interface
        (``len(layup_deg) - 2``).

    Examples
    --------
    >>> from bvidfe import AnalysisConfig, ImpactEvent, PanelGeometry
    >>> cfg = AnalysisConfig(
    ...     material="IM7/8552",
    ...     layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ...     ply_thickness_mm=0.152,
    ...     panel=PanelGeometry(150, 100),
    ...     impact=ImpactEvent(energy_J=10.0),
    ... )
    >>> cfg.tier, cfg.loading
    ('empirical', 'compression')
    """

    material: Union[str, OrthotropicMaterial]
    layup_deg: List[float]
    ply_thickness_mm: Union[float, Sequence[float]]
    panel: PanelGeometry
    loading: LoadingMode = "compression"
    tier: TierName = "empirical"
    impact: Optional[ImpactEvent] = None
    damage: Optional[DamageState] = None
    mesh: Optional[MeshParams] = None
    impact_panel: Optional[PanelGeometry] = None

    def __post_init__(self) -> None:
        _validate_choice(self.tier, _TIER_NAMES, "tier")
        _validate_choice(self.loading, _LOADING_MODES, "loading")
        if isinstance(self.material, str):
            if self.material not in MATERIAL_LIBRARY:
                raise ValueError(
                    f"material {self.material!r} is not in MATERIAL_LIBRARY; choose one of "
                    f"{sorted(MATERIAL_LIBRARY)} (case-sensitive) or pass an "
                    f"OrthotropicMaterial"
                )
        elif not isinstance(self.material, OrthotropicMaterial):
            raise ValueError(
                f"material must be a MATERIAL_LIBRARY name or an OrthotropicMaterial "
                f"(got {type(self.material).__name__})"
            )
        if isinstance(self.layup_deg, (str, bytes)):
            raise ValueError(
                f"layup_deg must be a list of ply angles in degrees, not a string "
                f"(got {self.layup_deg!r}); e.g. [0, 45, -45, 90]"
            )
        try:
            n_plies = len(self.layup_deg)
        except TypeError as exc:
            raise ValueError(
                f"layup_deg must be a list of ply angles in degrees "
                f"(got {type(self.layup_deg).__name__})"
            ) from exc
        if n_plies == 0:
            raise ValueError("layup_deg must contain at least one ply")
        for i, a in enumerate(self.layup_deg):
            require_finite(a, f"layup_deg[{i}]")
        for name in ("panel", "impact_panel"):
            value = getattr(self, name)
            if value is None and name == "impact_panel":
                continue
            if not isinstance(value, PanelGeometry):
                raise ValueError(
                    f"AnalysisConfig.{name} must be a PanelGeometry, e.g. "
                    f"PanelGeometry(Lx_mm=150, Ly_mm=100) (got {type(value).__name__}: "
                    f"{value!r})"
                )
        if self.impact is not None and not isinstance(self.impact, ImpactEvent):
            raise ValueError(
                f"AnalysisConfig.impact must be an ImpactEvent (got {type(self.impact).__name__})"
            )
        if self.damage is not None and not isinstance(self.damage, DamageState):
            raise ValueError(
                f"AnalysisConfig.damage must be a DamageState (got {type(self.damage).__name__})"
            )
        if self.mesh is not None and not isinstance(self.mesh, MeshParams):
            raise ValueError(
                f"AnalysisConfig.mesh must be a MeshParams (got {type(self.mesh).__name__})"
            )
        if (self.impact is None) == (self.damage is None):
            raise ValueError(
                "Provide exactly one of AnalysisConfig.impact or AnalysisConfig.damage"
            )
        if self.impact_panel is not None and self.impact is None:
            raise ValueError(
                "AnalysisConfig.impact_panel only applies to an impact; it cannot be "
                "combined with a given damage state"
            )
        # Validate ply_thickness_mm shape against layup_deg here so callers
        # get an early, descriptive error before Laminate construction.
        if isinstance(self.ply_thickness_mm, (list, tuple)):
            if len(self.ply_thickness_mm) != len(self.layup_deg):
                raise ValueError(
                    f"ply_thickness_mm sequence length "
                    f"({len(self.ply_thickness_mm)}) must equal the number "
                    f"of plies ({len(self.layup_deg)})."
                )
            for i, t in enumerate(self.ply_thickness_mm):
                require_positive(t, f"ply_thickness_mm[{i}]")
        else:
            require_positive(self.ply_thickness_mm, "ply_thickness_mm")
        if self.damage is not None:
            # An n-ply laminate has interfaces 0 .. n-2. A larger index used
            # to be silently ignored by the closed-form tiers and to raise
            # IndexError when fe3d placed it on the mesh.
            last = len(self.layup_deg) - 2
            for d in self.damage.delaminations:
                if d.interface_index > last:
                    raise ValueError(
                        f"delamination interface_index {d.interface_index} is outside "
                        f"the {len(self.layup_deg)}-ply layup (interfaces 0-{last})."
                    )
        if self.tier == "fe3d" and self.mesh is None:
            self.mesh = MeshParams()
