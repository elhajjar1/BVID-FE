"""Empirical residual-strength models for CAI (Soutis) and TAI (Whitney-Nuismer).

CAI knockdown (after Soutis & Curtis 1996):
    sigma_CAI / sigma_0 = 1 / (1 + k_s * (DPA / A_ref)^m)
with A_ref = SOUTIS_REFERENCE_AREA_MM2, the 150 x 100 mm ASTM D7137 coupon,
and k_s, m the material's ``soutis_k_s`` and ``soutis_m``.

TAI via an equivalent open hole for the fibre damage (``fibre_damage_tai``):
the hole diameter is ``2R = tai_notch_fraction * w``, w the width of the
delamination footprint transverse to the load (y), and the strength is the
Whitney & Nuismer (1974) point-stress value with Tan's (1988) orthotropic
finite-width correction ``K_T_inf / K_T``:
    sigma_N / sigma_0 = 2 / (2 + xi^2 + 3*xi^4 - (Kt_inf - 3)*(5*xi^6 - 7*xi^8))
                        * K_T_inf / K_T(2R / W)
where xi = R / (R + d0), d0 = material characteristic distance (``wn_d0_mm``),
Kt_inf = infinite-plate stress concentration (3.0 for isotropic; for
orthotropic laminates use ``lekhnitskii_kt_infinity``) and W the panel width.
Delaminations cut no fibres, so the hole is a small fraction of the
footprint. ``whitney_nuismer_tai`` is the earlier whole-footprint form
(``R = sqrt(DPA / pi)``, infinite plate), kept for API compatibility; the
analysis tiers no longer use it.
"""

from __future__ import annotations

import math
import warnings

from bvidfe.core.laminate import Laminate
from bvidfe.core.material import OrthotropicMaterial
from bvidfe.damage.state import DamageState


def lekhnitskii_kt_infinity(lam: Laminate) -> float:
    """Infinite-plate stress concentration factor for an orthotropic laminate.

    Whitney & Nuismer's orthotropic form (Lekhnitskii anisotropic elasticity):

        Kt_inf = 1 + sqrt( 2*(sqrt(Ex/Ey) - nuxy) + Ex/Gxy )

    reduces to the isotropic value 3.0 for a quasi-isotropic laminate. The
    engineering constants are taken from ``lam.effective_engineering_constants``.

    Raises ``ValueError`` when the constants make the radicand negative, which
    only happens for a material card that violates the orthotropic stability
    bound (nu_xy < sqrt(Ex/Ey)).
    """
    Ex, Ey, Gxy, nuxy = lam.effective_engineering_constants()
    radicand = 2.0 * (math.sqrt(Ex / Ey) - nuxy) + Ex / Gxy
    if radicand < 0.0:
        raise ValueError(
            "Lekhnitskii Kt_inf is undefined for these laminate constants "
            f"(Ex={Ex:.4g}, Ey={Ey:.4g}, Gxy={Gxy:.4g}, nu_xy={nuxy:.4g}): "
            "2*(sqrt(Ex/Ey) - nu_xy) + Ex/Gxy < 0. Check the material card."
        )
    return 1.0 + math.sqrt(radicand)


#: Fixed area (mm^2) the damage area is scaled by in the CAI knockdown: the
#: 150 x 100 mm ASTM D7137 coupon, on which ``soutis_k_s`` is calibrated.
#: It used to be the analysed panel's area, which sent the knockdown to 1 as
#: the panel grew although the damage is a local notch, and overpredicted
#: the knockdown of Hasebe et al.'s 80 x 50 mm coupons.
SOUTIS_REFERENCE_AREA_MM2 = 150.0 * 100.0


def soutis_cai(
    m: OrthotropicMaterial,
    dpa_mm2: float,
    A_ref_mm2: float,
    sigma_pristine_MPa: float,
) -> float:
    """Compression-after-impact residual strength via Soutis knockdown.

    Returns ``kd * sigma_pristine_MPa`` where
    ``kd = 1 / (1 + k_s * (DPA/A_ref)^m)``; ``BvidAnalysis`` passes
    ``SOUTIS_REFERENCE_AREA_MM2``. Reaches ``sigma_pristine_MPa`` exactly when
    ``dpa_mm2 == 0``. The returned value is the numerator used by
    ``BvidAnalysis.run()`` to compute ``AnalysisResults.knockdown``; the
    pristine reference is the same for all three tiers.

    Parameters
    ----------
    m : OrthotropicMaterial
        Supplies ``soutis_k_s`` and ``soutis_m``.
    dpa_mm2 : float
        Projected damage area (mm²); ``<= 0`` returns the pristine strength.
    A_ref_mm2 : float
        Reference area (mm²) the damage area is scaled by; pass
        ``SOUTIS_REFERENCE_AREA_MM2`` to match the calibration of
        ``soutis_k_s``.
    sigma_pristine_MPa : float
        Unnotched compressive strength (MPa).

    Examples
    --------
    >>> from bvidfe import MATERIAL_LIBRARY
    >>> m = MATERIAL_LIBRARY["IM7/8552"]  # soutis_k_s = 4.8, soutis_m = 0.5
    >>> round(soutis_cai(m, 928.3, SOUTIS_REFERENCE_AREA_MM2, 601.3), 1)
    274.1
    """
    if dpa_mm2 <= 0:
        return sigma_pristine_MPa
    kd = 1.0 / (1.0 + m.soutis_k_s * (dpa_mm2 / A_ref_mm2) ** m.soutis_m)
    return kd * sigma_pristine_MPa


def _point_stress_knockdown(R_mm: float, d0_mm: float, Kt_inf: float) -> float:
    """Whitney-Nuismer point-stress knockdown of an infinite plate with a
    circular hole of radius ``R_mm`` (``<= 0`` gives 1.0).

    ``denom / 2`` is the approximate net-section stress ratio at distance
    d0 from the hole edge, which is physically never below 1. For
    Kt_inf >~ 9.2 the polynomial dips below 1 (knockdown > 1, i.e. residual
    above pristine) and for Kt_inf >~ 20.3 it crosses zero (a pole). The
    material presets stay below Kt_inf ~ 7.8 (all-0 deg layups; about 3 for
    quasi-isotropic ones); custom high-modulus cards can reach these
    values. The ratio is then clamped to 1 (knockdown 1.0) with a
    ``UserWarning`` so the result stays bounded instead of exceeding
    pristine or diverging.
    """
    if R_mm <= 0:
        return 1.0
    xi = R_mm / (R_mm + d0_mm)
    denom = 2.0 + xi**2 + 3 * xi**4 - (Kt_inf - 3.0) * (5 * xi**6 - 7 * xi**8)
    if denom < 2.0:
        raw = "undefined (denominator <= 0)" if denom <= 0 else f"{2.0 / denom:.3g}"
        warnings.warn(
            f"Whitney-Nuismer stress-ratio approximation is outside its valid "
            f"range for Kt_inf={Kt_inf:.3g} at xi={xi:.3g} (raw knockdown "
            f"{raw}); clamping the knockdown to 1.0.",
            UserWarning,
            stacklevel=3,
        )
        return 1.0
    return 2.0 / denom


def whitney_nuismer_tai(
    m: OrthotropicMaterial,
    dpa_mm2: float,
    sigma_pristine_MPa: float,
    Kt_inf: float | None = None,
) -> float:
    """Whitney-Nuismer point-stress strength of an infinite plate with an
    equivalent circular hole of the whole damage area, diameter
    ``2*sqrt(DPA/pi)``.

    This was the TAI model of the analysis tiers until it was replaced by
    :func:`fibre_damage_tai`. It is kept for API compatibility and is not
    used by ``BvidAnalysis``: delaminations cut no fibres, so a hole the
    size of the footprint (``R >> d0``) put the knockdown near its
    large-hole limit ``1 / Kt_inf`` at every BVID energy, below the CAI
    knockdown, and the infinite plate ignored the panel width.

    ``Kt_inf`` is the infinite-plate stress concentration factor; callers
    should pass the orthotropic value from ``lekhnitskii_kt_infinity``. When
    left as ``None`` it falls back to the isotropic value 3.0 and a
    ``UserWarning`` is emitted, since that is unconservative for CFRP. The
    knockdown is clamped to 1.0 (with a ``UserWarning``) where the
    point-stress polynomial leaves its valid range (Kt_inf >~ 9.2), see
    ``_point_stress_knockdown``.
    """
    if Kt_inf is None:
        warnings.warn(
            "whitney_nuismer_tai called without Kt_inf; falling back to the "
            "isotropic value 3.0. This underpredicts the notch sensitivity of "
            "orthotropic CFRP laminates. Pass lekhnitskii_kt_infinity(lam).",
            UserWarning,
            stacklevel=2,
        )
        Kt_inf = 3.0
    if dpa_mm2 <= 0:
        return sigma_pristine_MPa
    R = math.sqrt(dpa_mm2 / math.pi)
    return _point_stress_knockdown(R, m.wn_d0_mm, Kt_inf) * sigma_pristine_MPa


def finite_width_factor(hole_diameter_over_width: float, Kt_inf: float) -> float:
    """Finite-width correction ``K_T_inf / K_T`` for a central circular hole
    in an orthotropic strip (Tan 1988):

        K_T_inf / K_T = 3 (1 - x) / (2 + (1 - x)^3)
                        + 1/2 (x M)^6 (Kt_inf - 3) (1 - (x M)^2),
        M^2 = (sqrt(1 - 8 (3 (1 - x) / (2 + (1 - x)^3) - 1)) - 1) / (2 x^2),

    with ``x = 2R / W``. For ``Kt_inf = 3`` it is the isotropic correction
    ``3 (1 - x) / (2 + (1 - x)^3)``, the inverse of
    ``[2 + (1 - 2R/W)^3] / [3 (1 - 2R/W)]``. The notched strength of the
    strip is the infinite-plate strength times this factor. Returns 1.0 at
    ``x <= 0`` and 0.0 at ``x >= 1`` (no net section left).

    >>> round(finite_width_factor(0.5, 3.0), 4)
    0.7059
    """
    x = hole_diameter_over_width
    if x <= 0.0:
        return 1.0
    if x >= 1.0:
        return 0.0
    iso = 3.0 * (1.0 - x) / (2.0 + (1.0 - x) ** 3)
    M2 = (math.sqrt(1.0 - 8.0 * (iso - 1.0)) - 1.0) / (2.0 * x * x)
    xm = x * math.sqrt(M2)
    return iso + 0.5 * xm**6 * (Kt_inf - 3.0) * (1.0 - xm**2)


def transverse_damage_width_mm(damage: DamageState, panel_width_mm: float | None = None) -> float:
    """Width (mm) of the delamination footprint across the load direction:
    the y extent of the union of the ellipses (x is the load direction),
    capped at ``panel_width_mm`` when given. Ellipses far apart count as one
    footprint spanning both. 0.0 with no delaminations.

    >>> from bvidfe.damage.state import DelaminationEllipse
    >>> d = DamageState([DelaminationEllipse(0, (75.0, 50.0), 20.0, 10.0, 90.0)])
    >>> transverse_damage_width_mm(d)
    40.0
    """
    if not damage.delaminations:
        return 0.0
    lo, hi = math.inf, -math.inf
    for e in damage.delaminations:
        _, half = e.half_extents_mm()
        lo = min(lo, e.centroid_mm[1] - half)
        hi = max(hi, e.centroid_mm[1] + half)
    width = hi - lo
    if panel_width_mm is not None:
        width = min(width, panel_width_mm)
    return width


def tai_hole_radius_mm(
    m: OrthotropicMaterial, damage: DamageState, panel_width_mm: float | None = None
) -> float:
    """Radius (mm) of the equivalent hole of :func:`fibre_damage_tai`:
    ``tai_notch_fraction * transverse_damage_width_mm(damage, W) / 2``, and
    at least ``damage.fiber_break_radius_mm`` when the damage state carries
    a fibre-break core (a measured one, or the fe3d fibre-break model's)."""
    R = 0.5 * m.tai_notch_fraction * transverse_damage_width_mm(damage, panel_width_mm)
    return max(R, damage.fiber_break_radius_mm)


def fibre_damage_tai(
    m: OrthotropicMaterial,
    damage: DamageState,
    sigma_pristine_MPa: float,
    Kt_inf: float,
    panel_width_mm: float | None = None,
) -> float:
    """Tension-after-impact residual strength (MPa) from an equivalent open
    hole for the fibre damage.

    Tension strength after impact is set by the broken fibres, which sit
    in a small region under the impactor; the delaminations around it
    leave the fibres intact. The hole diameter is therefore a fraction of
    the footprint width transverse to the load,
    ``2R = m.tai_notch_fraction * w`` (see :func:`tai_hole_radius_mm`),
    and its strength is the Whitney & Nuismer (1974) point-stress value at
    ``d0 = m.wn_d0_mm`` with the orthotropic finite-width correction of
    Tan (1988) for a panel ``panel_width_mm`` wide (``None``: infinite
    plate), see :func:`finite_width_factor`.

    ``tai_notch_fraction`` (0.05) is fitted, with ``wn_d0_mm`` = 1 mm, to
    one measured condition: Körbelin, Dreiner & Fiedler (2020), M21/T800S
    quasi-isotropic coupons at 25 J (knockdown 0.889, see
    ``validation/datasets/korbelin_m21_t800s_tai.json``). The TAI model is
    otherwise unvalidated.

    Used by both the ``empirical`` and ``semi_analytical`` tiers, so those
    two tiers report identical tension knockdowns. Returns
    ``sigma_pristine_MPa`` for an undamaged state.

    >>> from bvidfe import MATERIAL_LIBRARY
    >>> from bvidfe.damage.state import DelaminationEllipse
    >>> m = MATERIAL_LIBRARY["IM7/8552"]
    >>> d = DamageState([DelaminationEllipse(0, (75.0, 50.0), 20.0, 10.0, 90.0)])
    >>> round(fibre_damage_tai(m, d, 1000.0, 3.0, panel_width_mm=100.0), 1)
    820.2
    """
    R = tai_hole_radius_mm(m, damage, panel_width_mm)
    if R <= 0:
        return sigma_pristine_MPa
    kd = _point_stress_knockdown(R, m.wn_d0_mm, Kt_inf)
    if panel_width_mm is not None:
        kd *= finite_width_factor(2.0 * R / panel_width_mm, Kt_inf)
    return kd * sigma_pristine_MPa
