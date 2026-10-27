"""Olsson quasi-static impact threshold load and onset energy.

Closed-form threshold load from plate bending + delamination fracture balance
(quasi-static, large mass, one delamination):

    Pc = pi * sqrt(32 * D* * G_IIc / 3)
    D* = sqrt(D11 * D22 * (A + 1) / 2),   A = (D12 + 2 * D66) / sqrt(D11 * D22)

The onset energy is the energy the plate and the contact store when the
load reaches Pc. The plate follows the single-term large-deflection
(von Karman) relation of a centrally loaded plate, the contact a power law:

    P = k_b * w + k_m * w^3                  (plate, deflection w)
    P = K_c * delta^n                        (contact, indentation delta)
    E_onset = k_b w_c^2 / 2 + k_m w_c^4 / 4 + Pc * delta_c / (n + 1)

with ``w_c`` the real root of ``k_b w + k_m w^3 = Pc`` and ``delta_c`` the
indentation at Pc. For ``k_m = 0`` and a linear contact this is
``Pc^2 / (2 k_cb)``, ``k_cb = 1 / (1/k_b + 1/k_c)``, the model used up to the
membrane correction.

The bending stiffness ``k_b`` depends on the plate edge boundary condition
(simply-supported / clamped / free). We compute a Navier-series point-load
compliance for the simply-supported case and apply a closed-form boundary
multiplier for the other cases — this is a first-order correction that keeps
the code simple while making the GUI ``boundary`` selector physically
meaningful. Clamped plates are typically 2-3x stiffer at center than SSSS of
the same dimensions (Timoshenko & Woinowsky-Krieger §4.4); a "free" plate has
only inertial resistance to a point load, so we treat it as a soft pseudo-
support and emit a warning.

The membrane stiffness ``k_m`` (``_k_membrane``) is the single-term
stress-function Galerkin solution of the von Karman plate equations for an
orthotropic rectangular plate deflected in its fundamental mode (Levy 1942
for the isotropic simply supported plate; Chia 1980 for orthotropic
laminates): ``sin(pi x/a) sin(pi y/b)`` with in-plane movable, straight
edges for ``simply_supported``, ``sin^2(pi x/a) sin^2(pi y/b)`` with
in-plane immovable edges for ``clamped``, and none for ``free``.

The contact law depends on the impactor tip shape (Johnson 1985):
- hemispherical: Hertz sphere-on-flat, P = (4/3) E_eff sqrt(R) delta^1.5.
- flat (cylindrical punch): P = 2 E_eff R delta.
- conical: P = (2 E_eff tan(alpha) / pi) delta^2 (Love's cone).
Its energy at Pc is ``Pc delta_c / (n + 1)``: ``2 Pc delta_c / 5`` for Hertz.

References:
- Davies & Robinson (1992); Davies & Zhang (1995), Int. J. Impact Eng. 16,
  149 (isotropic form Pc^2 = 8 pi^2 E h^3 G_IIc / (9 (1 - nu^2))).
- Suemasu & Majima (1996), J. Compos. Mater. 30, 441 (n delaminations).
- Olsson (1992), AIAA J. 30, 1587 (orthotropic effective stiffness D*).
- Olsson (2001), Composites Part A, 32(9), 1207 (threshold with D*).
- Olsson, Donadon & Falzon (2006), Int. J. Solids Struct. 43, 3124
  (small-mass dynamic threshold).
- Gonzalez et al. (2011), Compos. Sci. Technol. 71, 805, Eq. (8);
  Salvetti et al. (2016), Procedia Eng. 167, 2, Eq. (4).
- Olsson (2010), Int. J. Solids Struct., 47(21).
- Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, §4 and
  ch. 13 (large deflections).
- Levy (1942), NACA Report 737, Bending of rectangular plates with large
  deflections.
- Chia (1980), Nonlinear Analysis of Plates, McGraw-Hill (orthotropic and
  laminated rectangular plates, single-term solutions).
- Shivakumar, Elber & Illg (1985), J. Appl. Mech. 52, 674 (the same
  ``k_b w + k_m w^3`` plate law for circular laminates under impact).
- Johnson (1985), Contact Mechanics, §3.4 (flat punch), §3.5 (cone),
  §4.2 (Hertz).
"""

from __future__ import annotations

import math
import warnings
from functools import lru_cache

import numpy as np

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate

NAVIER_N: int = 11  # Navier series truncation (N x N modes)

# Closed-form multipliers on the SSSS center-point bending stiffness for
# other boundary conditions. These are approximate but capture the right
# order-of-magnitude trend so the GUI selector is not inert.
#
# For an orthotropic thin plate under a central point load, exact Galerkin
# solutions give a clamped/SSSS stiffness ratio of ~2-3 for aspect ratios
# in [0.5, 2.0] — we use 2.5 as a representative value. "Free" has no
# restoring bending stiffness from edges; we treat it as 0.4 * SSSS (a
# crude approximation that nonetheless gives users a visible effect).
_BOUNDARY_BENDING_FACTOR: dict[str, float] = {
    "simply_supported": 1.0,
    "clamped": 2.5,
    "free": 0.4,
}

# Typical conical drop-weight tip half-angle (measured from the cone axis).
# A 60-degree included angle gives a 30-degree half-angle.
_CONICAL_HALF_ANGLE_DEG: float = 30.0

# Single-term membrane solution per edge condition: the assumed mode, the
# Fourier terms of its compatibility source and its edge shortening.
#
# For w = W phi(x, y) the von Karman compatibility equation
#   a22 F,xxxx + (2 a12 + a66) F,xxyy + a11 F,yyyy = w,xy^2 - w,xx w,yy
# (a = A^-1, the laminate membrane compliance) has a source
#   W^2 alpha^2 beta^2 * sum r cos(p alpha x) cos(q beta y),
# alpha = pi / a, beta = pi / b, with the (p, q, r) below. Its stress
# function F is the same sum with each term divided by
#   Delta_pq = a22 (p alpha)^4 + (2 a12 + a66) (p alpha)^2 (q beta)^2 + a11 (q beta)^4,
# which leaves no shear on the edges, no mean normal force and straight
# edges that shorten by c W^2 alpha^2 a (and c W^2 beta^2 b). Galerkin on the
# equilibrium equation (equivalently, the membrane energy of F) adds
#   K_mov = 2 a b alpha^4 beta^4 sum eps_pq r^2 / Delta_pq
# to P phi0 / W^3, eps_pq = 1/4 if p, q > 0 and 1/2 otherwise. Holding the
# edges in place (immovable) adds the uniform edge forces N0 = A (e_x, e_y),
# e = c W^2 (alpha^2, beta^2), which cancel the shortening:
#   K_imm = K_mov + 2 a b c^2 (A11 alpha^4 + 2 A12 alpha^2 beta^2 + A22 beta^4).
# For the sine mode this is Levy's (1942) single-term solution, orthotropic
# as in Chia (1980): K_mov = (a b / 64) (alpha^4 / a11 + beta^4 / a22).
_MEMBRANE_MODES: dict[str, tuple[tuple[tuple[int, int, float], ...], float]] = {
    # phi = sin(alpha x) sin(beta y)
    "sine": (((2, 0, 0.5), (0, 2, 0.5)), 1.0 / 8.0),
    # phi = sin^2(alpha x) sin^2(beta y)
    "clamped": (
        (
            (2, 0, 0.5),
            (0, 2, 0.5),
            (4, 0, -0.5),
            (0, 4, -0.5),
            (2, 2, -1.0),
            (2, 4, 0.5),
            (4, 2, 0.5),
        ),
        3.0 / 32.0,
    ),
}

# Edge condition -> (mode, in-plane immovable). A simply supported panel
# (an ASTM D7136 coupon resting on its support window under four toggle
# clamps) is free to slide in-plane: movable edges, the lower membrane
# stiffness. A clamped panel (gripped between two plates) is held in-plane.
# A free panel has no edge to stretch against.
_MEMBRANE_EDGE: dict[str, tuple[str, bool] | None] = {
    "simply_supported": ("sine", False),
    "clamped": ("clamped", True),
    "free": None,
}


@lru_cache(maxsize=64)
def _navier_basis_ssss(
    a: float, b: float, x0: float, y0: float, n_modes: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cached SSSS Navier modal basis for a point load at ``(x0, y0)``.

    Returns ``(sin2_m, sin2_n, kx, ky)`` where ``kx = m*pi/a``,
    ``ky = n*pi/b`` and ``sin2_* = sin(k*coord)**2`` for ``m, n`` in
    ``1..n_modes``. These depend only on panel size, impact location and
    mode count — not on the laminate ``D`` matrix — so a parametric sweep
    that holds geometry + location fixed (e.g. an energy sweep) reuses the
    ``~n_modes**2`` trig evaluations instead of recomputing them per point.
    """
    idx = np.arange(1, n_modes + 1, dtype=float)
    kx = idx * np.pi / a
    ky = idx * np.pi / b
    sin2_m = np.sin(kx * x0) ** 2
    sin2_n = np.sin(ky * y0) ** 2
    return sin2_m, sin2_n, kx, ky


def _k_bending_ssss(
    lam: Laminate,
    pan: PanelGeometry,
    x0: float,
    y0: float,
    n_modes: int = NAVIER_N,
) -> float:
    """Navier series point-load stiffness of a simply-supported rectangular
    orthotropic plate at (x0, y0). Returns N/mm.

    w(x0,y0)/P = (4 / (a*b)) * sum_{m,n} sin^2(m*pi*x0/a) sin^2(n*pi*y0/b) / D_mn
    D_mn = D11 * (m*pi/a)^4 + 2*(D12 + 2*D66) * (m*pi/a)^2 * (n*pi/b)^2
           + D22 * (n*pi/b)^4
    k_bending = 1 / (w/P)
    """
    _, _, D = lam.abd_matrices()
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]
    a, b = pan.Lx_mm, pan.Ly_mm
    if not (0.0 < x0 < a and 0.0 < y0 < b):
        raise ValueError(
            f"impact location ({x0}, {y0}) must lie strictly inside the panel "
            f"(0, {a}) x (0, {b}); SSSS bending compliance is singular at the boundary."
        )
    sin2_m, sin2_n, kx, ky = _navier_basis_ssss(a, b, x0, y0, n_modes)
    kx2 = kx**2
    ky2 = ky**2
    # D_mn = D11*kx^4 + 2*(D12+2*D66)*kx^2*ky^2 + D22*ky^4  (outer over m, n)
    Dmn = (
        D11 * kx2[:, None] ** 2
        + 2.0 * (D12 + 2.0 * D66) * kx2[:, None] * ky2[None, :]
        + D22 * ky2[None, :] ** 2
    )
    numer = sin2_m[:, None] * sin2_n[None, :]
    w_over_P = (4.0 / (a * b)) * float(np.sum(numer / Dmn))
    return 1.0 / w_over_P


def _k_bending(
    lam: Laminate,
    pan: PanelGeometry,
    x0: float,
    y0: float,
    n_modes: int = NAVIER_N,
) -> float:
    """Boundary-aware center-point bending stiffness (N/mm).

    Dispatches to the SSSS Navier solution and then applies a boundary-dependent
    multiplier. "free" is physically ill-defined for a central point impact so
    we warn the caller.
    """
    k_ss = _k_bending_ssss(lam, pan, x0, y0, n_modes=n_modes)
    factor = _BOUNDARY_BENDING_FACTOR.get(pan.boundary, 1.0)
    if pan.boundary == "free":
        warnings.warn(
            "PanelGeometry.boundary='free' — a fully unrestrained plate has no "
            "bending restoring force at a central point load. BVID-FE applies a "
            "0.4x SSSS-equivalent stiffness as a soft-support approximation; "
            "results should be interpreted qualitatively. Use 'clamped' or "
            "'simply_supported' for quantitative predictions.",
            UserWarning,
            stacklevel=3,
        )
    return k_ss * factor


def _k_membrane_single_term(A: np.ndarray, a: float, b: float, mode: str, immovable: bool) -> float:
    """Single-term membrane coefficient ``K`` (N/mm^3) of an ``a x b`` plate.

    ``P = k1 W + K W^3`` for a central point load ``P`` and mode amplitude
    ``W`` (the deflection under the load). ``A`` is the 3x3 CLT extensional
    stiffness (N/mm); its shear-extension terms A16, A26 are ignored, as
    D16, D26 are in ``k_b``. See ``_MEMBRANE_MODES`` for the derivation.
    """
    terms, c = _MEMBRANE_MODES[mode]
    A = np.asarray(A, dtype=float)
    A11, A12, A22, A66 = A[0, 0], A[0, 1], A[1, 1], A[2, 2]
    # Membrane compliance of the orthotropic part of A (A16 = A26 = 0).
    det = A11 * A22 - A12**2
    s11, s22, s12, s66 = A22 / det, A11 / det, -A12 / det, 1.0 / A66
    al, be = math.pi / a, math.pi / b
    total = 0.0
    for p, q, r in terms:
        pa, qb = p * al, q * be
        delta = s22 * pa**4 + (2.0 * s12 + s66) * pa**2 * qb**2 + s11 * qb**4
        eps = 0.25 if (p > 0 and q > 0) else 0.5
        total += eps * r**2 / delta
    K = 2.0 * a * b * al**4 * be**4 * total
    if immovable:
        K += 2.0 * a * b * c**2 * (A11 * al**4 + 2.0 * A12 * al**2 * be**2 + A22 * be**4)
    return float(K)


def _k_membrane(lam: Laminate, pan: PanelGeometry) -> float:
    """Membrane (large-deflection) stiffness ``k_m`` (N/mm^3) of the panel.

    The single-term coefficient of ``_k_membrane_single_term`` for the
    panel's edge condition (``_MEMBRANE_EDGE``); 0 for a ``"free"`` panel.
    It is the value for a central load and is used at every impact point:
    one assumed mode cannot represent a deflection concentrated near an
    edge, where ``K / phi(x0, y0)^4`` would grow without bound.
    """
    edge = _MEMBRANE_EDGE.get(pan.boundary)
    if edge is None:
        return 0.0
    mode, immovable = edge
    A, _, _ = lam.abd_matrices()
    return _k_membrane_single_term(A, pan.Lx_mm, pan.Ly_mm, mode, immovable)


def plate_deflection(k_b: float, k_m: float, P: float) -> float:
    """Deflection ``w`` (mm) at which ``k_b w + k_m w^3 = P``.

    ``k_b > 0`` (N/mm), ``k_m >= 0`` (N/mm^3), ``P >= 0`` (N). The cubic has
    one real root for ``k_m > 0``; it is evaluated in the hyperbolic form
    ``w = 2 sqrt(k_b / (3 k_m)) sinh(asinh(3 P / (2 k_b) sqrt(3 k_m / k_b)) / 3)``,
    which tends to ``P / k_b`` as ``k_m -> 0`` without cancellation.
    """
    if k_b <= 0.0 or k_m < 0.0 or P < 0.0:
        raise ValueError(f"need k_b > 0, k_m >= 0, P >= 0 (got {k_b}, {k_m}, {P})")
    if k_m == 0.0:
        return P / k_b
    s = math.sqrt(k_b / (3.0 * k_m))
    return 2.0 * s * math.sinh(math.asinh(1.5 * P / (k_b * s)) / 3.0)


def _contact_law(lam: Laminate, imp: ImpactorGeometry) -> tuple[float, float]:
    """Shape-aware contact law ``P = K_c * delta^n``; returns ``(K_c, n)``.

    - hemispherical: Hertz, ``K_c = (4/3) E_eff sqrt(R)``, ``n = 1.5``.
    - flat: cylindrical punch, ``K_c = 2 E_eff R``, ``n = 1`` (Johnson §3.4).
    - conical: Love's cone, ``K_c = 2 E_eff tan(alpha) / pi``, ``n = 2``
      (Johnson §3.5), ``alpha`` the half-angle from the cone axis.
    ``R`` = diameter / 2 (mm), ``K_c`` in N/mm^n.
    """
    E_eff = _contact_E_eff(lam)
    R = imp.diameter_mm / 2.0
    if imp.shape == "hemispherical":
        return (4.0 / 3.0) * math.sqrt(R) * E_eff, 1.5
    if imp.shape == "flat":
        return 2.0 * E_eff * R, 1.0
    if imp.shape == "conical":
        tan_alpha = math.tan(math.radians(_CONICAL_HALF_ANGLE_DEG))
        return 2.0 * E_eff * tan_alpha / math.pi, 2.0
    raise ValueError(f"impactor shape {imp.shape!r} not recognized")


def _contact_indentation(lam: Laminate, imp: ImpactorGeometry, P: float) -> float:
    """Contact indentation ``delta`` (mm) at load ``P`` (N)."""
    K_c, n = _contact_law(lam, imp)
    return (P / K_c) ** (1.0 / n)


def _contact_energy(lam: Laminate, imp: ImpactorGeometry, P: float) -> float:
    """Energy (N*mm) stored in the contact at load ``P``: ``P delta / (n + 1)``
    for ``P = K_c delta^n`` (``2 P delta / 5`` for Hertz)."""
    _, n = _contact_law(lam, imp)
    return P * _contact_indentation(lam, imp, P) / (n + 1.0)


def _contact_E_eff(lam: Laminate) -> float:
    """Effective Hertz-contact modulus between steel impactor and CFRP plate."""
    E_steel = 200e3  # MPa
    nu_steel = 0.3
    E_plate = lam.material.E22
    nu_plate = 0.3
    inv_E = (1 - nu_steel**2) / E_steel + (1 - nu_plate**2) / E_plate
    return 1.0 / inv_E


def _olsson_D_star(lam: Laminate) -> float:
    """Olsson's effective plate stiffness D* (N*mm) for an orthotropic laminate.

    ``D* = sqrt(D11 * D22 * (A + 1) / 2)`` with ``A = (D12 + 2 * D66) /
    sqrt(D11 * D22)`` (Olsson 1992; 2001). ``A`` is 1 for an isotropic
    plate (so ``D* = D``), about 1.0-1.25 for quasi-isotropic tape laminates
    and about 0.15 for a cross-ply, where ``D*`` is about 0.76
    ``sqrt(D11 * D22)``. It is the harmonic mean over directions of the
    directional bending stiffness ``D(theta)``, which sets a plate's static
    point-load deflection: exact for ``D11 = D22``, within a few percent
    otherwise.
    """
    _, _, D = lam.abd_matrices()
    D12, D66 = float(D[0, 1]), float(D[2, 2])
    D_geo = lam.flexural_rigidity_Deff()  # sqrt(D11 * D22), validated positive
    A = (D12 + 2.0 * D66) / D_geo
    return D_geo * math.sqrt((A + 1.0) / 2.0)


def threshold_load(lam: Laminate, pan: PanelGeometry, imp: ImpactorGeometry) -> float:
    """Quasi-static delamination threshold load Pc (Newtons).

    Load at which a single delamination (n = 1) starts to grow unstably in a
    centrally loaded plate, from the energy release rate of a plate split
    into two sublaminates (Davies & Robinson 1992; Suemasu & Majima 1996),
    with Olsson's (2001) effective stiffness for orthotropic laminates:

        Pc = pi * sqrt(32 * D* * G_IIc / 3)
        D* = sqrt(D11 * D22 * (A + 1) / 2),  A = (D12 + 2 * D66) / sqrt(D11 * D22)

    ``G_IIc`` is the mode-II interlaminar fracture toughness (N/mm, from the
    material card) and ``D11, D22, D12, D66`` the CLT bending stiffnesses
    (N*mm), so Pc is in N. For an isotropic plate ``D* = D = E h^3 / (12 (1 -
    nu^2))`` and this is Davies' ``Pc^2 = 8 pi^2 E h^3 G_IIc / (9 (1 -
    nu^2))``. It is the n = 1 case of ``pi * sqrt(32 D* G_IIc / (n + 2))``
    (Gonzalez et al. 2011, Eq. 8; Salvetti et al. 2016, Eq. 4). The threshold
    is independent of the panel size, boundary and impactor shape; those
    enter the *onset energy* via ``onset_energy``, not here.

    No dynamic factor is applied. Olsson, Donadon & Falzon (2006) found the
    threshold about 21% higher (factor 1 / sqrt(1 - 7 pi^2 / 216) = 1.213)
    for small-mass, inertia-dominated impacts; ASTM D7136 drop-weight tests
    (5.5 kg on a 150 x 100 mm coupon) and the validation datasets are
    large-mass, quasi-static impacts, where the quasi-static value applies.

    On unclustered quasi-isotropic tape laminates this value is 0.86-0.99x
    the measured onset loads (Gonzalez et al. 2011 AS4/8552; Salvetti et al.
    2016; NIAR IM7/8552; see ``docs/physics_models.md``). Up to 0.2.0 the
    code used ``pi * sqrt(8 G_IIc sqrt(D11 D22) / 9)``, Davies' formula with
    D in place of ``E h^3 / (1 - nu^2) = 12 D``, which is sqrt(12) = 3.5x low.

    Parameters
    ----------
    lam : Laminate
        Laminate carrying the material card and CLT D matrix.
    pan : PanelGeometry
        Accepted for API symmetry with `onset_energy`; not used.
    imp : ImpactorGeometry
        Accepted for API symmetry with `onset_energy`; not used.

    Returns
    -------
    float
        Threshold load Pc in Newtons. Always positive.
    """
    D_star = _olsson_D_star(lam)  # N*mm
    G_IIc = lam.material.G_IIc  # N/mm
    return math.pi * math.sqrt(32.0 * D_star * G_IIc / 3.0)


def onset_energy(
    lam: Laminate,
    pan: PanelGeometry,
    imp: ImpactorGeometry,
    location_xy_mm: tuple[float, float] | None = None,
) -> float:
    """Impact energy (J) at which BVID damage onsets.

    The energy stored in the plate and the contact when the load reaches
    the threshold load ``Pc`` (``threshold_load``)::

        E_onset = k_b w_c^2 / 2 + k_m w_c^4 / 4 + Pc delta_c / (n + 1)

    ``w_c`` solves ``k_b w + k_m w^3 = Pc`` (``plate_deflection``): ``k_b``
    is the boundary-aware Navier bending stiffness at the impact point and
    ``k_m`` the single-term membrane stiffness of the panel (``_k_membrane``),
    which stiffens a plate that deflects by more than about its thickness
    before ``Pc`` (thin or large panels). ``delta_c`` is the indentation of
    the shape-aware contact law ``P = K_c delta^n`` at ``Pc``
    (``_contact_law``). See the module docstring for the formulas.
    """
    if location_xy_mm is None:
        location_xy_mm = (pan.Lx_mm / 2, pan.Ly_mm / 2)
    x0, y0 = location_xy_mm
    Pc = threshold_load(lam, pan, imp)  # N
    k_b = _k_bending(lam, pan, x0, y0, n_modes=NAVIER_N)  # N/mm (boundary-aware)
    k_m = _k_membrane(lam, pan)  # N/mm^3 (0 for a free panel)
    w_c = plate_deflection(k_b, k_m, Pc)  # mm
    E_plate = 0.5 * k_b * w_c**2 + 0.25 * k_m * w_c**4  # N*mm
    E_contact = _contact_energy(lam, imp, Pc)  # N*mm
    return (E_plate + E_contact) * 1e-3  # J
