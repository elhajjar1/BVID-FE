"""Semi-analytical tier: sublaminate buckling and delamination growth.

The delaminated sublaminate is a plate **clamped** along the delamination
front, which is the actual ellipse (semi-axes ``major_mm`` x ``minor_mm``
rotated by ``orientation_deg``, #165): the front is held by the substrate,
so its edges are near-clamped (Chai & Babcock 1985; Shivakumar & Whitcomb
1985), which is also the edge condition the thin-film energy release rate
of :func:`sublaminate_growth_stress` assumes. Its buckling load under the
far-field compression along x comes from a Rayleigh-Ritz solution
(:func:`_clamped_ellipse_load_factor`; trial functions
``(1 - (x'/a)^2 - (y'/b)^2)^2 P_i(x'/a) P_j(y'/b)`` in the ellipse's own
frame, Legendre ``P_i``) with the full anisotropic bending stiffness,
D16 and D26 included, reduced for the sublaminate's bending-extension
coupling: ``D* = D - B A^-1 B`` (a sublaminate cut from a symmetric
laminate is in general unsymmetric, B != 0). The panel's edge condition
does not enter: the sublaminate's edges are the delamination front, not
the panel edges (#132; it used to be a simply supported enclosing
rectangle times a panel-boundary factor).

The whole-panel buckling load (:func:`panel_buckling_load`) is a
Rayleigh-Ritz solution for the simply supported rectangle with D16 and D26
(trial functions ``x (a - x) y (b - y) P_i P_j``, exact on the Navier value
when D16 = D26 = 0), times the panel boundary factor for clamped and free
panels.

Sublaminate selection: at interface ``i`` the laminate splits into plies
``0..i`` and ``i+1..n-1``; the thinner stack (by thickness, not ply count)
is the buckling sublaminate because it buckles first.

Buckling of the sublaminate does not fail the laminate; the compressive
strength is the far-field stress at which the buckled sublaminate grows its
delamination (:func:`sublaminate_growth_stress`: the energy release rate of
the post-buckled elliptical film, a 2D blister, from a nonlinear
Rayleigh-Ritz solution, :class:`_PostbuckledFilm`), minimised over the
delaminated interfaces (:func:`weakest_sublaminate_growth`).
:func:`semi_analytical_cai` takes the lower of that and the Soutis bound
of the empirical tier. In tension (:func:`semi_analytical_tai`) the tier
is the empirical fibre-damage open-hole model unchanged.
"""

from __future__ import annotations

import functools
import math
import warnings
from dataclasses import dataclass
from typing import Optional, Sequence, Union

import numpy as np
import scipy.linalg
from numpy.polynomial import legendre

from bvidfe._types import BoundaryKind
from bvidfe._validation import require_choice
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import OrthotropicMaterial
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.soutis_openhole import (
    SOUTIS_REFERENCE_AREA_MM2,
    fibre_damage_tai,
    lekhnitskii_kt_infinity,
    soutis_cai,
)


@dataclass(frozen=True)
class SemiAnalyticalResult:
    """Structured result of :func:`semi_analytical_cai`.

    Attributes
    ----------
    residual_strength_MPa : float
        Compression-after-impact residual strength in MPa. The minimum of
        the Soutis empirical knockdown and the delamination growth stress
        when damage is present; equal to the pristine strength when the
        damage state is empty.
    critical_interface_index : int | None
        Interface whose delamination grows at the lowest far-field stress
        (:func:`weakest_sublaminate_growth`). ``None`` when the damage
        state contains no delaminations.
    critical_buckling_load_N : float | None
        Sublaminate buckling force per unit width (N/mm) at the critical
        interface. ``None`` when there is no damage or when the buckling
        sublaminate is degenerate (zero plies).
    """

    residual_strength_MPa: float
    critical_interface_index: int | None = None
    critical_buckling_load_N: float | None = None


# Whole-panel buckling coefficient multiplier on the simply supported
# Rayleigh-Ritz result (:func:`panel_buckling_load`) for clamped and free
# panels: ratios of fundamental compression buckling coefficients (k) from
# Timoshenko & Gere (1961) Theory of Elastic Stability §9.2 for square
# plates, transferred approximately to the rectangular orthotropic case.
# The sublaminate buckling load does not use it (#132): the sublaminate's
# edges are the delamination front, which the substrate clamps whatever
# the panel's edge condition.
_BOUNDARY_BUCKLING_FACTOR: dict[BoundaryKind, float] = {
    "simply_supported": 1.0,
    "clamped": 1.9,
    "free": 0.5,
}

# Maximum sublaminate aspect ratio (longer / shorter semi-axis) the
# clamped-ellipse Rayleigh-Ritz basis is converged on (< 0.1% at 10:1 with
# degree 40 along the longer axis). The peanut template never comes near it
# (its aspect ratio is at most 3.25, see ``bvidfe.impact.shape_templates``),
# so the clip only fires on user-supplied C-scan ellipses: a slender one
# (e.g. b ~ 0.01 mm, a ~ 5 mm) would buckle far above the Soutis empirical
# bound, so the buckling channel would go silently inert. We clip the aspect
# ratio to this value (widening the shorter semi-axis) and emit a
# ``DegenerateThinSublaminateWarning`` so the degenerate-thin condition is
# visible rather than masquerading as a clean empirical-tier result. (It was
# 50 for the enclosing-rectangle closed form before #132.)
_MAX_SUBLAMINATE_ASPECT: float = 10.0


class DegenerateThinSublaminateWarning(UserWarning):
    """The buckling sublaminate is so slender that its aspect ratio was
    clipped to the trusted domain; the buckling channel may be inactive."""


# Mode mix of the delamination toughness used for growth. 52.1 deg is the
# Hutchinson & Suo (1992) phase angle omega for a thin film on a substrate of
# the same material loaded by an edge force alone. At the edge of a buckled
# film the phase angle moves from about 38 deg at buckling onset towards
# mode II as the film stress rises (Hutchinson & Suo give it for the strip
# and the circular blister), so a fixed value is a representative mix, not
# a tracked one. The Benzeggagh-Kenane (1996) exponent is the
# typical value for carbon/epoxy.
_GROWTH_PHASE_ANGLE_DEG: float = 52.1
_BK_EXPONENT: float = 2.0


@dataclass(frozen=True)
class SublaminateGrowth:
    """Buckling-driven growth of one delamination.

    Attributes
    ----------
    interface_index : int
        Delaminated interface.
    buckling_load_N_per_mm : float
        Buckling force per unit width of the thinner sublaminate (N/mm),
        from :func:`sublaminate_buckling_load`.
    film_buckling_stress_MPa : float
        Sublaminate stress at buckling, ``buckling_load_N_per_mm / h``.
    film_growth_stress_MPa : float
        Sublaminate stress at which the delamination grows.
    growth_stress_MPa : float
        Far-field laminate stress at which the delamination grows. ``inf``
        when the sublaminate is degenerate.
    growth_axis : str | None
        The semi-axis of the ellipse whose growth releases ``G_c`` first:
        ``"major"`` (``major_mm``, along ``orientation_deg``) or
        ``"minor"``. ``None`` when the sublaminate is degenerate.
    """

    interface_index: int
    buckling_load_N_per_mm: float
    film_buckling_stress_MPa: float
    film_growth_stress_MPa: float
    growth_stress_MPa: float
    growth_axis: Optional[str] = None


def mixed_mode_toughness(material: OrthotropicMaterial) -> float:
    """Delamination toughness (N/mm) at the thin-film mode mix.

    Benzeggagh-Kenane: ``G_c = G_Ic + (G_IIc - G_Ic) * (G_II / G)^eta``
    with ``G_II / G = sin^2(psi)``, ``psi = 52.1 deg`` and ``eta = 2``.
    0.478 N/mm for IM7/8552.
    """
    mode_ii_fraction = math.sin(math.radians(_GROWTH_PHASE_ANGLE_DEG)) ** 2
    return material.G_Ic + (material.G_IIc - material.G_Ic) * mode_ii_fraction**_BK_EXPONENT


def _thinner_sublaminate(lam: Laminate, interface_index: int) -> tuple[list[float], list[float]]:
    """Layup and per-ply thicknesses of the thinner sublaminate at an interface.

    Chooses between "above" (plies 0..i) and "below" (plies i+1..) by
    through-thickness (sum of per-ply thicknesses), not ply count — for
    non-uniform laminates the side with fewer plies may be the geometrically
    thicker one, and the *thinner* stack is what buckles first.
    """
    i = interface_index
    thicknesses = lam.ply_thicknesses_mm
    if sum(thicknesses[: i + 1]) <= sum(thicknesses[i + 1 :]):
        return list(lam.layup_deg[: i + 1]), thicknesses[: i + 1]
    return list(lam.layup_deg[i + 1 :]), thicknesses[i + 1 :]


def sublaminate_buckling_load(
    lam: Laminate,
    ellipse: DelaminationEllipse,
    boundary: BoundaryKind = "simply_supported",
) -> float:
    """Critical buckling force per unit width N_cr (N/mm) for the thinner
    sublaminate at the given delamination interface.

    The sublaminate is a plate **clamped** along the delamination front,
    the ellipse itself (semi-axes ``major_mm`` and ``minor_mm``, rotated by
    ``orientation_deg`` from the load axis x), loaded by the far-field
    compression along x. The front is held by the substrate, so it is
    near-clamped (Chai & Babcock 1985; Shivakumar & Whitcomb 1985), the
    edge condition the thin-film energy release rate of
    :func:`sublaminate_growth_stress` assumes. The Rayleigh-Ritz solution
    (:func:`_clamped_ellipse_load_factor`) uses the sublaminate's full
    anisotropic bending stiffness in the ellipse's frame, D16 and D26
    included, reduced for its bending-extension coupling:

        D* = D - B A^-1 B

    with ``A``, ``B``, ``D`` of the sublaminate on its own (a sublaminate
    cut from a symmetric laminate is in general unsymmetric; ``D* = D``
    when ``B = 0``). ``D*`` does not depend on the reference plane.

    The selected sublaminate is the *thinner* (by through-thickness, not
    ply count) of the "above" and "below" stacks at the interface (it
    buckles first).

    Before #132 this was the simply supported rectangle enclosing the
    ellipse (its bounding box in the panel frame, Navier closed form, D16,
    D26 and B ignored) times a factor from the panel's edge condition
    (clamped 1.9, free 0.5). The sublaminate's edges are the delamination
    front, not the panel edges, so ``boundary`` no longer changes the
    result; it is still validated.

    Parameters
    ----------
    lam : Laminate
        Full panel laminate; supplies material, layup, and ply thickness.
    ellipse : DelaminationEllipse
        Delamination at which the sublaminate forms. Its
        ``interface_index`` selects the sublaminate; ``major_mm`` and
        ``minor_mm`` are its semi-axes (``area_mm2 = pi * major * minor``)
        and ``orientation_deg`` the angle of ``major_mm`` from x.
    boundary : BoundaryKind
        One of ``"simply_supported"``, ``"clamped"``, ``"free"``; any other
        value raises ``ValueError``. Does not change the result (see
        above).

    The ellipse aspect ratio (longer over shorter semi-axis) is clipped to
    ``_MAX_SUBLAMINATE_ASPECT`` (10), the domain the Ritz basis is
    converged on, by widening the shorter semi-axis. A slender ellipse
    (e.g. a 0.01 mm wide C-scan slice) would otherwise buckle at a stress
    far above the Soutis bound, so the buckling channel would go silently
    inert. When the clip fires a :class:`DegenerateThinSublaminateWarning`
    is emitted so the degenerate condition is visible.

    Returns
    -------
    float
        Critical buckling force per unit sublaminate width N_cr in N/mm.
        Returns ``inf`` when the sublaminate is degenerate (zero plies,
        zero ellipse area).
    """
    # A typo used to fall back silently to the simply-supported factor.
    require_choice(boundary, _BOUNDARY_BUCKLING_FACTOR, "boundary")
    film = _sublaminate_film(lam, ellipse, warn_stacklevel=3)
    if film is None:
        return float("inf")
    D_star, _, a, b, (n_xx, n_yy, n_xy) = film
    return _clamped_ellipse_load_factor(D_star, a, b, n_xx, n_yy, n_xy)


def _sublaminate_film(
    lam: Laminate,
    ellipse: DelaminationEllipse,
    warn_stacklevel: Optional[int] = None,
) -> Optional[tuple[np.ndarray, np.ndarray, float, float, tuple[float, float, float]]]:
    """The thinner sublaminate over ``ellipse`` as a clamped elliptical
    film, in the ellipse's frame (x' along ``major_mm``): ``(D*, A, a, b,
    (n_xx, n_yy, n_xy))`` with the reduced bending stiffness ``D*``, the
    extensional stiffness ``A`` (both N, mm), the semi-axes after the
    aspect-ratio clip and the unit compression along x. ``None`` when the
    sublaminate or the ellipse is degenerate. The clip warns only when
    ``warn_stacklevel`` is given (the stack level of the warning)."""
    sub_layup, sub_thicknesses = _thinner_sublaminate(lam, ellipse.interface_index)
    if len(sub_layup) == 0:
        return None

    a = float(ellipse.major_mm)
    b = float(ellipse.minor_mm)
    if not (a > 0 and b > 0 and math.isfinite(a) and math.isfinite(b)):
        return None

    # Clip the aspect ratio to the domain the Ritz basis is converged on,
    # and surface the clip.
    aspect = max(a, b) / min(a, b)
    if aspect > _MAX_SUBLAMINATE_ASPECT:
        if warn_stacklevel is not None:
            warnings.warn(
                f"Sublaminate ellipse aspect ratio ({aspect:.1f}) exceeds the "
                f"trusted Rayleigh-Ritz domain ({_MAX_SUBLAMINATE_ASPECT:.0f}); "
                f"clipping to {_MAX_SUBLAMINATE_ASPECT:.0f} by widening the shorter "
                f"semi-axis. The sublaminate-buckling channel is degenerate for this "
                f"thin slice and may not constrain the residual strength; the "
                f"empirical Soutis tier likely governs.",
                DegenerateThinSublaminateWarning,
                stacklevel=warn_stacklevel,
            )
        if a < b:
            a = b / _MAX_SUBLAMINATE_ASPECT
        else:
            b = a / _MAX_SUBLAMINATE_ASPECT

    # Stiffness in the ellipse's frame (x' along major_mm): the plies at
    # angle phi from x are at phi - theta from x'.
    theta = float(ellipse.orientation_deg)
    local_layup = [float(phi) - theta for phi in sub_layup]
    D_star = _reduced_bending_stiffness(lam.material, local_layup, sub_thicknesses)
    A_sub, _, _ = Laminate(lam.material, local_layup, sub_thicknesses).abd_matrices()
    # The compression along x, in the ellipse's frame: x = (cos t, -sin t).
    t = math.radians(theta)
    c, s = math.cos(t), math.sin(t)
    return D_star, np.asarray(A_sub, dtype=float), a, b, (c * c, s * s, -c * s)


def _reduced_bending_stiffness(
    material: OrthotropicMaterial,
    sub_layup_deg: list[float],
    ply_thickness_mm: Union[float, Sequence[float]],
) -> np.ndarray:
    """``D* = D - B A^-1 B`` (3, 3) of a sublaminate, in N mm: the bending
    stiffness of a plate whose membrane force resultant stays zero while
    it bends (the coupled sublaminate stiffness; Yin 1985). Equal to ``D``
    for a symmetric stack, and independent of the reference plane."""
    A, B, D = Laminate(material, sub_layup_deg, ply_thickness_mm).abd_matrices()
    D_star = D - B @ np.linalg.solve(A, B)
    return np.asarray(0.5 * (D_star + D_star.T), dtype=float)


def _legendre_table(s: np.ndarray, n: int) -> np.ndarray:
    """``P_i(s)``, ``P_i'(s)`` and ``P_i''(s)`` for ``i < n``: shape
    ``(3, n, s.size)``."""
    out = np.zeros((3, n, s.size))
    for i in range(n):
        c = np.zeros(i + 1)
        c[i] = 1.0
        out[0, i] = legendre.legval(s, c)
        if i >= 1:
            out[1, i] = legendre.legval(s, legendre.legder(c, 1))
        if i >= 2:
            out[2, i] = legendre.legval(s, legendre.legder(c, 2))
    return out


# Ritz basis sizes. The clamped ellipse uses total degree 8 in its own
# frame for a circle and up to _ELLIPSE_MAX_DEGREE along the longer axis of a
# slender one (a strip buckles in several half-waves along it); with these
# sizes the eigenvalue is converged to < 0.1% up to the aspect-ratio clip
# (_MAX_SUBLAMINATE_ASPECT). The simply supported panel uses
# _PANEL_BASE_TERMS Legendre terms in y and three more in x per Navier
# half-wave along the load.
_ELLIPSE_BASE_DEGREE: int = 8
_ELLIPSE_MAX_DEGREE: int = 40
_PANEL_BASE_TERMS: int = 12
_PANEL_MAX_TERMS: int = 60

# Field names of the ellipse basis: w and its derivatives in the unit-disk
# coordinates (xi, eta) = (x'/a, y'/b).
_ELLIPSE_FIELDS = ("x", "y", "xx", "yy", "xy")


@functools.lru_cache(maxsize=32)
def _clamped_disk_integrals(p_x: int, p_y: int) -> dict[tuple[str, str], np.ndarray]:
    """Integrals over the unit disk of products of the derivative fields of
    the clamped trial functions ``(1 - xi^2 - eta^2)^2 P_i(xi) P_j(eta)``
    with ``i / p_x + j / p_y <= 1``. Keys are pairs of
    :data:`_ELLIPSE_FIELDS`; the quadrature (Gauss-Legendre in r, the
    trapezoidal rule in the polar angle) is exact for these polynomials."""
    idx = [(i, j) for i in range(p_x + 1) for j in range(p_y + 1) if i * p_y + j * p_x <= p_x * p_y]
    degree = max(i + j for i, j in idx) + 4
    rg, rw = legendre.leggauss(degree + 1)
    r = 0.5 * (rg + 1.0)
    n_t = 2 * degree + 2
    t = np.linspace(0.0, 2.0 * np.pi, n_t, endpoint=False)
    R, T = np.meshgrid(r, t, indexing="ij")
    xi = (R * np.cos(T)).ravel()
    eta = (R * np.sin(T)).ravel()
    weight = (np.outer(0.5 * rw, np.full(n_t, 2.0 * np.pi / n_t)) * R).ravel()

    Px = _legendre_table(xi, p_x + 1)
    Py = _legendre_table(eta, p_y + 1)
    ii = np.array([i for i, _ in idx])
    jj = np.array([j for _, j in idx])
    m0 = Px[0, ii] * Py[0, jj]
    mx = Px[1, ii] * Py[0, jj]
    my = Px[0, ii] * Py[1, jj]
    g = 1.0 - xi**2 - eta**2
    G0 = g * g
    Gx = -4.0 * g * xi
    Gy = -4.0 * g * eta
    Gxx = 8.0 * xi**2 - 4.0 * g
    Gyy = 8.0 * eta**2 - 4.0 * g
    Gxy = 8.0 * xi * eta
    fields = {
        "x": Gx * m0 + G0 * mx,
        "y": Gy * m0 + G0 * my,
        "xx": Gxx * m0 + 2.0 * Gx * mx + G0 * Px[2, ii] * Py[0, jj],
        "yy": Gyy * m0 + 2.0 * Gy * my + G0 * Px[0, ii] * Py[2, jj],
        "xy": Gxy * m0 + Gx * my + Gy * mx + G0 * Px[1, ii] * Py[1, jj],
    }
    return {(f, h): (fields[f] * weight) @ fields[h].T for f in fields for h in fields}


def _clamped_ellipse_load_factor(
    D: np.ndarray,
    a: float,
    b: float,
    n_xx: float,
    n_yy: float,
    n_xy: float,
) -> float:
    """Rayleigh-Ritz buckling load factor of a clamped elliptical plate.

    The plate has semi-axes ``a`` along its x' axis and ``b`` along y'
    (mm), bending stiffness ``D`` (3, 3, N mm, in the x'y' frame, D16 and
    D26 included) and the in-plane compressive force resultant
    ``lambda * (n_xx, n_yy, n_xy)`` (N/mm, compression positive, in the
    x'y' frame); returns the lowest ``lambda``. With
    ``(n_xx, n_yy, n_xy) = (1, 0, 0)`` it is the uniaxial ``N_cr``.
    Isotropic clamped circle of radius R: ``14.68 D / R^2`` under uniform
    radial compression (exact: the square of the first zero of J1) and
    ``27.2 D / R^2`` under uniaxial compression.

    Trial functions ``(1 - (x'/a)^2 - (y'/b)^2)^2 P_i(x'/a) P_j(y'/b)``
    (Legendre ``P_i``) satisfy the clamped conditions on the whole ellipse.
    For a circle they span the polynomials of total degree up to 8 times
    the clamp factor, so the result does not depend on the frame; along
    the longer axis of an elongated ellipse the degree grows with the
    aspect ratio (up to 40 at 10:1).
    """
    aspect = max(a, b) / min(a, b)
    p_long = min(_ELLIPSE_MAX_DEGREE, _ELLIPSE_BASE_DEGREE + 4 * math.ceil(aspect - 1.0 - 1e-9))
    p_x, p_y = (p_long, _ELLIPSE_BASE_DEGREE) if a >= b else (_ELLIPSE_BASE_DEGREE, p_long)
    J = _clamped_disk_integrals(p_x, p_y)
    D11, D22, D12, D66, D16, D26 = D[0, 0], D[1, 1], D[0, 1], D[2, 2], D[0, 2], D[1, 2]
    # d/dx' = (1/a) d/dxi, d/dy' = (1/b) d/deta, dA = a b dxi deta; the
    # common factor a b cancels in the eigenvalue.
    K = (
        D11 / a**4 * J["xx", "xx"]
        + D22 / b**4 * J["yy", "yy"]
        + D12 / (a * b) ** 2 * (J["xx", "yy"] + J["yy", "xx"])
        + 4.0 * D66 / (a * b) ** 2 * J["xy", "xy"]
        + 2.0 * D16 / (a**3 * b) * (J["xx", "xy"] + J["xy", "xx"])
        + 2.0 * D26 / (a * b**3) * (J["yy", "xy"] + J["xy", "yy"])
    )
    G = (
        n_xx / a**2 * J["x", "x"]
        + n_yy / b**2 * J["y", "y"]
        + n_xy / (a * b) * (J["x", "y"] + J["y", "x"])
    )
    # Largest mu of G v = mu K v (K is positive definite; G is only
    # semi-definite for a general load), lambda = 1 / mu.
    n = K.shape[0]
    mu = scipy.linalg.eigh(G, K, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0]
    return float(1.0 / mu) if mu > 0 else float("inf")


@functools.lru_cache(maxsize=32)
def _ss_beam_integrals(n_terms: int) -> np.ndarray:
    """``J[r, t] = int_-1^1 f_i^(r) f_j^(t) ds`` for ``r, t in 0..2`` and
    the simply supported trial functions ``f_i(s) = (1 - s^2) P_i(s)``,
    ``i < n_terms``: shape ``(3, 3, n_terms, n_terms)``. Gauss-Legendre
    quadrature, exact for these polynomials."""
    s, w = legendre.leggauss(n_terms + 3)
    P = _legendre_table(s, n_terms)
    q = 1.0 - s * s
    F = np.stack(
        [
            q * P[0],
            -2.0 * s * P[0] + q * P[1],
            -2.0 * P[0] - 4.0 * s * P[1] + q * P[2],
        ]
    )
    return np.einsum("rik,tjk,k->rtij", F, F, w)


def _ssss_ritz_N_cr(D: np.ndarray, a: float, b: float) -> float:
    """Rayleigh-Ritz ``N_cr`` (N/mm) of a simply supported ``a`` x ``b``
    rectangle (mm) with bending stiffness ``D`` (3, 3, N mm; D16 and D26
    included) under uniaxial compression along x (the side ``a``).

    Trial functions ``x (a - x) y (b - y) P_i(2x/a - 1) P_j(2y/b - 1)``
    satisfy ``w = 0`` on the edges; the zero edge moment is natural, so the
    Ritz value converges to the simply supported plate's. With D16 = D26 =
    0 it is the Navier value ``min over m, n of (pi^2 / a^2) [D11 m^4 +
    2 (D12 + 2 D66)(m a/b)^2 n^2 + D22 (n a/b)^4] / m^2`` (Timoshenko &
    Gere §9.2) to 1e-8; D16 and D26 only lower it (it is the minimum over a
    space containing the Navier mode, whose D16/D26 energy is zero), by
    3.0% on the canonical 200 x 150 mm 8-ply QI IM7/8552 panel (12.58 ->
    12.20 MPa). The x basis grows with the Navier half-wave count
    ``(a/b) (D22/D11)^(1/4)``.
    """
    D11, D22, D12, D66, D16, D26 = D[0, 0], D[1, 1], D[0, 1], D[2, 2], D[0, 2], D[1, 2]
    half_waves = (a / b) * (D22 / D11) ** 0.25
    n_x = min(_PANEL_MAX_TERMS, _PANEL_BASE_TERMS + 3 * math.ceil(half_waves))
    Jx = _ss_beam_integrals(n_x)
    Jy = _ss_beam_integrals(_PANEL_BASE_TERMS)

    def Ix(r: int, t: int) -> np.ndarray:
        # dx = (a/2) ds and d/dx = (2/a) d/ds.
        return Jx[r, t] * (a / 2.0) * (2.0 / a) ** (r + t)

    def Iy(r: int, t: int) -> np.ndarray:
        return Jy[r, t] * (b / 2.0) * (2.0 / b) ** (r + t)

    k = np.kron
    K = (
        D11 * k(Ix(2, 2), Iy(0, 0))
        + D22 * k(Ix(0, 0), Iy(2, 2))
        + D12 * (k(Ix(2, 0), Iy(0, 2)) + k(Ix(0, 2), Iy(2, 0)))
        + 4.0 * D66 * k(Ix(1, 1), Iy(1, 1))
        + 2.0 * D16 * (k(Ix(2, 1), Iy(0, 1)) + k(Ix(1, 2), Iy(1, 0)))
        + 2.0 * D26 * (k(Ix(0, 1), Iy(2, 1)) + k(Ix(1, 0), Iy(1, 2)))
    )
    G = k(Ix(1, 1), Iy(0, 0))
    n = K.shape[0]
    mu = scipy.linalg.eigh(G, K, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0]
    return float(1.0 / mu) if mu > 0 else float("inf")


def panel_buckling_load(
    lam: Laminate,
    Lx_mm: float,
    Ly_mm: float,
    boundary: BoundaryKind = "simply_supported",
) -> float:
    """Full-panel buckling force per unit width N_cr (N/mm) for the intact
    laminate under uniaxial compression along x.

    The Rayleigh-Ritz solution for the simply supported ``Lx_mm`` x
    ``Ly_mm`` rectangle with the full laminate's CLT D matrix, D16 and D26
    included (:func:`_ssss_ritz_N_cr`; it equals the Navier closed form
    when D16 = D26 = 0 and is 3.0% below it on the canonical 200 x 150 mm
    8-ply QI IM7/8552 panel, 12.20 against 12.58 MPa), times the boundary
    factor for clamped (1.9) and free (0.5) panels. Used by the fe3d
    buckling channel after the 3D K_g eigensolve was retired (#129).

    The load direction is fixed at +x (the compression axis); a/b > 1 simply
    means the panel is longer along the load direction than wide.
    """
    if Lx_mm <= 0 or Ly_mm <= 0:
        return float("inf")
    # A typo used to fall back silently to the simply-supported factor.
    require_choice(boundary, _BOUNDARY_BUCKLING_FACTOR, "boundary")
    _, _, D = lam.abd_matrices()
    return _ssss_ritz_N_cr(D, Lx_mm, Ly_mm) * _BOUNDARY_BUCKLING_FACTOR[boundary]


def find_critical_interface(damage: DamageState, lam: Laminate) -> Optional[int]:
    """Return the interface index that would fail first under compression.

    A heuristic score, used by ``fe3d_cai_buckling`` to pick the sublaminate
    for its buckling-onset value; the interface reported as
    ``AnalysisResults.critical_sublaminate`` comes from
    :func:`weakest_sublaminate_growth` instead and can differ.

    Scoring: max_area_i * max(|z_upper_i|, |z_lower_i|), where z is distance
    from interface to the top/bottom laminate surface. Largest wins.

    For non-uniform laminates the ``z_upper`` / ``z_lower`` distances are
    cumulative sums of the actual per-ply thicknesses rather than
    ``n_plies * uniform_thickness``.
    """
    if not damage.delaminations:
        return None
    thicknesses = lam.ply_thicknesses_mm
    total_h = float(sum(thicknesses))
    # cum_z[i] is the through-thickness distance from the bottom face to the
    # top of ply i — i.e. the z position of interface i (interface k separates
    # ply k from ply k+1).
    cum_z = [0.0] * (len(thicknesses) + 1)
    for k, t in enumerate(thicknesses):
        cum_z[k + 1] = cum_z[k] + t
    per_iface_max_area: dict[int, float] = {}
    for e in damage.delaminations:
        per_iface_max_area[e.interface_index] = max(
            per_iface_max_area.get(e.interface_index, 0.0), e.area_mm2
        )

    best_idx: Optional[int] = None
    best_score = -1.0
    for idx, area in per_iface_max_area.items():
        # Distance from the top/bottom laminate surface to the interface.
        z_upper = cum_z[idx + 1]
        z_lower = total_h - cum_z[idx + 1]
        score = area * max(z_upper, z_lower)
        if score > best_score:
            best_score = score
            best_idx = idx
    return best_idx


# --- Post-buckled clamped elliptical film (2D blister) ---------------------
#
# Numerical sizes of the post-buckling Rayleigh-Ritz, not physical constants.
# The in-plane displacement basis has this factor times the bending basis
# degree in each direction (the energy release rate is converged to 0.5% on
# a 3.25:1 ellipse; with the same degree it is up to 7% off). The
# continuation from buckling onset stops at this ratio of film stress to
# film buckling stress; on the exact circular blister the Ritz energy
# release rate is within 0.2% of the exact one up to there.
_MEMBRANE_DEGREE_FACTOR: float = 1.5
_MAX_POSTBUCKLING_RATIO: float = 100.0
_NEWTON_MAX_ITER: int = 50
_NEWTON_RTOL: float = 1e-9
# Continuation steps in N - N_cr: the first to 5% above onset (from
# Koiter's initial post-buckling amplitude), then at most this factor
# (halved while Newton fails).
_FIRST_STEP: float = 0.2
_CONTINUATION_STEP: float = 3.0


@dataclass(frozen=True)
class _FilmBasis:
    """Rayleigh-Ritz fields of the clamped film on the unit disk
    ``(xi, eta) = (x'/a, y'/b)``, at the quadrature points of one half of
    the disk (every integrand is even under ``(xi, eta) -> -(xi, eta)``).

    ``bend`` holds the first and second derivatives (``"x"``, ``"y"``,
    ``"xx"``, ``"yy"``, ``"xy"``) of the deflection trial functions
    ``(1 - xi^2 - eta^2)^2 P_i(xi) P_j(eta)`` with ``i + j`` of one parity,
    ``memb_x`` and ``memb_y`` the first derivatives of the in-plane
    displacement trial functions ``(1 - xi^2 - eta^2) P_i(xi) P_j(eta)``
    with ``i + j`` odd, each of shape ``(n_functions, n_points)``."""

    weights: np.ndarray
    bend: dict[str, np.ndarray]
    memb_x: np.ndarray
    memb_y: np.ndarray


def _ritz_indices(p_x: int, p_y: int) -> list[tuple[int, int]]:
    """Legendre index pairs of the clamped-ellipse basis (as in
    :func:`_clamped_disk_integrals`): ``i / p_x + j / p_y <= 1``."""
    return [
        (i, j) for i in range(p_x + 1) for j in range(p_y + 1) if i * p_y + j * p_x <= p_x * p_y
    ]


@functools.lru_cache(maxsize=64)
def _clamped_film_basis(p_x: int, p_y: int, parity: int) -> _FilmBasis:
    """:class:`_FilmBasis` for the bending basis of
    :func:`_clamped_ellipse_load_factor` (``p_x``, ``p_y``) restricted to
    deflections of the given parity under inversion (0 even, 1 odd), and
    the in-plane basis ``_MEMBRANE_DEGREE_FACTOR`` times its degree. The
    quadrature (Gauss-Legendre in r, the trapezoidal rule in the polar
    angle) is exact for the bending, prestress and membrane stiffness
    matrices; the quartic membrane energy of the deflection is of higher
    degree, and integrating it exactly changes growth stresses by less
    than 1e-8 for twice the points."""
    m_x = math.ceil(_MEMBRANE_DEGREE_FACTOR * p_x)
    m_y = math.ceil(_MEMBRANE_DEGREE_FACTOR * p_y)
    degree = 3 * max(p_x, p_y) + 4
    rg, rw = legendre.leggauss(degree // 2 + 2)
    r = 0.5 * (rg + 1.0)
    n_t = 2 * (degree // 2 + 1)
    t = np.arange(n_t // 2) * (2.0 * np.pi / n_t)
    R, T = np.meshgrid(r, t, indexing="ij")
    xi = (R * np.cos(T)).ravel()
    eta = (R * np.sin(T)).ravel()
    weights = (np.outer(0.5 * rw, np.full(n_t // 2, 4.0 * np.pi / n_t)) * R).ravel()
    g = 1.0 - xi**2 - eta**2

    def tables(
        px: int, py: int, keep: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        idx = [(i, j) for i, j in _ritz_indices(px, py) if (i + j) % 2 == keep]
        Px = _legendre_table(xi, px + 1)
        Py = _legendre_table(eta, py + 1)
        return Px, Py, np.array([i for i, _ in idx]), np.array([j for _, j in idx])

    Px, Py, ii, jj = tables(p_x, p_y, parity)
    m0 = Px[0, ii] * Py[0, jj]
    mx = Px[1, ii] * Py[0, jj]
    my = Px[0, ii] * Py[1, jj]
    G0, Gx, Gy = g * g, -4.0 * g * xi, -4.0 * g * eta
    bend = {
        "x": Gx * m0 + G0 * mx,
        "y": Gy * m0 + G0 * my,
        "xx": (8.0 * xi**2 - 4.0 * g) * m0 + 2.0 * Gx * mx + G0 * Px[2, ii] * Py[0, jj],
        "yy": (8.0 * eta**2 - 4.0 * g) * m0 + 2.0 * Gy * my + G0 * Px[0, ii] * Py[2, jj],
        "xy": 8.0 * xi * eta * m0 + Gx * my + Gy * mx + G0 * Px[1, ii] * Py[1, jj],
    }
    Px, Py, ii, jj = tables(m_x, m_y, 1)
    n0 = Px[0, ii] * Py[0, jj]
    return _FilmBasis(
        weights=weights,
        bend=bend,
        memb_x=-2.0 * xi * n0 + g * Px[1, ii] * Py[0, jj],
        memb_y=-2.0 * eta * n0 + g * Px[0, ii] * Py[1, jj],
    )


class _ContinuationError(RuntimeError):
    """Newton's method lost the post-buckled branch of the film."""


class _PostbuckledFilm:
    """A thin film clamped along an ellipse on a rigid substrate, under a
    uniform in-plane compression ``N (n_xx, n_yy, n_xy)`` (N/mm,
    compression positive), past buckling.

    The film is the thin-film model of Chai, Babcock & Knauss (1981) in two
    dimensions (Chai & Babcock 1985): a von Karman plate with bending
    stiffness ``D`` and extensional stiffness ``A`` (both in the ellipse's
    frame, x' along the semi-axis ``a``), its deflection and in-plane
    displacement zero on the front. The post-buckled state minimises the
    potential energy

        Pi = 1/2 int k.D.k - N/2 int n.(grad w grad w) + 1/2 int e.A.e,
        e = sym grad u + 1/2 grad w grad w

    over a Rayleigh-Ritz basis: the deflection basis of the buckling-onset
    solution (:func:`_clamped_ellipse_load_factor`) restricted to the
    parity of the lowest mode, so the buckling load is the same, and in-plane
    displacements ``(1 - xi^2 - eta^2) P_i P_j``. Newton's method follows
    the branch from onset. The energy release rate of growth of one
    semi-axis with the elliptical shape kept (Chai & Babcock 1985) is
    ``G_a = -(dPi/da) / (pi b)`` (and ``G_b`` likewise), the change of the
    film's energy per unit of new delaminated area at fixed load, from the
    explicit dependence of ``Pi`` on ``a`` at the converged state (the
    state is stationary in every Ritz coefficient).
    """

    def __init__(
        self,
        D: np.ndarray,
        A: np.ndarray,
        a: float,
        b: float,
        n_xx: float,
        n_yy: float,
        n_xy: float,
    ) -> None:
        aspect = max(a, b) / min(a, b)
        p_long = min(_ELLIPSE_MAX_DEGREE, _ELLIPSE_BASE_DEGREE + 4 * math.ceil(aspect - 1.0 - 1e-9))
        p_x, p_y = (p_long, _ELLIPSE_BASE_DEGREE) if a >= b else (_ELLIPSE_BASE_DEGREE, p_long)
        self.a, self.b = a, b
        self.n = (n_xx, n_yy, n_xy)
        self.D, self.A = np.asarray(D, dtype=float), np.asarray(A, dtype=float)
        self.N_cr = math.inf
        best_mu = 0.0
        for parity in (0, 1):
            basis = _clamped_film_basis(p_x, p_y, parity)
            fields = self._fields(basis)
            K, G = self._bending_matrices(fields, basis.weights * a * b)
            m = K.shape[0]
            mu, v = scipy.linalg.eigh(G, K, subset_by_index=[m - 1, m - 1])
            if mu[0] > best_mu:
                best_mu = float(mu[0])
                self._basis, self._f, self._K, self._G = basis, fields, K, G
                self._mode = v[:, 0] / math.sqrt(best_mu)  # mode . G . mode = 1
        if best_mu <= 0.0:
            return
        self.N_cr = 1.0 / best_mu
        self._wa = self._basis.weights * a * b
        self._cx, self._cy = self._basis.memb_x / a, self._basis.memb_y / b
        self._Km = scipy.linalg.cho_factor(self._membrane_matrix())
        self._nm = self._cx.shape[0]
        # Solved states on the branch, (N, coefficients), by increasing N.
        self._states: list[tuple[float, np.ndarray]] = [(self.N_cr, np.zeros_like(self._mode))]
        # Koiter's quartic coefficient of the mode: Pi = 1/2 xi^2 (N_cr - N)
        # + C4 xi^4 for w = xi * mode, which predicts the first state.
        *_, eps = self._membrane(self._mode)
        self._C4 = 0.5 * float(np.sum(eps * (self.A @ eps) * self._wa))

    def _fields(self, basis: _FilmBasis) -> dict[str, np.ndarray]:
        a, b = self.a, self.b
        f = basis.bend
        return {
            "x": f["x"] / a,
            "y": f["y"] / b,
            "xx": f["xx"] / a**2,
            "yy": f["yy"] / b**2,
            "xy": f["xy"] / (a * b),
        }

    def _bending_matrices(
        self, f: dict[str, np.ndarray], wa: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        D = self.D
        n_xx, n_yy, n_xy = self.n

        def J(p: str, q: str) -> np.ndarray:
            return np.asarray((f[p] * wa) @ f[q].T)

        K = (
            D[0, 0] * J("xx", "xx")
            + D[1, 1] * J("yy", "yy")
            + D[0, 1] * (J("xx", "yy") + J("yy", "xx"))
            + 4.0 * D[2, 2] * J("xy", "xy")
            + 2.0 * D[0, 2] * (J("xx", "xy") + J("xy", "xx"))
            + 2.0 * D[1, 2] * (J("yy", "xy") + J("xy", "yy"))
        )
        G = n_xx * J("x", "x") + n_yy * J("y", "y") + n_xy * (J("x", "y") + J("y", "x"))
        return K, G

    def _membrane_matrix(self) -> np.ndarray:
        A, wa, cx, cy = self.A, self._wa, self._cx, self._cy
        xx, yy = (cx * wa) @ cx.T, (cy * wa) @ cy.T
        xy = (cx * wa) @ cy.T
        yx = xy.T
        UU = A[0, 0] * xx + A[2, 2] * yy + A[0, 2] * (xy + yx)
        VV = A[1, 1] * yy + A[2, 2] * xx + A[1, 2] * (yx + xy)
        UV = A[0, 1] * xy + A[2, 2] * yx + A[0, 2] * xx + A[1, 2] * yy
        return np.asarray(np.block([[UU, UV], [UV.T, VV]]))

    def _membrane(self, c: np.ndarray) -> tuple[np.ndarray, ...]:
        """Slopes, in-plane displacement gradients and membrane strain of the
        deflection ``c``, with the in-plane displacement that minimises the
        membrane energy (it is linear in the deflection's strain)."""
        wx, wy = c @ self._f["x"], c @ self._f["y"]
        ew = np.stack([0.5 * wx**2, 0.5 * wy**2, wx * wy])
        s = (self.A @ ew) * self._wa
        rhs = np.concatenate([self._cx @ s[0] + self._cy @ s[2], self._cy @ s[1] + self._cx @ s[2]])
        q = -scipy.linalg.cho_solve(self._Km, rhs)
        qu, qv = q[: self._nm], q[self._nm :]
        ux, uy, vx, vy = qu @ self._cx, qu @ self._cy, qv @ self._cx, qv @ self._cy
        eps = ew + np.stack([ux, vy, uy + vx])
        return wx, wy, ux, uy, vx, vy, eps

    def _newton(self, N: float, c: np.ndarray) -> Optional[np.ndarray]:
        """The post-buckled state at load ``N`` from the guess ``c``; ``None``
        if Newton's method does not converge onto a buckled state."""
        A, wa, fx, fy = self.A, self._wa, self._f["x"], self._f["y"]
        cx, cy = self._cx, self._cy
        K_N = self._K - N * self._G
        lu: Optional[tuple[np.ndarray, np.ndarray]] = None
        norm_prev = math.inf
        for it in range(_NEWTON_MAX_ITER):
            wx, wy, *_, eps = self._membrane(c)
            Nm = (A @ eps) * wa  # membrane force resultant (tension +) x weight
            grad = K_N @ c + fx @ (Nm[0] * wx + Nm[2] * wy) + fy @ (Nm[1] * wy + Nm[2] * wx)
            # Residual small against the bending force it balances (the step
            # test alone can stall at round-off on an ill-conditioned basis).
            norm = float(np.linalg.norm(grad))
            if it > 0 and norm <= _NEWTON_RTOL * np.linalg.norm(self._K @ c):
                break
            # The Hessian is kept while it cuts the residual at least
            # fivefold per step (such a chord step costs a tenth of a Newton
            # step), and rebuilt otherwise.
            if lu is None or norm > 0.2 * norm_prev:
                # Bending and prestress, the membrane force's geometric
                # stiffness, and the membrane energy of the slope increments
                # with the in-plane displacement condensed out.
                E = (fx * wx, fy * wy, fx * wy + fy * wx)
                AE = [A[i, 0] * E[0] + A[i, 1] * E[1] + A[i, 2] * E[2] for i in range(3)]
                L = np.concatenate(
                    [
                        cx @ (AE[0] * wa).T + cy @ (AE[2] * wa).T,
                        cy @ (AE[1] * wa).T + cx @ (AE[2] * wa).T,
                    ]
                )
                H = (
                    K_N
                    + (fx * Nm[0]) @ fx.T
                    + (fy * Nm[1]) @ fy.T
                    + (fx * Nm[2]) @ fy.T
                    + (fy * Nm[2]) @ fx.T
                    + sum((AE[i] * wa) @ E[i].T for i in range(3))
                    - L.T @ scipy.linalg.cho_solve(self._Km, L)
                )
                lu = scipy.linalg.lu_factor(H)
            norm_prev = norm
            dc = scipy.linalg.lu_solve(lu, -grad)
            c = c + dc
            if not np.all(np.isfinite(c)):
                return None
            if np.linalg.norm(dc) <= _NEWTON_RTOL * np.linalg.norm(c):
                break
        else:
            return None
        # On the buckled branch the energy is below the flat film's.
        *_, eps = self._membrane(c)
        energy = 0.5 * c @ K_N @ c + 0.5 * np.sum(eps * (A @ eps) * wa)
        return c if energy < 0.0 else None

    def _state(self, N: float) -> np.ndarray:
        """Ritz coefficients of the post-buckled state at ``N > N_cr``,
        continued from the nearest solved state below it, multiplying
        ``N - N_cr`` by at most ``_CONTINUATION_STEP`` per step (fewer if
        Newton fails)."""
        N_from, c_from = max((s for s in self._states if s[0] <= N), key=lambda s: s[0])
        while N_from < N:
            N_to = min(N, self.N_cr + _CONTINUATION_STEP * (N_from - self.N_cr))
            if N_from <= self.N_cr:
                N_to = min(N, (1.0 + _FIRST_STEP) * self.N_cr)
            while True:
                # Koiter's amplitude from onset; past it, its scaling with
                # sqrt(N - N_cr).
                if N_from <= self.N_cr:
                    guess = math.sqrt((N_to - self.N_cr) / (4.0 * self._C4)) * self._mode
                else:
                    guess = c_from * math.sqrt((N_to - self.N_cr) / (N_from - self.N_cr))
                c = self._newton(N_to, guess)
                if c is not None:
                    break
                # Halve the step: in N - N_cr, geometrically from a buckled
                # state, arithmetically from onset.
                if N_from > self.N_cr:
                    N_mid = self.N_cr + math.sqrt((N_from - self.N_cr) * (N_to - self.N_cr))
                else:
                    N_mid = self.N_cr + 0.5 * (N_to - self.N_cr)
                if (N_mid - N_from) <= 1e-6 * self.N_cr:
                    raise _ContinuationError(
                        "post-buckling continuation of the delaminated sublaminate failed "
                        f"at {N_to / self.N_cr:.4g} times its buckling load"
                    )
                N_to = N_mid
            self._states.append((N_to, c))
            N_from, c_from = N_to, c
        return c_from

    def release_rates(self, N: float) -> tuple[float, float]:
        """``(G_a, G_b)`` (N/mm) at the load ``N`` (N/mm): the energy
        released per unit of new delaminated area when the semi-axis ``a``
        (``b``) grows, the ellipse kept. Zero up to buckling."""
        if not N > self.N_cr:
            return 0.0, 0.0
        _, dPi_da, dPi_db = self._energy(self._state(N), N)
        return -dPi_da / (math.pi * self.b), -dPi_db / (math.pi * self.a)

    def potential_energy(self, N: float) -> float:
        """The post-buckled film's potential energy at the load ``N``
        relative to the flat film (N mm); zero up to buckling."""
        if not N > self.N_cr:
            return 0.0
        return self._energy(self._state(N), N)[0]

    def _energy(self, c: np.ndarray, N: float) -> tuple[float, float, float]:
        """``Pi`` and its derivatives with respect to ``a`` and ``b`` at
        fixed fields in ``(xi, eta)``: there an x' derivative scales as
        ``1/a``, a y' derivative as ``1/b`` and the area as ``a b``. The
        state is stationary in every Ritz coefficient, so these are the
        derivatives of the energy of the post-buckled state."""
        A, D, wa, a, b = self.A, self.D, self._wa, self.a, self.b
        n_xx, n_yy, n_xy = self.n
        wx, wy, ux, uy, vx, vy, eps = self._membrane(c)
        f = self._f
        wxx, wyy, wxy = c @ f["xx"], c @ f["yy"], c @ f["xy"]
        kappa = np.stack([wxx, wyy, 2.0 * wxy])
        M, Nm = D @ kappa, A @ eps
        prestress = n_xx * wx**2 + n_yy * wy**2 + 2.0 * n_xy * wx * wy
        density = (
            0.5 * np.sum(kappa * M, axis=0) - 0.5 * N * prestress + 0.5 * np.sum(eps * Nm, axis=0)
        )
        Pi = float(np.sum(density * wa))
        d_da = (
            -(2.0 * wxx * M[0] + 2.0 * wxy * M[2])
            + N * (n_xx * wx**2 + n_xy * wx * wy)
            - ((ux + wx**2) * Nm[0] + (vx + wx * wy) * Nm[2])
        )
        d_db = (
            -(2.0 * wyy * M[1] + 2.0 * wxy * M[2])
            + N * (n_yy * wy**2 + n_xy * wx * wy)
            - ((vy + wy**2) * Nm[1] + (uy + wx * wy) * Nm[2])
        )
        return Pi, Pi / a + float(np.sum(d_da * wa)) / a, Pi / b + float(np.sum(d_db * wa)) / b

    def growth_load(self, G_c: float, N_cap: float = math.inf) -> tuple[float, int]:
        """The load (N/mm) at which ``max(G_a, G_b)`` reaches ``G_c``, and
        which semi-axis grows (0: ``a``, 1: ``b``). The continuation stops at
        ``_MAX_POSTBUCKLING_RATIO`` times the buckling load, or below it if
        Newton's method loses the branch; above that the energy release
        rate is continued with the straight-sided blister's shape
        ``(N - N_cr)(N + 3 N_cr)``, scaled to match at the last state.
        ``inf`` when the load is above ``N_cap`` (it is then not solved
        for)."""
        if not (math.isfinite(self.N_cr) and N_cap > self.N_cr):
            return math.inf, 0
        r_cap = N_cap / self.N_cr

        def f(r: float) -> tuple[float, int]:
            G = self.release_rates(r * self.N_cr)
            axis = 0 if G[0] >= G[1] else 1
            return math.sqrt(max(G[axis], 0.0)) - math.sqrt(G_c), axis

        def strip_continuation(r_last: float, f_last: float, axis: int) -> tuple[float, int]:
            G_last = (f_last + math.sqrt(G_c)) ** 2
            q = G_c / G_last * (r_last - 1.0) * (r_last + 3.0)
            r = -1.0 + math.sqrt(4.0 + q)
            return (r * self.N_cr if r <= r_cap else math.inf), axis

        lo, f_lo, axis_lo = 1.0, -math.sqrt(G_c), 0
        hi = min(1.0 + _FIRST_STEP, r_cap)
        while True:
            try:
                f_hi, axis = f(hi)
            except _ContinuationError:
                if lo <= 1.0:
                    raise
                return strip_continuation(lo, f_lo, axis_lo)
            if f_hi >= 0.0:
                break
            if hi >= r_cap:
                return math.inf, axis
            lo, f_lo, axis_lo = hi, f_hi, axis
            if hi >= _MAX_POSTBUCKLING_RATIO:
                return strip_continuation(lo, f_lo, axis_lo)
            hi = min(_MAX_POSTBUCKLING_RATIO, r_cap, 1.0 + _CONTINUATION_STEP * (hi - 1.0))
        # Illinois regula falsi on sqrt(G) - sqrt(G_c), bracketed.
        side = 0
        for _ in range(100):
            r = hi - f_hi * (hi - lo) / (f_hi - f_lo)
            f_r, axis_r = f(r)
            if f_r == 0.0 or hi - lo <= 1e-12 * hi:
                hi, axis = r, axis_r
                break
            if f_r > 0.0:
                hi, f_hi, axis = r, f_r, axis_r
                if side == 1:
                    f_lo *= 0.5
                side = 1
            else:
                lo, f_lo = r, f_r
                if side == -1:
                    f_hi *= 0.5
                side = -1
            if abs(f_r) <= 1e-7 * math.sqrt(G_c):
                hi, axis = r, axis_r
                break
        return hi * self.N_cr, axis


def _strip_release_rate(sigma: float, sigma_c: float, E_f: float, h: float) -> float:
    """Energy release rate (N/mm) of the straight-sided (1D) thin-film
    blister with clamped edges (Chai, Babcock & Knauss 1981; Hutchinson &
    Suo 1992): ``h / (2 E_f) (sigma - sigma_c)(sigma + 3 sigma_c)`` at film
    stress ``sigma`` above the buckling stress ``sigma_c`` (MPa), zero
    below. The growth law before the 2D blister; kept as a reference."""
    if sigma <= sigma_c:
        return 0.0
    return h / (2.0 * E_f) * (sigma - sigma_c) * (sigma + 3.0 * sigma_c)


def _strip_growth_film_stress(sigma_c: float, E_f: float, G_c: float, h: float) -> float:
    """Film stress (MPa) at which :func:`_strip_release_rate` reaches
    ``G_c``: ``-sigma_c + sqrt(4 sigma_c^2 + 2 E_f G_c / h)``."""
    return -sigma_c + math.sqrt(4.0 * sigma_c**2 + 2.0 * E_f * G_c / h)


@dataclass(frozen=True)
class _GrowthTerms:
    """What :func:`sublaminate_growth_stress` needs besides the 2D solve:
    the buckling load (N/mm), film thickness (mm), film stress per unit
    laminate stress ``k``, and the film's uniaxial modulus along x (MPa)."""

    interface_index: int
    N_cr: float
    h: float
    k: float
    E_x_film: float


def _growth_terms(
    lam: Laminate, ellipse: DelaminationEllipse, boundary: BoundaryKind
) -> _GrowthTerms:
    i = ellipse.interface_index
    sub_layup, sub_thicknesses = _thinner_sublaminate(lam, i)
    N_cr = sublaminate_buckling_load(lam, ellipse, boundary=boundary)
    h = float(sum(sub_thicknesses))
    if len(sub_layup) == 0 or h <= 0 or not math.isfinite(N_cr) or N_cr <= 0:
        return _GrowthTerms(i, N_cr, h, math.nan, math.nan)
    A_lam, _, _ = lam.abd_matrices()
    a_star = np.linalg.solve(A_lam, np.eye(3)) * lam.thickness_mm  # 1/MPa
    A_sub, _, _ = Laminate(lam.material, sub_layup, sub_thicknesses).abd_matrices()
    k = float(A_sub[0, :] @ a_star[:, 0]) / h
    if not math.isfinite(k):
        k = 0.0
    E_x_film = 1.0 / (float(np.linalg.solve(A_sub, np.eye(3))[0, 0]) * h)
    return _GrowthTerms(i, N_cr, h, k, E_x_film)


def _growth_stress_lower_bound(terms: _GrowthTerms, G_c: float) -> float:
    """A lower bound on the laminate growth stress (MPa) of
    :func:`sublaminate_growth_stress`, used to skip interfaces that cannot
    be the weakest: the straight-sided blister's growth stress at the same
    buckling stress, with the film's uniaxial modulus along x (the softest
    of its moduli). The 2D blister releases less energy at the same film
    stress: stiffening across the load holds its deflection down."""
    if not (math.isfinite(terms.k) and terms.k > 0):
        return math.inf
    sigma_c = terms.N_cr / terms.h
    return _strip_growth_film_stress(sigma_c, terms.E_x_film, G_c, terms.h) / terms.k


def _growth_from_terms(
    lam: Laminate,
    ellipse: DelaminationEllipse,
    terms: _GrowthTerms,
    cap_MPa: float = math.inf,
) -> SublaminateGrowth:
    """:func:`sublaminate_growth_stress` from its :class:`_GrowthTerms`.
    A growth stress above ``cap_MPa`` is not solved for: it is ``inf``."""
    i, N_cr, h, k = terms.interface_index, terms.N_cr, terms.h, terms.k
    if not math.isfinite(k):
        return SublaminateGrowth(i, N_cr, math.inf, math.inf, math.inf)
    if k <= 0:
        return SublaminateGrowth(i, N_cr, N_cr / h, math.inf, math.inf)
    film_terms = _sublaminate_film(lam, ellipse)
    assert film_terms is not None  # N_cr is finite
    D_star, A_sub, a, b, (n_xx, n_yy, n_xy) = film_terms
    film = _PostbuckledFilm(D_star, A_sub, a, b, n_xx, n_yy, n_xy)
    N_g, axis = film.growth_load(mixed_mode_toughness(lam.material), N_cap=cap_MPa * k * h)
    sigma_g = N_g / h
    return SublaminateGrowth(
        interface_index=i,
        buckling_load_N_per_mm=N_cr,
        film_buckling_stress_MPa=N_cr / h,
        film_growth_stress_MPa=sigma_g,
        growth_stress_MPa=sigma_g / k,
        growth_axis=("major", "minor")[axis] if math.isfinite(sigma_g) else None,
    )


def sublaminate_growth_stress(
    lam: Laminate,
    ellipse: DelaminationEllipse,
    boundary: BoundaryKind = "simply_supported",
) -> SublaminateGrowth:
    """Far-field compressive stress at which a delamination grows (MPa).

    The thinner sublaminate over the delamination buckles at the film
    stress ``sigma_c = N_cr / h`` (:func:`sublaminate_buckling_load`), and
    keeps carrying load after buckling. Past buckling it is a 2D blister: a
    thin film clamped along the elliptical front on a rigid substrate
    (Chai, Babcock & Knauss 1981; Chai & Babcock 1985), whose post-buckled
    state comes from a nonlinear (von Karman) Rayleigh-Ritz solution with
    the onset basis (:class:`_PostbuckledFilm`): the reduced bending
    stiffness ``D*`` and the extensional stiffness ``A`` of the
    sublaminate, in the ellipse's frame, under the uniform film stress
    ``sigma`` along x. The delamination grows when the energy released per
    unit of new delaminated area, as one semi-axis of the ellipse grows
    with its shape kept (``G_a`` for ``major_mm``, ``G_b`` for
    ``minor_mm``), reaches the mixed-mode toughness ``G_c``
    (:func:`mixed_mode_toughness`) on either axis. ``G`` is zero at
    buckling and rises with the film stress; for an isotropic circle under
    equibiaxial compression it is Hutchinson & Suo's (1992) circular
    blister, ``G / G_0 = c_1 (1 - (sigma_c/sigma)^2)`` with ``c_1 = 1 /
    (1 + 0.9021 (1 - nu))`` near onset, and the exact axisymmetric
    solution above it. It is below the straight-sided strip's
    ``h / (2 E) (sigma - sigma_c)(sigma + 3 sigma_c)`` at the same
    ``sigma_c``: the film stretches across the load as it deflects.

    Film and laminate share the far-field strain: under a uniaxial
    laminate stress ``sigma_x`` (strains from the laminate's extensional
    compliance) the film x-stress is ``k * sigma_x`` with
    ``k = A_sub[0, :] . a*[:, 0] / h``, and the laminate stress at growth
    is ``sigma_g / k``.

    The thin-film model assumes the sublaminate is thin compared with the
    rest of the laminate; it is least accurate for sublaminates near the
    mid-plane.

    Returns a :class:`SublaminateGrowth` with ``growth_stress_MPa = inf``
    when the sublaminate is degenerate (no plies, zero-area ellipse).
    """
    return _growth_from_terms(lam, ellipse, _growth_terms(lam, ellipse, boundary))


def weakest_sublaminate_growth(
    lam: Laminate,
    damage: DamageState,
    boundary: BoundaryKind = "simply_supported",
) -> Optional[SublaminateGrowth]:
    """The delamination that grows at the lowest far-field stress.

    Evaluates :func:`sublaminate_growth_stress` for the largest ellipse at
    each delaminated interface and returns the minimum (lowest interface
    index on a tie). ``None`` when there are no delaminations. Interfaces
    are taken in order of a lower bound on their growth stress
    (:func:`_growth_stress_lower_bound`), and those whose bound is above
    the lowest growth stress found are not solved; nor is the growth
    stress of one past the lowest found so far.
    """
    largest: dict[int, DelaminationEllipse] = {}
    for e in damage.delaminations:
        current = largest.get(e.interface_index)
        if current is None or e.area_mm2 > current.area_mm2:
            largest[e.interface_index] = e
    if not largest:
        return None
    G_c = mixed_mode_toughness(lam.material)
    terms = {i: _growth_terms(lam, largest[i], boundary) for i in sorted(largest)}
    order = sorted(terms, key=lambda i: (_growth_stress_lower_bound(terms[i], G_c), i))
    best: Optional[SublaminateGrowth] = None
    for i in order:
        if best is not None and _growth_stress_lower_bound(terms[i], G_c) > best.growth_stress_MPa:
            break
        cap = best.growth_stress_MPa if best is not None else math.inf
        g = _growth_from_terms(lam, largest[i], terms[i], cap_MPa=cap)
        if best is None or (g.growth_stress_MPa, i) < (
            best.growth_stress_MPa,
            best.interface_index,
        ):
            best = g
    return best


def semi_analytical_cai(
    lam: Laminate,
    damage: DamageState,
    sigma_pristine_MPa: float,
    A_ref_mm2: float = SOUTIS_REFERENCE_AREA_MM2,
    boundary: BoundaryKind = "simply_supported",
) -> SemiAnalyticalResult:
    """Semi-analytical compression-after-impact residual strength (MPa).

    Takes the minimum of:
      (a) Soutis empirical knockdown at total DPA, and
      (b) the far-field stress at which the weakest delamination grows by
          sublaminate buckling (:func:`weakest_sublaminate_growth`).

    Sublaminate buckling onset alone is not used: a buckled sublaminate
    keeps carrying load until its delamination grows. At the governing
    interface the film buckling stress is typically 20-430 MPa (IM7/8552
    150 x 100 mm panels at 10-30 J).

    Because of (b), the returned residual stress is always less-than-or-equal
    to the empirical-tier result on the same input — so the resulting
    ``knockdown`` (as computed downstream by ``BvidAnalysis.run()``) is
    guaranteed to be less-than-or-equal to the ``empirical`` tier's (with
    the default ``A_ref_mm2``, which is what both tiers use).

    Parameters
    ----------
    lam : Laminate
        The laminate.
    damage : DamageState
        Damage in the panel frame.
    sigma_pristine_MPa : float
        Unnotched compressive strength (MPa) the Soutis knockdown scales.
    A_ref_mm2 : float
        Reference area of the Soutis knockdown, default
        ``SOUTIS_REFERENCE_AREA_MM2`` (15,000 mm², the coupon ``soutis_k_s``
        is calibrated on).
    boundary : BoundaryKind
        Panel edge condition (``"simply_supported"``, ``"clamped"`` or
        ``"free"``); ``BvidAnalysis`` passes ``panel.boundary``. Validated,
        but it does not change the result: the sublaminate is clamped by
        the delamination front whatever the panel's edges (#132; it used
        to scale the sublaminate buckling load by 1.9 or 0.5).

    Returns
    -------
    SemiAnalyticalResult
        Structured result with ``residual_strength_MPa``,
        ``critical_interface_index``, and ``critical_buckling_load_N``.
        If the damage state is empty, ``residual_strength_MPa`` equals
        ``sigma_pristine_MPa`` and the two interface fields are ``None``.

    Examples
    --------
    >>> from bvidfe import MATERIAL_LIBRARY, load_cscan_json
    >>> from bvidfe.core.laminate import Laminate
    >>> lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    >>> r = semi_analytical_cai(lam, load_cscan_json("examples/sample_cscan.json"), 601.3)
    >>> round(r.residual_strength_MPa), r.critical_interface_index
    (266, 5)
    """
    if not damage.delaminations:
        return SemiAnalyticalResult(
            residual_strength_MPa=sigma_pristine_MPa,
            critical_interface_index=None,
            critical_buckling_load_N=None,
        )

    # Soutis bound
    dpa = damage.projected_damage_area_mm2
    sigma_soutis = soutis_cai(lam.material, dpa, A_ref_mm2, sigma_pristine_MPa)

    # Delamination growth bound (damage.delaminations is non-empty here)
    weakest = weakest_sublaminate_growth(lam, damage, boundary=boundary)
    assert weakest is not None
    N_cr = weakest.buckling_load_N_per_mm
    return SemiAnalyticalResult(
        residual_strength_MPa=float(min(sigma_soutis, weakest.growth_stress_MPa)),
        critical_interface_index=weakest.interface_index,
        critical_buckling_load_N=N_cr if math.isfinite(N_cr) else None,
    )


def semi_analytical_tai(
    lam: Laminate,
    damage: DamageState,
    sigma_pristine_MPa: float,
    panel_width_mm: Optional[float] = None,
) -> float:
    """Semi-analytical tension-after-impact residual strength.

    Delegates to the empirical tier's equivalent open hole for the fibre
    damage (:func:`~bvidfe.failure.soutis_openhole.fibre_damage_tai`) on a
    panel ``panel_width_mm`` wide across the load (``None``: infinite
    plate). A cohesive-zone notch model with in-situ ply strength is not
    implemented.
    """
    Kt_inf = lekhnitskii_kt_infinity(lam)
    return fibre_damage_tai(lam.material, damage, sigma_pristine_MPa, Kt_inf, panel_width_mm)
