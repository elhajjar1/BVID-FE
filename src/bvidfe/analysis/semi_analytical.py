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
delamination (:func:`sublaminate_growth_stress`), minimised over the
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
# strip the phase angle moves from about 38 deg at buckling onset towards
# mode II as the film stress rises, so a fixed value is a representative
# mix, not a tracked one. The Benzeggagh-Kenane (1996) exponent is the
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
    """

    interface_index: int
    buckling_load_N_per_mm: float
    film_buckling_stress_MPa: float
    film_growth_stress_MPa: float
    growth_stress_MPa: float


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
    sub_layup, sub_thicknesses = _thinner_sublaminate(lam, ellipse.interface_index)
    if len(sub_layup) == 0:
        return float("inf")

    a = float(ellipse.major_mm)
    b = float(ellipse.minor_mm)
    if not (a > 0 and b > 0 and math.isfinite(a) and math.isfinite(b)):
        return float("inf")

    # Clip the aspect ratio to the domain the Ritz basis is converged on,
    # and surface the clip.
    aspect = max(a, b) / min(a, b)
    if aspect > _MAX_SUBLAMINATE_ASPECT:
        warnings.warn(
            f"Sublaminate ellipse aspect ratio ({aspect:.1f}) exceeds the "
            f"trusted Rayleigh-Ritz domain ({_MAX_SUBLAMINATE_ASPECT:.0f}); "
            f"clipping to {_MAX_SUBLAMINATE_ASPECT:.0f} by widening the shorter "
            f"semi-axis. The sublaminate-buckling channel is degenerate for this "
            f"thin slice and may not constrain the residual strength; the "
            f"empirical Soutis tier likely governs.",
            DegenerateThinSublaminateWarning,
            stacklevel=2,
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
    # The compression along x, in the ellipse's frame: x = (cos t, -sin t).
    t = math.radians(theta)
    c, s = math.cos(t), math.sin(t)
    return _clamped_ellipse_load_factor(D_star, a, b, c * c, s * s, -c * s)


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


def sublaminate_growth_stress(
    lam: Laminate,
    ellipse: DelaminationEllipse,
    boundary: BoundaryKind = "simply_supported",
) -> SublaminateGrowth:
    """Far-field compressive stress at which a delamination grows (MPa).

    The thinner sublaminate over the delamination buckles at the film
    stress ``sigma_c = N_cr / h`` (:func:`sublaminate_buckling_load`), and
    keeps carrying load after buckling. The delamination grows when the
    energy release rate of the buckled film (Chai, Babcock & Knauss 1981;
    Hutchinson & Suo 1992, straight-sided blister)

        G = h / (2 E_f) * (sigma - sigma_c) * (sigma + 3 sigma_c)

    reaches the mixed-mode toughness ``G_c`` (:func:`mixed_mode_toughness`),
    i.e. at the film stress

        sigma_g = -sigma_c + sqrt(4 sigma_c^2 + 2 E_f G_c / h)

    which tends to ``sqrt(2 E_f G_c / h)`` for a film that buckles at a
    negligible stress. Film and laminate share the far-field strain: under
    a uniaxial laminate stress ``sigma_x`` (strains from the laminate's
    extensional compliance) the film x-stress is ``k * sigma_x`` with
    ``k = A_sub[0, :] . a*[:, 0] / h``, ``E_f = k * Ex`` is the film stress
    per unit laminate strain, and the laminate stress at growth is
    ``sigma_g / k``.

    The thin-film energy release rate assumes the sublaminate is thin
    compared with the rest of the laminate; it is least accurate for
    sublaminates near the mid-plane.

    Returns a :class:`SublaminateGrowth` with ``growth_stress_MPa = inf``
    when the sublaminate is degenerate (no plies, zero-area ellipse).
    """
    i = ellipse.interface_index
    sub_layup, sub_thicknesses = _thinner_sublaminate(lam, i)
    N_cr = sublaminate_buckling_load(lam, ellipse, boundary=boundary)
    h = float(sum(sub_thicknesses))
    if len(sub_layup) == 0 or h <= 0 or not math.isfinite(N_cr) or N_cr <= 0:
        return SublaminateGrowth(i, N_cr, math.inf, math.inf, math.inf)

    A_lam, _, _ = lam.abd_matrices()
    a_star = np.linalg.solve(A_lam, np.eye(3)) * lam.thickness_mm  # 1/MPa
    A_sub, _, _ = Laminate(lam.material, sub_layup, sub_thicknesses).abd_matrices()
    k = float(A_sub[0, :] @ a_star[:, 0]) / h
    if not math.isfinite(k) or k <= 0:
        return SublaminateGrowth(i, N_cr, N_cr / h, math.inf, math.inf)
    E_f = k / float(a_star[0, 0])

    sigma_c = N_cr / h
    G_c = mixed_mode_toughness(lam.material)
    sigma_g = -sigma_c + math.sqrt(4.0 * sigma_c**2 + 2.0 * E_f * G_c / h)
    return SublaminateGrowth(
        interface_index=i,
        buckling_load_N_per_mm=N_cr,
        film_buckling_stress_MPa=sigma_c,
        film_growth_stress_MPa=sigma_g,
        growth_stress_MPa=sigma_g / k,
    )


def weakest_sublaminate_growth(
    lam: Laminate,
    damage: DamageState,
    boundary: BoundaryKind = "simply_supported",
) -> Optional[SublaminateGrowth]:
    """The delamination that grows at the lowest far-field stress.

    Evaluates :func:`sublaminate_growth_stress` for the largest ellipse at
    each delaminated interface and returns the minimum (lowest interface
    index on a tie). ``None`` when there are no delaminations.
    """
    largest: dict[int, DelaminationEllipse] = {}
    for e in damage.delaminations:
        current = largest.get(e.interface_index)
        if current is None or e.area_mm2 > current.area_mm2:
            largest[e.interface_index] = e
    if not largest:
        return None
    growths = [sublaminate_growth_stress(lam, largest[i], boundary) for i in sorted(largest)]
    return min(growths, key=lambda g: g.growth_stress_MPa)


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
    interface the film buckling stress is typically 100-430 MPa (IM7/8552
    150 x 100 mm panels at 10-30 J); a one-ply film over the largest,
    back-face delamination buckles at 7-50 MPa.

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
    (266, 4)
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
