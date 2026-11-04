"""Stress field around a circular hole in an infinite orthotropic plate.

Lekhnitskii's complex-potential solution (Lekhnitskii, *Anisotropic
Plates*, 1968; Tan, *Stress Concentrations in Laminated Composites*, 1994)
for a traction-free circular hole of radius R in an infinite plate of an
orthotropic membrane (``a16 = a26 = 0``) under a uniform far-field
``(sigma_x, sigma_y, tau_xy)``. With the compliance ``a = h * A^-1``
(``a11, a22, a12, a66``) the characteristic equation

    a11 mu^4 + (2 a12 + a66) mu^2 + a22 = 0

has two roots ``mu_1, mu_2`` with positive imaginary part. For an in-plane
isotropic laminate (quasi-isotropic, ``mu_1 = mu_2 = i``) the general
solution is singular and Kirsch's solution is used instead.

The axial (x) stress on the net section (0, y) at ``y = R`` is the
stress concentration ``1 + n - k * beta`` per unit far-field sigma_x, with
``beta = sigma_y / sigma_x``, ``k = sqrt(Ex / Ey)`` and ``n =
sqrt(2 (k - nu_xy) + Ex / Gxy)``: at ``beta = 0`` the
``lekhnitskii_kt_infinity`` of :mod:`bvidfe.failure.soutis_openhole`.

:func:`biaxial_point_stress_factor` is the biaxial correction of the
empirical TAI hole (issue #115): the ratio of the axial point stress at
``(0, R + d0)`` under uniaxial load to that under the biaxial load, the
factor by which a lateral stress changes the Whitney-Nuismer point-stress
knockdown. Like ``lekhnitskii_kt_infinity`` it uses the laminate's
``effective_engineering_constants``, so any ``A16`` / ``A26`` coupling of
an unbalanced layup is neglected.
"""

from __future__ import annotations

import cmath
import math
from typing import Tuple

import numpy as np

from bvidfe.core.laminate import Laminate

#: Relative separation of the two characteristic roots below which the
#: laminate is taken as in-plane isotropic and Kirsch's solution is used.
#: A computed quasi-isotropic laminate has roots ``i`` to round-off; the
#: general solution divides by ``mu_1 - mu_2``.
REPEATED_ROOT_TOLERANCE = 1e-6

Compliance = Tuple[float, float, float, float]


def laminate_compliance(lam: Laminate) -> Compliance:
    """``(a11, a22, a12, a66)`` of ``h * A^-1`` from the laminate's
    effective engineering constants: ``(1/Ex, 1/Ey, -nu_xy/Ex, 1/Gxy)``."""
    Ex, Ey, Gxy, nu_xy = lam.effective_engineering_constants()
    return (1.0 / Ex, 1.0 / Ey, -nu_xy / Ex, 1.0 / Gxy)


def characteristic_roots(compliance: Compliance) -> Tuple[complex, complex]:
    """The two roots ``mu`` (positive imaginary part, ``mu_1`` the smaller)
    of ``a11 mu^4 + (2 a12 + a66) mu^2 + a22 = 0``.

    Raises ``ValueError`` for a compliance that has no such pair (not
    positive definite)."""
    a11, a22, a12, a66 = (float(c) for c in compliance)
    if not (a11 > 0 and a22 > 0 and a66 > 0 and a12 * a12 < a11 * a22):
        raise ValueError(
            f"compliance (a11, a22, a12, a66) = {compliance!r} is not positive definite"
        )
    # mu^2 = t solves a11 t^2 + (2 a12 + a66) t + a22 = 0.
    b = 2.0 * a12 + a66
    disc = cmath.sqrt(b * b - 4.0 * a11 * a22)
    t1 = (-b + disc) / (2.0 * a11)
    t2 = (-b - disc) / (2.0 * a11)
    roots = []
    for t in (t1, t2):
        mu = cmath.sqrt(t)
        roots.append(mu if mu.imag > 0 else -mu)
    roots.sort(key=lambda z: (z.imag, z.real))
    return roots[0], roots[1]


def _is_isotropic(compliance: Compliance, mu: Tuple[complex, complex]) -> bool:
    """True when the two roots coincide (to ``REPEATED_ROOT_TOLERANCE``).

    Repeated roots with ``a11 != a22`` are not in-plane isotropic and have
    no closed form here: ``ValueError``."""
    if abs(mu[0] - mu[1]) >= REPEATED_ROOT_TOLERANCE * abs(mu[0]):
        return False
    a11, a22, _, _ = compliance
    if abs(a11 - a22) > 1e-4 * a11:
        raise ValueError(
            "the orthotropic hole solution has repeated characteristic roots for this "
            f"compliance {compliance!r} but the laminate is not in-plane isotropic "
            "(a11 != a22); this degenerate case is not supported"
        )
    return True


def _kirsch_field(
    far: Tuple[float, float, float], R: float, x: np.ndarray, y: np.ndarray
) -> np.ndarray:
    """Kirsch's isotropic solution in Cartesian components, ``(3, ...)``."""
    sx, sy, txy = far
    r2 = x * x + y * y
    th = np.arctan2(y, x)
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    q2 = R * R / r2
    q4 = q2 * q2
    p = 0.5 * (sx + sy)  # mean far-field stress
    d = 0.5 * (sx - sy)  # deviatoric part along x
    srr = p * (1 - q2) + (d * c2 + txy * s2) * (1 - 4 * q2 + 3 * q4)
    stt = p * (1 + q2) - (d * c2 + txy * s2) * (1 + 3 * q4)
    srt = (-d * s2 + txy * c2) * (1 + 2 * q2 - 3 * q4)
    c, s = np.cos(th), np.sin(th)
    out = np.empty((3,) + np.shape(x))
    out[0] = srr * c * c + stt * s * s - 2 * srt * s * c
    out[1] = srr * s * s + stt * c * c + 2 * srt * s * c
    out[2] = (srr - stt) * s * c + srt * (c * c - s * s)
    return out


def _lekhnitskii_field(
    mu: Tuple[complex, complex],
    far: Tuple[float, float, float],
    R: float,
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Lekhnitskii's complex-potential solution, ``(3, ...)``.

    ``Phi_k(z_k) = G_k z_k + A_k / zeta_k`` with ``z_k = x + mu_k y`` and
    the conformal map ``z_k = p_k zeta_k + q_k / zeta_k``, ``p_k = R (1 - i
    mu_k) / 2``, ``q_k = R (1 + i mu_k) / 2`` (``|zeta_k| = 1`` on the
    hole). ``G_k`` give the uniform far field; ``A_k`` cancel its tractions
    on the hole."""
    sx, sy, txy = far
    m = np.array(mu, dtype=complex)
    # 2 Re sum mu^2 G = sx, 2 Re sum G = sy, -2 Re sum mu G = txy, Im G_1 = 0
    # (the free rigid rotation): unknowns Re G_1, Re G_2, Im G_2.
    M = np.array(
        [
            [2 * (m[0] ** 2).real, 2 * (m[1] ** 2).real, -2 * (m[1] ** 2).imag],
            [2.0, 2.0, 0.0],
            [-2 * m[0].real, -2 * m[1].real, 2 * m[1].imag],
        ]
    )
    g = np.linalg.solve(M, np.array([sx, sy, txy], dtype=float))
    G = np.array([g[0], g[1] + 1j * g[2]])
    p = R * (1 - 1j * m) / 2
    q = R * (1 + 1j * m) / 2
    # Traction-free hole: sum (A_k + conj ...) cancels the uniform part.
    r1 = -np.sum(np.conj(G) * np.conj(p) + G * q)
    r2 = -np.sum(np.conj(m * G) * np.conj(p) + m * G * q)
    A = np.linalg.solve(np.array([[1, 1], [m[0], m[1]]]), np.array([r1, r2]))
    dphi = []
    for k in range(2):
        z = x + m[k] * y
        root = np.sqrt(z * z - 4 * p[k] * q[k])
        zeta_a = (z + root) / (2 * p[k])
        zeta_b = (z - root) / (2 * p[k])
        zeta = np.where(np.abs(zeta_a) >= np.abs(zeta_b), zeta_a, zeta_b)  # outside the hole
        dz_dzeta = p[k] - q[k] / zeta**2
        dphi.append(G[k] - A[k] / (zeta**2 * dz_dzeta))
    out = np.empty((3,) + np.shape(x))
    out[0] = 2 * np.real(m[0] ** 2 * dphi[0] + m[1] ** 2 * dphi[1])
    out[1] = 2 * np.real(dphi[0] + dphi[1])
    out[2] = -2 * np.real(m[0] * dphi[0] + m[1] * dphi[1])
    return out


def hole_stress_field(
    compliance: Compliance,
    far_field_MPa: Tuple[float, float, float],
    radius_mm: float,
    x_mm: np.ndarray,
    y_mm: np.ndarray,
) -> np.ndarray:
    """``(3, ...)`` array of ``sigma_x, sigma_y, tau_xy`` (MPa) at the points
    ``(x_mm, y_mm)`` (hole centre at the origin, on or outside the hole)
    around a traction-free circular hole of radius ``radius_mm`` in an
    infinite orthotropic plate with compliance ``(a11, a22, a12, a66)``
    (``a16 = a26 = 0``) under the uniform far field ``far_field_MPa =
    (sigma_x, sigma_y, tau_xy)``.

    Lekhnitskii's complex potentials; repeated roots (in-plane isotropic,
    e.g. quasi-isotropic laminates) use the Kirsch solution.

    >>> import numpy as np
    >>> iso = (1.0, 1.0, -0.3, 2.6)  # isotropic, nu = 0.3
    >>> s = hole_stress_field(iso, (1.0, 0.0, 0.0), 1.0, np.array(0.0), np.array(1.0))
    >>> round(float(s[0]), 6)  # Kirsch: 3 at the hole edge on the net section
    3.0
    """
    if not radius_mm > 0:
        raise ValueError(f"radius_mm must be > 0 (got {radius_mm!r})")
    x = np.asarray(x_mm, dtype=float)
    y = np.asarray(y_mm, dtype=float)
    far = tuple(float(f) for f in far_field_MPa)
    if len(far) != 3:
        raise ValueError(
            f"far_field_MPa must be (sigma_x, sigma_y, tau_xy) (got {far_field_MPa!r})"
        )
    mu = characteristic_roots(compliance)
    if _is_isotropic(compliance, mu):
        return _kirsch_field(far, float(radius_mm), x, y)  # type: ignore[arg-type]
    return _lekhnitskii_field(mu, far, float(radius_mm), x, y)  # type: ignore[arg-type]


def biaxial_point_stress_factor(
    lam: Laminate, radius_mm: float, d0_mm: float, transverse_stress_ratio: float
) -> float:
    """``F(beta) = sigma_x(0, R + d0; 1, 0) / sigma_x(0, R + d0; 1, beta)``:
    the axial stress at the Whitney-Nuismer characteristic point on the net
    section under a unit uniaxial far field over that under the far field
    ``(1, beta, 0)``, ``beta = transverse_stress_ratio``.

    The point-stress strength of the hole scales by ``F``: ``F < 1`` where
    the lateral stress raises the axial stress at the point. Exactly
    ``1.0`` at ``beta == 0`` or for no hole (``radius_mm <= 0``). At the
    BVID TAI hole sizes (R = 0.5-2 mm, d0 = 1 mm) F is within about 5% of 1
    for ``|beta| <= 1``; its sign of effect flips near ``R = d0 / (sqrt(3)
    - 1)`` for a quasi-isotropic laminate (Kirsch: the sigma_y-induced
    sigma_x at r is ``sigma_y (R^2/r^2 - 3 R^4/r^4) / 2``).

    Raises ``ValueError`` when the axial stress at the point is not
    positive under the biaxial far field (no axial point-stress failure is
    defined); this needs ``sqrt(Ex/Ey) * beta`` well above 1 and does not
    occur for the presets at ``|beta| <= 1``.

    >>> from bvidfe import MATERIAL_LIBRARY
    >>> from bvidfe.core.laminate import Laminate
    >>> qi = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    >>> biaxial_point_stress_factor(qi, 0.75, 1.0, 0.0)
    1.0
    >>> round(biaxial_point_stress_factor(qi, 0.75, 1.0, 1.0), 3)
    0.965
    """
    beta = float(transverse_stress_ratio)
    if beta == 0 or radius_mm <= 0:
        return 1.0
    a = laminate_compliance(lam)
    x, y = np.array(0.0), np.array(radius_mm + d0_mm)
    gx = float(hole_stress_field(a, (1.0, 0.0, 0.0), radius_mm, x, y)[0])
    gy = float(hole_stress_field(a, (0.0, 1.0, 0.0), radius_mm, x, y)[0])
    denom = gx + beta * gy
    if not denom > 0 or not math.isfinite(denom):
        raise ValueError(
            f"the axial stress at the characteristic point (0, R + d0) is not tensile "
            f"under sigma_y/sigma_x = {beta:g} (sigma_x = {denom:.3g} per unit far-field "
            f"sigma_x), so the axial point-stress criterion is undefined"
        )
    return gx / denom
