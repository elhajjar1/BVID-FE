"""Energy release rate of the post-buckled 2D blister (issue: 2D G law).

The semi-analytical growth stress used to come from the straight-sided (1D)
thin-film blister, ``G = h / (2 E) (sigma - sigma_c)(sigma + 3 sigma_c)``
(Chai, Babcock & Knauss 1981; Hutchinson & Suo 1992). The delamination is an
ellipse, so the film is now a clamped elliptical plate on a rigid substrate,
post-buckled by a nonlinear (von Karman) Rayleigh-Ritz solution
(``_PostbuckledFilm``), and ``G`` is the energy released per unit of new
delaminated area as one semi-axis grows (Chai & Babcock 1985). These tests
check it against the known limits:

- the exact axisymmetric circular blister under equibiaxial compression
  (Evans & Hutchinson 1984), solved here independently as a boundary value
  problem, and Hutchinson & Suo's closed form near onset;
- ``G = 0`` at onset and rising with the stress;
- the straight-sided strip, approached by an ellipse elongated across the
  load and, for any shape, far above onset;
- growth across the load (Chai & Babcock 1985; Whitcomb 1989) when the
  ellipse is not elongated across it;
- ``G`` is the derivative of the post-buckled energy.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_bvp
from scipy.special import j1, jn_zeros

from bvidfe.analysis.semi_analytical import (
    _growth_stress_lower_bound,
    _growth_terms,
    _PostbuckledFilm,
    _strip_release_rate,
    _sublaminate_film,
    mixed_mode_toughness,
    sublaminate_growth_stress,
)
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DelaminationEllipse

_E, _H = 1.0, 0.02  # the von Karman solution depends on sigma / sigma_c only


def _isotropic(nu: float) -> tuple[np.ndarray, np.ndarray]:
    Q = np.array([[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 0.5 * (1.0 - nu)]])
    Q *= _E / (1.0 - nu**2)
    return Q * _H**3 / 12.0, Q * _H


def _exact_circular_blister(nu: float, ratios: list[float]) -> list[float]:
    """``G / G_0`` (``G_0 = (1 - nu) h sigma^2 / E``) of the clamped circular
    film of unit radius under equibiaxial compression ``sigma`` (Evans &
    Hutchinson 1984), from the axisymmetric von Karman equations solved as
    a boundary value problem, and Hutchinson & Suo's (1992) edge formula
    ``G = (1 - nu^2) / (2 E h) (dN^2 + 12 M^2 / h^2)``."""
    D = _E * _H**3 / (12.0 * (1.0 - nu**2))
    C = _E * _H / (1.0 - nu**2)
    k1 = jn_zeros(1, 1)[0]
    sigma_c = D * k1**2 / _H
    r = np.linspace(1e-6, 1.0, 400)
    out, guess = {}, None
    # Continued from just above onset, so that the solver stays on the
    # buckled branch (the flat film solves the problem too).
    for ratio in sorted({1.02, 1.2, *ratios}):
        sigma = ratio * sigma_c
        eps0 = sigma * (1.0 - nu) / _E

        def rhs(x, y):
            u, phi, Nr, Q = y
            du = Nr / C - phi**2 / 2.0 - nu * u / x
            Nth = C * (u / x + nu * (du + phi**2 / 2.0))
            return np.vstack([du, Q / D - phi / x, (Nth - Nr) / x, Nr * phi])

        def bc(ya, yb):
            return np.array([ya[0], ya[1], yb[0] + eps0, yb[1]])

        if guess is None:
            phi = _H * math.sqrt(ratio - 1.0) * k1 * j1(k1 * r)  # the buckling mode
            Q = D * (np.gradient(phi, r) + phi / r)
            y0 = np.vstack([-eps0 * r, phi, -sigma * _H * np.ones_like(r), Q])
        else:
            y0 = guess.sol(r)
        guess = solve_bvp(rhs, bc, r, y0, tol=1e-9, max_nodes=200000)
        assert guess.status == 0
        _, phi1, Nr1, Q1 = guess.sol(1.0)
        M = D * (Q1 / D - phi1) + nu * D * phi1
        G = (1.0 - nu**2) / (2.0 * _E * _H) * ((Nr1 + sigma * _H) ** 2 + 12.0 * M**2 / _H**2)
        out[ratio] = G / ((1.0 - nu) * _H * sigma**2 / _E)
    return [out[ratio] for ratio in ratios]


def test_circular_blister_matches_the_exact_axisymmetric_solution():
    ratios = [1.5, 2.0, 3.0, 5.0, 10.0]
    D, A = _isotropic(0.3)
    film = _PostbuckledFilm(D, A, 1.0, 1.0, 1.0, 1.0, 0.0)
    exact = _exact_circular_blister(0.3, ratios)
    for ratio, G_exact in zip(ratios, exact):
        sigma = ratio * film.N_cr / _H
        G_a, G_b = film.release_rates(ratio * film.N_cr)
        assert G_a == pytest.approx(G_b, rel=1e-9)  # the circle grows evenly
        assert G_a / ((1.0 - 0.3) * _H * sigma**2 / _E) == pytest.approx(G_exact, rel=2e-4)
    # The exact values at 2 and 5 times onset, to the digits printed.
    assert exact[1] == pytest.approx(0.48087, abs=2e-5)
    assert exact[3] == pytest.approx(0.67253, abs=2e-5)


@pytest.mark.parametrize("nu", [0.2, 0.3, 0.4])
def test_circular_blister_near_onset_is_hutchinson_suo(nu):
    """``G / G_0 = c_1 (1 - (sigma_c / sigma)^2)``, ``c_1 = 1 / (1 + 0.9021
    (1 - nu))`` (Hutchinson & Suo 1992), and the buckling stress ``1.2235 E
    / (1 - nu^2) (h / R)^2``."""
    D, A = _isotropic(nu)
    film = _PostbuckledFilm(D, A, 1.0, 1.0, 1.0, 1.0, 0.0)
    sigma_c = film.N_cr / _H
    assert sigma_c == pytest.approx(1.2235 * _E / (1.0 - nu**2) * _H**2, rel=1e-4)
    c1 = 1.0 / (1.0 + 0.9021 * (1.0 - nu))
    ratio = 1.002
    G_a, _ = film.release_rates(ratio * film.N_cr)
    G_0 = (1.0 - nu) * _H * (ratio * sigma_c) ** 2 / _E
    assert G_a / G_0 == pytest.approx(c1 * (1.0 - ratio**-2), rel=1e-3)


def test_release_rate_is_zero_at_onset_and_rises_with_the_load():
    D, A = _isotropic(0.3)
    film = _PostbuckledFilm(D, A, 1.0, 2.0, 0.75, 0.25, -math.sqrt(3.0) / 4.0)
    assert film.release_rates(film.N_cr) == (0.0, 0.0)
    assert film.release_rates(0.5 * film.N_cr) == (0.0, 0.0)
    G = [max(film.release_rates(r * film.N_cr)) for r in (1.01, 1.1, 1.5, 2.0, 4.0, 8.0)]
    assert G[0] > 0.0
    assert all(g1 > g0 for g0, g1 in zip(G, G[1:]))


def _laminate_film(theta=30.0):
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    ellipse = DelaminationEllipse(1, (0.0, 0.0), 14.0, 8.0, theta)
    film = _sublaminate_film(lam, ellipse)
    assert film is not None
    D, A, a, b, n = film
    return D, A, a, b, n


def test_release_rates_are_the_derivatives_of_the_postbuckled_energy():
    """``G_a = -(dPi/da) / (pi b)`` from re-solved films at ``a (1 +- d)``,
    on an unsymmetric [0/45] sublaminate (D16, D26) at 30 deg."""
    D, A, a, b, n = _laminate_film()
    film = _PostbuckledFilm(D, A, a, b, *n)
    N = 2.5 * film.N_cr
    G_a, G_b = film.release_rates(N)
    d = 1e-3
    dPi_da = (
        _PostbuckledFilm(D, A, a * (1 + d), b, *n).potential_energy(N)
        - _PostbuckledFilm(D, A, a * (1 - d), b, *n).potential_energy(N)
    ) / (2 * d * a)
    dPi_db = (
        _PostbuckledFilm(D, A, a, b * (1 + d), *n).potential_energy(N)
        - _PostbuckledFilm(D, A, a, b * (1 - d), *n).potential_energy(N)
    ) / (2 * d * b)
    assert G_a == pytest.approx(-dPi_da / (math.pi * b), rel=1e-4)
    assert G_b == pytest.approx(-dPi_db / (math.pi * a), rel=1e-4)
    assert film.potential_energy(N) < 0.0


def test_buckling_load_is_the_onset_one():
    """The post-buckling basis has the onset basis's buckling load
    (:func:`sublaminate_buckling_load`)."""
    from bvidfe.analysis.semi_analytical import _clamped_ellipse_load_factor

    D, A, a, b, n = _laminate_film()
    assert _PostbuckledFilm(D, A, a, b, *n).N_cr == pytest.approx(
        _clamped_ellipse_load_factor(D, a, b, *n), rel=1e-9
    )


def _ratio_to_strip(b_over_a: float, ratio: float) -> tuple[float, float]:
    """``(G_a, G_b) / G_strip`` for an isotropic ellipse with ``a`` along the
    load, at ``ratio`` times its buckling stress; the strip has the same
    buckling stress and the plane-strain modulus."""
    nu = 0.3
    D, A = _isotropic(nu)
    film = _PostbuckledFilm(D, A, 1.0, b_over_a, 1.0, 0.0, 0.0)
    sigma_c = film.N_cr / _H
    G_strip = _strip_release_rate(ratio * sigma_c, sigma_c, _E / (1.0 - nu**2), _H)
    G_a, G_b = film.release_rates(ratio * film.N_cr)
    return G_a / G_strip, G_b / G_strip


def test_2d_blister_releases_less_than_the_strip_and_more_as_it_elongates_across_the_load():
    """At three times onset the circle under uniaxial load releases 0.47 of
    the strip's energy (growing across the load); elongated across the
    load, growth along the load approaches the strip: 0.50 at 2:1, 0.72
    at 3.25:1."""
    circle = _ratio_to_strip(1.0, 3.0)
    assert circle[1] > circle[0]  # across the load
    assert max(circle) == pytest.approx(0.470, abs=0.005)
    two = _ratio_to_strip(2.0, 3.0)
    peanut = _ratio_to_strip(3.25, 3.0)
    assert two[0] == pytest.approx(0.501, abs=0.005)
    assert peanut[0] == pytest.approx(0.725, abs=0.005)
    assert max(circle) < 1.0 and two[0] < peanut[0] < 1.0


@pytest.mark.slow
def test_elongated_ellipse_across_the_load_tends_to_the_strip():
    """At 10:1 (the aspect-ratio clip) the ellipse's buckling stress is
    within 5% of the clamped strip's, and at ten times onset growth along
    the load releases the strip's energy to 1%."""
    G_a, _ = _ratio_to_strip(10.0, 10.0)
    assert G_a == pytest.approx(1.0, abs=0.01)


def test_far_above_onset_the_blister_releases_the_strips_energy():
    """For any shape, far above onset the energy release rate of the growing
    axis tends to the straight-sided strip's steady state: the circle at
    100 times onset releases 0.97 of it (across the load)."""
    G = _ratio_to_strip(1.0, 100.0)
    assert G[1] == pytest.approx(0.976, abs=0.01)


def test_growth_is_across_the_load_unless_the_ellipse_is_elongated_across_it():
    """Chai & Babcock (1985) and Whitcomb (1989): a buckled thin-film
    delamination under uniaxial compression grows across the load."""
    for b_over_a in (1.0, 0.5, 1 / 3.25):  # circle, and elongated along the load
        G_a, G_b = _ratio_to_strip(b_over_a, 2.0)
        assert G_b > G_a
    G_a, G_b = _ratio_to_strip(3.25, 2.0)
    assert G_a > G_b


_M = MATERIAL_LIBRARY["IM7/8552"]
_QI16 = [0, 45, -45, 90] * 4


@pytest.mark.parametrize("interface", [0, 1, 3, 7, 14])
def test_growth_stress_is_above_onset_and_its_lower_bound(interface):
    lam = Laminate(_M, _QI16, 0.152)
    ellipse = DelaminationEllipse(interface, (0.0, 0.0), 18.0, 9.0, 45.0)
    g = sublaminate_growth_stress(lam, ellipse)
    assert g.film_growth_stress_MPa > g.film_buckling_stress_MPa
    terms = _growth_terms(lam, ellipse, "simply_supported")
    assert g.growth_stress_MPa >= _growth_stress_lower_bound(terms, mixed_mode_toughness(_M))
    assert g.growth_axis in ("major", "minor")


def test_growth_releases_exactly_the_mixed_mode_toughness():
    lam = Laminate(_M, _QI16, 0.152)
    ellipse = DelaminationEllipse(3, (0.0, 0.0), 20.0, 15.0, 0.0)
    g = sublaminate_growth_stress(lam, ellipse)
    film_terms = _sublaminate_film(lam, ellipse)
    assert film_terms is not None
    D, A, a, b, n = film_terms
    film = _PostbuckledFilm(D, A, a, b, *n)
    h = 4 * 0.152
    G = film.release_rates(g.film_growth_stress_MPa * h)
    assert max(G) == pytest.approx(mixed_mode_toughness(_M), rel=1e-6)
    assert g.growth_axis == ("major", "minor")[int(G[1] > G[0])]
