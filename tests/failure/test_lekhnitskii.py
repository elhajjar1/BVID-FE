"""Orthotropic open-hole stress field and the biaxial TAI point-stress factor
(``bvidfe.failure.lekhnitskii``, issue #115)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.lekhnitskii import (
    biaxial_point_stress_factor,
    characteristic_roots,
    hole_stress_field,
    laminate_compliance,
)
from bvidfe.failure.soutis_openhole import _point_stress_knockdown, lekhnitskii_kt_infinity

_IM7 = MATERIAL_LIBRARY["IM7/8552"]
QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
CROSS = [0, 90, 0, 90, 90, 0, 90, 0]
UD = [0] * 8
LAYUPS = {
    "QI8": QI8,
    "cross-ply": CROSS,
    "[0]8": UD,
    "hard": [0, 45, 0, -45, 90, 90, -45, 0, 45, 0],
    "[90]8": [90] * 8,
    "[+-45]2s": [45, -45, 45, -45, -45, 45, -45, 45],
}
_ISO = (1.0 / 70000.0, 1.0 / 70000.0, -0.3 / 70000.0, 2.6 / 70000.0)  # aluminium-like


def _lam(layup):
    return Laminate(_IM7, list(layup), 0.152)


def _kirsch_polar(far, R, r, th):
    """Kirsch in polar components for a far field (sx, sy, txy)."""
    sx, sy, txy = far
    p, d = (sx + sy) / 2, (sx - sy) / 2
    q2, q4 = R * R / r**2, R**4 / r**4
    w = d * np.cos(2 * th) + txy * np.sin(2 * th)
    srr = p * (1 - q2) + w * (1 - 4 * q2 + 3 * q4)
    stt = p * (1 + q2) - w * (1 + 3 * q4)
    srt = (-d * np.sin(2 * th) + txy * np.cos(2 * th)) * (1 + 2 * q2 - 3 * q4)
    return srr, stt, srt


def _to_polar(s, th):
    sx, sy, txy = s
    c, si = np.cos(th), np.sin(th)
    srr = sx * c * c + sy * si * si + 2 * txy * si * c
    stt = sx * si * si + sy * c * c - 2 * txy * si * c
    srt = (sy - sx) * si * c + txy * (c * c - si * si)
    return srr, stt, srt


@pytest.mark.parametrize("layup", sorted(LAYUPS))
def test_uniaxial_load_reproduces_lekhnitskii_kt_infinity(layup):
    lam = _lam(LAYUPS[layup])
    s = hole_stress_field(laminate_compliance(lam), (1.0, 0.0, 0.0), 2.0, 0.0, 2.0)
    assert float(s[0]) == pytest.approx(lekhnitskii_kt_infinity(lam), rel=1e-9)


@pytest.mark.parametrize("layup", ["cross-ply", "[0]8", "hard", "[+-45]2s"])
@pytest.mark.parametrize("beta", [-1.0, -0.5, 0.5, 1.0])
def test_biaxial_concentration_at_the_net_section_edge_is_1_plus_n_minus_k_beta(layup, beta):
    lam = _lam(LAYUPS[layup])
    Ex, Ey, _, _ = lam.effective_engineering_constants()
    k = math.sqrt(Ex / Ey)
    s = hole_stress_field(laminate_compliance(lam), (1.0, beta, 0.0), 1.0, 0.0, 1.0)
    assert float(s[0]) == pytest.approx(lekhnitskii_kt_infinity(lam) - k * beta, rel=1e-9)


@pytest.mark.parametrize(
    "layup,expected",
    [("QI8", (4.0, 3.5, 3.0, 2.5, 2.0)), ("[0]8", (11.568, 9.352, 7.136, 4.920, 2.704))],
)
def test_probe_2_3_table(layup, expected):
    a = laminate_compliance(_lam(LAYUPS[layup]))
    got = [
        float(hole_stress_field(a, (1.0, b, 0.0), 1.0, 0.0, 1.0)[0])
        for b in (-1.0, -0.5, 0.0, 0.5, 1.0)
    ]
    assert got == pytest.approx(list(expected), abs=5e-4)


@pytest.mark.parametrize("far", [(1.0, 0.0, 0.0), (0.3, -0.8, 0.0), (0.0, 0.0, 1.0), (1, 1, 0.5)])
def test_quasi_isotropic_laminate_is_kirsch(far):
    a = laminate_compliance(_lam(QI8))
    mu = characteristic_roots(a)
    assert abs(mu[0] - mu[1]) < 1e-6  # repeated roots i, i
    r = np.array([1.0, 1.3, 2.0, 5.0])[:, None]
    th = np.linspace(0.0, 2 * np.pi, 13)[None, :]
    R = 1.0
    s = hole_stress_field(a, far, R, r * np.cos(th), r * np.sin(th))
    got = _to_polar(s, th)
    exp = _kirsch_polar(far, R, r, th)
    for g, e in zip(got, exp):
        np.testing.assert_allclose(g, e, atol=1e-12)


def test_near_isotropic_general_solution_converges_to_kirsch():
    """The Lekhnitskii branch (roots apart) approaches the Kirsch branch."""
    a11, a22, a12, a66 = _ISO
    near = (a11, a22 * (1 + 1e-4), a12, a66)
    mu = characteristic_roots(near)
    assert abs(mu[0] - mu[1]) > 1e-6 * abs(mu[0])  # general branch
    x = np.array([0.0, 1.5, 0.7, -2.0])
    y = np.array([1.0, 0.0, 1.2, 0.4])
    for far in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]:
        np.testing.assert_allclose(
            hole_stress_field(near, far, 1.0, x, y),
            hole_stress_field(_ISO, far, 1.0, x, y),
            atol=2e-3,
        )


@pytest.mark.parametrize("layup", sorted(LAYUPS))
@pytest.mark.parametrize("far", [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1, -0.6, 0.3)])
def test_the_hole_is_traction_free(layup, far):
    a = laminate_compliance(_lam(LAYUPS[layup]))
    th = np.linspace(0.0, 2 * np.pi, 73)
    R = 1.7
    sx, sy, txy = hole_stress_field(a, far, R, R * np.cos(th), R * np.sin(th))
    nx, ny = np.cos(th), np.sin(th)
    tx = sx * nx + txy * ny
    ty = txy * nx + sy * ny
    assert np.max(np.abs(tx)) < 1e-9
    assert np.max(np.abs(ty)) < 1e-9


@pytest.mark.parametrize("layup", sorted(LAYUPS))
@pytest.mark.parametrize("far", [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1, -0.6, 0.3)])
def test_the_far_field_is_recovered_at_100_radii(layup, far):
    """The hole's disturbance decays as (R / r)^2: a few 1e-3 at 100 R (up
    to 3e-3 for [0]8) and 100 times less at 1000 R."""
    a = laminate_compliance(_lam(LAYUPS[layup]))
    th = np.linspace(0.0, 2 * np.pi, 37)

    def err(r):
        s = hole_stress_field(a, far, 1.0, r * np.cos(th), r * np.sin(th))
        return max(float(np.max(np.abs(s[k] - far[k]))) for k in range(3))

    assert err(100.0) < 5e-3
    assert err(1000.0) == pytest.approx(err(100.0) / 100, rel=0.05)


def test_superposition_is_linear():
    a = laminate_compliance(_lam(CROSS))
    x, y = np.array([0.3, 0.0]), np.array([1.4, 2.0])
    combined = hole_stress_field(a, (2.0, -0.5, 0.7), 1.0, x, y)
    parts = (
        2.0 * hole_stress_field(a, (1.0, 0.0, 0.0), 1.0, x, y)
        - 0.5 * hole_stress_field(a, (0.0, 1.0, 0.0), 1.0, x, y)
        + 0.7 * hole_stress_field(a, (0.0, 0.0, 1.0), 1.0, x, y)
    )
    np.testing.assert_allclose(combined, parts, atol=1e-12)


@pytest.mark.parametrize("layup", sorted(LAYUPS))
def test_factor_is_exactly_one_without_a_lateral_stress(layup):
    lam = _lam(LAYUPS[layup])
    assert biaxial_point_stress_factor(lam, 0.75, 1.0, 0.0) == 1.0
    assert biaxial_point_stress_factor(lam, 0.75, 1.0, -0.0) == 1.0
    assert biaxial_point_stress_factor(lam, 0.0, 1.0, 0.7) == 1.0  # no hole


@pytest.mark.parametrize("R", [0.5, 0.75, 2.0, 20.0])
def test_exact_point_stress_is_konish_whitney_for_quasi_isotropic(R):
    """For QI the Konish-Whitney approximation behind the TAI knockdown is the
    exact (Kirsch) axial stress at (0, R + d0)."""
    lam = _lam(QI8)
    sx = float(hole_stress_field(laminate_compliance(lam), (1.0, 0.0, 0.0), R, 0.0, R + 1.0)[0])
    assert 1.0 / sx == pytest.approx(
        _point_stress_knockdown(R, 1.0, lekhnitskii_kt_infinity(lam)), rel=1e-9
    )


# Probe 2.4 of the #115 design (d0 = 1 mm): F at beta = -1, -0.5, +0.5, +1.
_PROBE_F = {
    ("QI8", 0.5): (1.036, 1.018, 0.983, 0.967),
    ("QI8", 0.75): (1.037, 1.018, 0.982, 0.965),
    ("QI8", 2.0): (0.953, 0.976, 1.025, 1.051),
    ("QI8", 20.0): (0.775, 0.873, 1.170, 1.409),
    ("cross-ply", 0.75): (1.007, 1.004, 0.996, 0.993),
    ("cross-ply", 2.0): (0.964, 0.982, 1.019, 1.039),
    ("[0]8", 0.75): (1.061, 1.029, 0.972, 0.946),
}


@pytest.mark.parametrize("key", sorted(_PROBE_F))
def test_factor_matches_the_probe_table(key):
    layup, R = key
    lam = _lam(LAYUPS[layup])
    got = [biaxial_point_stress_factor(lam, R, 1.0, b) for b in (-1.0, -0.5, 0.5, 1.0)]
    assert got == pytest.approx(list(_PROBE_F[key]), abs=6e-4)


def test_factor_is_kirsch_closed_form_for_quasi_isotropic():
    """F = gx / (gx + beta gy), gx = 1 + R^2/2r^2 + 3R^4/2r^4 and the
    sigma_y-induced gy = (R^2/r^2 - 3R^4/r^4) / 2 at r = R + d0."""
    lam = _lam(QI8)
    for R in (0.5, 1.4, 3.0):
        r = R + 1.0
        gx = 1 + R**2 / (2 * r**2) + 3 * R**4 / (2 * r**4)
        gy = (R**2 / r**2 - 3 * R**4 / r**4) / 2
        for beta in (-1.0, 0.4, 1.0):
            assert biaxial_point_stress_factor(lam, R, 1.0, beta) == pytest.approx(
                gx / (gx + beta * gy), rel=1e-9
            )


@pytest.mark.parametrize("name", sorted(MATERIAL_LIBRARY))
@pytest.mark.parametrize("layup", sorted(LAYUPS))
def test_factor_stays_near_one_at_bvid_hole_sizes(name, layup):
    """At the TAI hole sizes (R = 0.5-2 mm) the factor is within about 10%
    of 1 for every preset and layup at |beta| <= 1 and never undefined."""
    lam = Laminate(MATERIAL_LIBRARY[name], list(LAYUPS[layup]), 0.15)
    d0 = MATERIAL_LIBRARY[name].wn_d0_mm
    for R in (0.5, 1.0, 2.0):
        for beta in (-1.0, 1.0):
            assert 0.9 < biaxial_point_stress_factor(lam, R, d0, beta) < 1.1


def test_factor_is_undefined_when_the_point_stress_is_not_tensile():
    """A strongly orthotropic card under biaxial tension: sqrt(Ex/Ey) * beta
    exceeds the uniaxial concentration at the point."""
    a = (1.0 / 200000.0, 1.0 / 2000.0, -0.3 / 200000.0, 1.0 / 150000.0)
    gx = float(hole_stress_field(a, (1.0, 0.0, 0.0), 5.0, 0.0, 5.0)[0])
    gy = float(hole_stress_field(a, (0.0, 1.0, 0.0), 5.0, 0.0, 5.0)[0])
    assert gx + gy < 0  # the setup of the case below, at the hole edge

    class _FakeLam:
        def effective_engineering_constants(self):
            return (1 / a[0], 1 / a[1], 1 / a[3], -a[2] / a[0])

    with pytest.raises(ValueError, match="not tensile"):
        biaxial_point_stress_factor(_FakeLam(), 5.0, 1e-9, 1.0)  # type: ignore[arg-type]


def test_bad_inputs_are_rejected():
    a = laminate_compliance(_lam(QI8))
    with pytest.raises(ValueError, match="radius"):
        hole_stress_field(a, (1.0, 0.0, 0.0), 0.0, 1.0, 1.0)
    with pytest.raises(ValueError, match="positive definite"):
        hole_stress_field((1.0, 1.0, 2.0, 1.0), (1.0, 0.0, 0.0), 1.0, 0.0, 2.0)
    with pytest.raises(ValueError, match="far_field"):
        hole_stress_field(a, (1.0, 0.0), 1.0, 0.0, 2.0)  # type: ignore[arg-type]
