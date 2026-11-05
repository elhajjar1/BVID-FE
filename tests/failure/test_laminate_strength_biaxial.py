"""Pristine strength along a biaxial load ray (issue #115, PR 1).

``unnotched_strength(lam, loading, transverse_stress_ratio)`` and
``unnotched_strength_along(lam, direction)`` run the ply-discount model with
the load ``h * direction`` on the right-hand side of the CLT solve. At
``transverse_stress_ratio == 0`` it must be the uniaxial model bit for bit;
the reference below is the uniaxial code as it was before the change.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from bvidfe.core.laminate import Laminate, _reduced_stiffness, _transform_reduced_stiffness
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.laminate_strength import (
    PLY_DISCOUNT_FACTOR,
    _load_to,
    unnotched_strength,
    unnotched_strength_along,
)

_IM7 = MATERIAL_LIBRARY["IM7/8552"]
QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
CROSS = [0, 90, 0, 90, 90, 0, 90, 0]
HARD = [0, 45, 0, -45, 90, 90, -45, 0, 45, 0]
LAYUPS = {
    "QI8": QI8,
    "cross-ply": CROSS,
    "hard": HARD,
    "[0]4": [0, 0, 0, 0],
    "[90]4": [90, 90, 90, 90],
    "[+-45]s": [45, -45, -45, 45],
    "QI 30/60": [0, 30, 60, 90, -60, -30, -30, -60, 90, 60, 30, 0],
    "unbalanced": [0, 30, 90, 30, 0],
}


def _reference_uniaxial(lam: Laminate, loading: str) -> float:
    """``unnotched_strength`` as it was before biaxial loading (e80b811)."""
    m = lam.material
    sign = 1.0 if loading == "tension" else -1.0

    def ply_stresses(cracked):
        Q_intact = _reduced_stiffness(m)
        Q_cracked = Q_intact.copy()
        Q_cracked[1, :2] *= PLY_DISCOUNT_FACTOR
        Q_cracked[0, 1] *= PLY_DISCOUNT_FACTOR
        Q_cracked[2, 2] *= PLY_DISCOUNT_FACTOR
        Qs = [Q_cracked if c else Q_intact for c in cracked]
        A = sum(
            _transform_reduced_stiffness(Q, math.radians(theta)) * t
            for Q, theta, t in zip(Qs, lam.layup_deg, lam.ply_thicknesses_mm)
        )
        eps_x, eps_y, gamma_xy = np.linalg.solve(A, np.array([sign * lam.thickness_mm, 0.0, 0.0]))
        out = np.empty((lam.n_plies, 3))
        for k, (Q, theta) in enumerate(zip(Qs, lam.layup_deg)):
            c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
            eps_1 = c * c * eps_x + s * s * eps_y + c * s * gamma_xy
            eps_2 = s * s * eps_x + c * c * eps_y - c * s * gamma_xy
            gamma_12 = 2 * c * s * (eps_y - eps_x) + (c * c - s * s) * gamma_xy
            out[k] = Q @ np.array([eps_1, eps_2, gamma_12])
        return out

    cracked = [False] * lam.n_plies
    load = 0.0
    while True:
        s = ply_stresses(cracked)
        fibre = _load_to(m.Xt, m.Xc, s[:, 0])
        matrix = np.minimum(_load_to(m.Yt, m.Yc, s[:, 1]), _load_to(m.S12, m.S12, s[:, 2]))
        matrix[cracked] = np.inf
        if fibre.min() <= matrix.min():
            return max(float(fibre.min()), load)
        load = max(float(matrix.min()), load)
        for k in np.flatnonzero(matrix <= load * (1.0 + 1e-9)):
            cracked[k] = True
        if all(cracked):
            return load


def _lam(layup, material=_IM7, t=0.152):
    return Laminate(material, list(layup), t)


@pytest.mark.parametrize("material", sorted(MATERIAL_LIBRARY))
@pytest.mark.parametrize("layup", sorted(LAYUPS))
@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_zero_ratio_is_the_uniaxial_model_bit_for_bit(material, layup, loading):
    lam = _lam(LAYUPS[layup], MATERIAL_LIBRARY[material])
    ref = _reference_uniaxial(lam, loading)
    assert unnotched_strength(lam, loading) == ref
    assert unnotched_strength(lam, loading, 0.0) == ref
    assert unnotched_strength(lam, loading, -0.0) == ref
    assert unnotched_strength(lam, loading, 0) == ref
    sign = 1.0 if loading == "tension" else -1.0
    assert unnotched_strength_along(lam, (sign, 0.0, 0.0)) == ref


def test_per_ply_thicknesses_are_bit_identical_too():
    lam = Laminate(_IM7, QI8, [0.1, 0.2, 0.15, 0.1, 0.1, 0.15, 0.2, 0.1])
    for loading in ("compression", "tension"):
        assert unnotched_strength(lam, loading, 0.0) == _reference_uniaxial(lam, loading)


# Probe 2.1 of the #115 design: QI8 IM7/8552, |sigma_x| at last-ply failure.
_QI8_COMPRESSION = {
    -1.0: 411.6,
    -0.5: 511.8,
    -0.25: 560.5,
    0.0: 601.3,
    0.25: 648.5,
    0.5: 703.8,
    1.0: 848.3,
}
_QI8_TENSION = {
    -1.0: 411.6,
    -0.5: 645.3,
    -0.25: 692.9,
    0.0: 753.2,
    0.25: 886.0,
    0.5: 889.6,
    1.0: 580.8,  # every ply cracked first: the D5 artefact, see below
}


@pytest.mark.parametrize("beta", sorted(_QI8_COMPRESSION))
def test_qi8_envelope_is_pinned(beta):
    lam = _lam(QI8)
    assert unnotched_strength(lam, "compression", beta) == pytest.approx(
        _QI8_COMPRESSION[beta], abs=0.05
    )
    assert unnotched_strength(lam, "tension", beta) == pytest.approx(_QI8_TENSION[beta], abs=0.05)


def test_cross_ply_compression_is_governed_by_the_0_deg_fibres():
    """[0/90]2s: 822.2-848.3 MPa over the ratio range (probe 2.1)."""
    lam = _lam(CROSS)
    values = [unnotched_strength(lam, "compression", b) for b in (-1, -0.5, 0, 0.5, 1)]
    assert values[0] == pytest.approx(822.2, abs=0.05)
    assert values[2] == pytest.approx(835.0, abs=0.05)
    assert values[-1] == pytest.approx(848.3, abs=0.05)


def test_equibiaxial_tension_ends_when_every_ply_has_cracked():
    """The D5 artefact, kept on purpose and documented under Limitations:
    under equibiaxial tension every QI ply cracks at the same load, and the
    laminate is taken to fail there (580.8 MPa), as a [90]n laminate does
    uniaxially, although the fibre network could carry about Xt/2 = 1,119
    MPa. So the tension envelope drops from beta = 0.5 to 1."""
    lam = _lam(QI8)
    equi = unnotched_strength(lam, "tension", 1.0)
    assert equi == pytest.approx(580.8, abs=0.05)
    assert equi < unnotched_strength(lam, "tension", 0.5)
    assert equi < 0.6 * _IM7.Xt  # far below the netting estimate Xt / 2
    # Uniaxial [90]n is the same rule, and right there (no fibre path).
    assert unnotched_strength(_lam([90] * 4), "tension") == pytest.approx(_IM7.Yt, rel=1e-9)


def test_pure_shear_at_45_deg_is_the_same_for_both_senses():
    """beta = -1 is (s, -s): pure shear on the +-45 deg planes, the same
    load for either sign of s on a symmetric balanced layup."""
    lam = _lam(QI8)
    c = unnotched_strength(lam, "compression", -1.0)
    t = unnotched_strength(lam, "tension", -1.0)
    assert c == pytest.approx(t, rel=1e-12)


@pytest.mark.parametrize("layup", ["QI8", "cross-ply", "hard", "[+-45]s", "QI 30/60"])
@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_equibiaxial_strength_does_not_change_when_the_layup_turns_90_deg(layup, loading):
    lam = _lam(LAYUPS[layup])
    turned = _lam([a + 90 for a in LAYUPS[layup]])
    assert unnotched_strength(lam, loading, 1.0) == pytest.approx(
        unnotched_strength(turned, loading, 1.0), rel=1e-9
    )


@pytest.mark.parametrize("layup", ["QI8", "cross-ply", "hard"])
@pytest.mark.parametrize("beta", [-0.7, 0.3, 0.9])
def test_a_load_along_y_is_the_load_along_x_of_the_layup_turned_90_deg(layup, beta):
    """The rotation workaround of the config's error message: the strength
    along (s, s*beta) of the layup equals the strength along (s*beta, s) of
    the layup turned by 90 deg."""
    lam = _lam(LAYUPS[layup])
    turned = _lam([a + 90 for a in LAYUPS[layup]])
    for s in (-1.0, 1.0):
        along_x = unnotched_strength_along(lam, (s, s * beta, 0.0))
        along_y = unnotched_strength_along(turned, (s * beta, s, 0.0))
        assert along_x == pytest.approx(along_y, rel=1e-9)


def test_the_multiplier_scales_inversely_with_the_direction():
    lam = _lam(HARD)
    d = (-1.0, -0.4, 0.0)
    assert unnotched_strength_along(lam, (2 * d[0], 2 * d[1], 0.0)) == pytest.approx(
        unnotched_strength_along(lam, d) / 2, rel=1e-12
    )
    assert unnotched_strength(lam, "compression", 0.4) == unnotched_strength_along(lam, d)


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("layup", ["QI8", "cross-ply", "hard"])
def test_a_tiny_ratio_is_continuous_with_the_uniaxial_strength(loading, layup):
    lam = _lam(LAYUPS[layup])
    base = unnotched_strength(lam, loading)
    for beta in (1e-9, -1e-9):
        assert unnotched_strength(lam, loading, beta) == pytest.approx(base, rel=1e-6)


def test_in_plane_shear_is_a_direction_too():
    """tau_xy is accepted; pure shear on QI equals the beta = -1 load turned
    by 45 deg, so the same multiplier."""
    lam = _lam(QI8)
    shear = unnotched_strength_along(lam, (0.0, 0.0, 1.0))
    assert shear == pytest.approx(unnotched_strength(lam, "tension", -1.0), rel=1e-6)


@pytest.mark.parametrize(
    "direction", [(0.0, 0.0, 0.0), (1.0, math.nan, 0.0), (math.inf, 0.0, 0.0), (1.0, 0.0)]
)
def test_unnotched_strength_along_rejects_a_bad_direction(direction):
    with pytest.raises(ValueError, match="direction"):
        unnotched_strength_along(_lam(QI8), direction)
