"""The closed-form buckling kernels of #132 PR0: the simply supported panel
Rayleigh-Ritz with the bending-twisting terms D16 and D26, and the clamped
ellipse Rayleigh-Ritz used for the delaminated sublaminate."""

import math

import numpy as np
import pytest

from bvidfe.analysis.semi_analytical import (
    _clamped_ellipse_load_factor,
    _ssss_ritz_N_cr,
    panel_buckling_load,
)
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY

_IM7 = MATERIAL_LIBRARY["IM7/8552"]
_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
_CP8 = [0, 90, 90, 0, 0, 90, 90, 0]


def _D(layup):
    return Laminate(_IM7, layup, 0.152).abd_matrices()[2]


def _navier(D, a, b):
    """Navier SSSS uniaxial N_cr (D16 = D26 = 0), minimised over m, n."""
    r = a / b
    return min(
        (math.pi / a) ** 2
        * (
            D[0, 0] * m**4
            + 2.0 * (D[0, 1] + 2.0 * D[2, 2]) * (m * r) ** 2 * n**2
            + D[1, 1] * (r * n) ** 4
        )
        / m**2
        for m in range(1, 40)
        for n in range(1, 4)
    )


def _orthotropic(D):
    D0 = np.array(D, dtype=float)
    D0[0, 2] = D0[2, 0] = D0[1, 2] = D0[2, 1] = 0.0
    return D0


def _isotropic_D(nu=0.3):
    d = 1.0 / (12.0 * (1.0 - nu**2))
    return np.array([[d, nu * d, 0.0], [nu * d, d, 0.0], [0.0, 0.0, 0.5 * (1.0 - nu) * d]]), d


# --- simply supported panel with D16 / D26 ---


@pytest.mark.parametrize("layup", [_CP8, [0] * 8, [90] * 8, [0, 90, 0, 90, 90, 0, 90, 0]])
@pytest.mark.parametrize("a, b", [(200.0, 150.0), (150.0, 100.0), (100.0, 150.0), (300.0, 50.0)])
def test_panel_ritz_equals_navier_without_bending_twisting(layup, a, b):
    """With D16 = D26 = 0 (cross-ply and unidirectional laminates) the
    Ritz value is the exact Navier closed form, multi-half-wave modes
    included (300 x 50 mm: m = 6 on the cross-ply)."""
    D = _D(layup)
    assert abs(D[0, 2]) < 1e-9 * D[0, 0] and abs(D[1, 2]) < 1e-9 * D[0, 0]
    assert _ssss_ritz_N_cr(D, a, b) == pytest.approx(_navier(D, a, b), rel=1e-7)


def test_panel_canonical_quasi_isotropic_value():
    """The canonical 200 x 150 mm 8-ply QI IM7/8552 SSSS panel (#129)
    converges to 12.204 MPa with D16 / D26 (#132 design probe: a
    geometric-BC polynomial Ritz, 12.2058 -> 12.2038 MPa for M = 6..14),
    3.0% below the Navier 12.582 MPa that ignores them."""
    lam = Laminate(_IM7, _QI8, 0.152)
    h = lam.thickness_mm
    sigma = panel_buckling_load(lam, 200.0, 150.0) / h
    assert sigma == pytest.approx(12.204, abs=1e-3)
    navier = _navier(lam.abd_matrices()[2], 200.0, 150.0) / h
    assert navier == pytest.approx(12.582, abs=1e-3)
    assert sigma / navier - 1.0 == pytest.approx(-0.030, abs=0.001)


@pytest.mark.parametrize(
    "layup", [_QI8, [45, -45, -45, 45], [45, -45, 45, -45, -45, 45, -45, 45], [30, -60, -60, 30]]
)
@pytest.mark.parametrize("a, b", [(200.0, 150.0), (150.0, 100.0), (100.0, 150.0), (300.0, 50.0)])
def test_bending_twisting_only_lowers_the_panel_load(layup, a, b):
    """The Navier mode is in the Ritz space and its D16 / D26 energy is
    zero, so D16 / D26 can only lower N_cr. The mirror-image layup
    (D16, D26 -> -D16, -D26) buckles at the same load."""
    D = _D(layup)
    with_twist = _ssss_ritz_N_cr(D, a, b)
    without = _ssss_ritz_N_cr(_orthotropic(D), a, b)
    assert without == pytest.approx(_navier(D, a, b), rel=1e-7)
    assert with_twist < without
    mirrored = _ssss_ritz_N_cr(_D([-p for p in layup]), a, b)
    assert mirrored == pytest.approx(with_twist, rel=1e-9)


def test_bending_twisting_effect_grows_with_the_twist_coupling():
    """On 200 x 150 mm the QI laminate (D16/D11 = 0.05) loses 3.0%; the
    [45/-45]s angle-ply (D16/D11 = 0.6) loses more."""
    qi, angle = _D(_QI8), _D([45, -45, -45, 45])
    loss_qi = 1.0 - _ssss_ritz_N_cr(qi, 200.0, 150.0) / _navier(qi, 200.0, 150.0)
    loss_angle = 1.0 - _ssss_ritz_N_cr(angle, 200.0, 150.0) / _navier(angle, 200.0, 150.0)
    assert loss_qi == pytest.approx(0.030, abs=0.001)
    assert loss_angle > 3.0 * loss_qi


def test_panel_boundary_factor_is_kept():
    lam = Laminate(_IM7, _QI8, 0.152)
    ss = panel_buckling_load(lam, 150.0, 100.0)
    assert panel_buckling_load(lam, 150.0, 100.0, "clamped") == pytest.approx(1.9 * ss, rel=1e-12)
    assert panel_buckling_load(lam, 150.0, 100.0, "free") == pytest.approx(0.5 * ss, rel=1e-12)


# --- clamped ellipse ---


def test_clamped_circle_under_radial_compression_is_exact():
    """Isotropic clamped circle under uniform radial compression:
    N R^2 / D = j_{1,1}^2 = 14.682, j_{1,1} the first zero of J1
    (Timoshenko & Gere 1961)."""
    D, d = _isotropic_D()
    j11 = 3.831705970207512
    for R in (1.0, 12.0):
        lam = _clamped_ellipse_load_factor(D, R, R, 1.0, 1.0, 0.0)
        assert lam * R**2 / d == pytest.approx(j11**2, rel=1e-6)


def test_clamped_circle_under_uniaxial_compression():
    """Uniaxial: between the radial value (it drops the y work) and twice it
    (the axisymmetric radial mode has equal x and y work), at 27.22 D / R^2
    (#132 design probe: 27.2)."""
    D, d = _isotropic_D()
    radial = _clamped_ellipse_load_factor(D, 1.0, 1.0, 1.0, 1.0, 0.0)
    uniaxial = _clamped_ellipse_load_factor(D, 1.0, 1.0, 1.0, 0.0, 0.0)
    assert radial < uniaxial < 2.0 * radial
    assert uniaxial / d == pytest.approx(27.22, abs=0.01)
    # In any direction, for an isotropic circle.
    for t in (0.3, 1.1, 2.0):
        c, s = math.cos(t), math.sin(t)
        assert _clamped_ellipse_load_factor(D, 1.0, 1.0, c * c, s * s, -c * s) == pytest.approx(
            uniaxial, rel=1e-9
        )


def test_clamped_ellipse_scales_with_size():
    D = _D([0, 45, -45])
    base = _clamped_ellipse_load_factor(D, 12.0, 7.0, 1.0, 0.0, 0.0)
    assert _clamped_ellipse_load_factor(D, 24.0, 14.0, 1.0, 0.0, 0.0) == pytest.approx(
        base / 4.0, rel=1e-10
    )


@pytest.mark.parametrize("material", ["isotropic", "IM7/8552 [0]"])
def test_elongated_across_the_load_tends_to_the_clamped_column(material):
    """Loaded across its long axis, an elongated clamped ellipse of
    half-width a along the load tends to the clamped-clamped wide column of
    length 2a, N = pi^2 D11 / a^2 (exact for cylindrical bending), from
    above (domain monotonicity: the ellipse lies inside the strip)."""
    D = _isotropic_D()[0] if material == "isotropic" else _D([0])
    column = math.pi**2 * D[0, 0]
    values = [_clamped_ellipse_load_factor(D, 1.0, k, 1.0, 0.0, 0.0) for k in (1, 2, 3.25, 5, 10)]
    assert all(v > column for v in values)
    assert all(v1 > v2 for v1, v2 in zip(values, values[1:]))
    assert values[-1] / column == pytest.approx(1.0, abs=0.06)


def test_elongated_along_the_load_tends_to_the_clamped_strip():
    """Loaded along its long axis, an isotropic clamped ellipse of
    half-width b tends to the long plate with clamped edges,
    N = 6.97 pi^2 D / (2b)^2 (Timoshenko & Gere 1961, k = 6.97), from above."""
    D, d = _isotropic_D()
    strip = 6.97 * math.pi**2 * d / 4.0
    values = [_clamped_ellipse_load_factor(D, k, 1.0, 1.0, 0.0, 0.0) for k in (1, 2, 3.25, 5, 10)]
    assert all(v > strip for v in values)
    assert all(v1 > v2 for v1, v2 in zip(values, values[1:]))
    assert values[-1] / strip == pytest.approx(1.0, abs=0.06)


def test_clamped_ellipse_is_frame_independent():
    """The same plate described with its axes swapped (x' <-> y'): the
    semi-axes, the stiffness indices and the load components swap."""
    D = _D([0, 45, -45])
    swap = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    D_swapped = swap @ D @ swap
    for n in ((1.0, 0.0, 0.0), (0.75, 0.25, -0.433)):
        assert _clamped_ellipse_load_factor(D, 15.0, 6.0, *n) == pytest.approx(
            _clamped_ellipse_load_factor(D_swapped, 6.0, 15.0, n[1], n[0], n[2]), rel=1e-9
        )
