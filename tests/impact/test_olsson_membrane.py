"""Large-deflection plate law and exact energies in the Olsson onset energy.

The plate follows ``P = k_b w + k_m w^3`` (single-term von Karman solution,
``k_m`` from ``_k_membrane_single_term``) and the contact ``P = K_c
delta^n``; the onset energy is the energy both store at the threshold load.
These tests pin the closed forms (the isotropic simply supported plate,
an independent Ritz minimisation of the membrane energy for orthotropic
plates and the clamped mode), the cubic solve, the small-deflection limit
and the exact contact energies.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact import olsson
from bvidfe.impact.olsson import (
    _contact_energy,
    _contact_indentation,
    _k_bending,
    _k_membrane,
    _k_membrane_single_term,
    onset_energy,
    plate_deflection,
    threshold_load,
)

_M = MATERIAL_LIBRARY["IM7/8552"]
_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
_QI16 = [0, 45, -45, 90] * 4


def _iso_A(E: float, nu: float, h: float) -> np.ndarray:
    c = E * h / (1.0 - nu**2)
    return c * np.array([[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, (1.0 - nu) / 2.0]])


def _ritz_membrane(A, a, b, mode, immovable, n_terms=5, n_quad=40):
    """Independent check: minimise the membrane energy of w = W phi over a
    Fourier in-plane field, ``U_m = C W^4``, and return ``K = 4 C`` (so
    ``P = k1 W + K W^3`` for a central load, ``phi(centre) = 1``)."""
    xg, wg = np.polynomial.legendre.leggauss(n_quad)
    xs, wx = (xg + 1) * a / 2, wg * a / 2
    ys, wy = (xg + 1) * b / 2, wg * b / 2
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    wq = np.outer(wx, wy)
    al, be = math.pi / a, math.pi / b
    if mode == "sine":
        wx_ = al * np.cos(al * X) * np.sin(be * Y)
        wy_ = be * np.sin(al * X) * np.cos(be * Y)
    else:
        wx_ = al * np.sin(2 * al * X) * np.sin(be * Y) ** 2
        wy_ = be * np.sin(al * X) ** 2 * np.sin(2 * be * Y)
    e0 = np.stack([0.5 * wx_**2, 0.5 * wy_**2, wx_ * wy_])
    z = np.zeros_like(X)
    cols = []
    for m in range(1, n_terms + 1):
        for n in range(0, n_terms + 1):
            # u = sin(m al x) cos(n be y): no normal displacement at x = 0, a
            ux = m * al * np.cos(m * al * X) * np.cos(n * be * Y)
            uy = -n * be * np.sin(m * al * X) * np.sin(n * be * Y)
            cols.append(np.stack([ux, z, uy]))
            # v = cos(n al x) sin(m be y): none at y = 0, b
            vy = m * be * np.cos(n * al * X) * np.cos(m * be * Y)
            vx = -n * al * np.sin(n * al * X) * np.sin(m * be * Y)
            cols.append(np.stack([z, vy, vx]))
    if not immovable:  # edges free to move, kept straight: u = e_x x, v = e_y y
        cols.append(np.stack([np.ones_like(X), z, z]))
        cols.append(np.stack([z, np.ones_like(X), z]))
    B = np.stack(cols, axis=-1)
    AB = np.einsum("ij,jxyk->ixyk", A, B)
    Ae0 = np.einsum("ij,jxy->ixy", A, e0)
    H = np.einsum("ixyk,ixyl,xy->kl", B, AB, wq)
    g = np.einsum("ixyk,ixy,xy->k", B, Ae0, wq)
    c = np.linalg.lstsq(H, -g, rcond=None)[0]
    u_min = 0.5 * np.einsum("ixy,ixy,xy->", e0, Ae0, wq) + 0.5 * g @ c
    return 4.0 * u_min


# ---------- membrane coefficient ----------


@pytest.mark.parametrize("a, b", [(100.0, 100.0), (125.0, 75.0)])
def test_isotropic_square_sine_mode_closed_form(a, b):
    """Levy (1942) single-term solution of a simply supported isotropic
    plate: ``K = (a b / 64) E h (alpha^4 + beta^4)`` with edges movable
    (kept straight); holding them adds
    ``(a b / 32) E h / (1 - nu^2) (alpha^4 + 2 nu alpha^2 beta^2 + beta^4)``.
    For a square plate the totals are ``pi^4 E h / (32 a^2)`` and
    ``pi^4 E h (3 - nu) / (32 a^2 (1 - nu))``."""
    E, nu, h = 70e3, 0.3, 2.0
    A = _iso_A(E, nu, h)
    al, be = math.pi / a, math.pi / b
    k_mov = a * b / 64 * E * h * (al**4 + be**4)
    k_imm = k_mov + a * b / 32 * E * h / (1 - nu**2) * (al**4 + 2 * nu * al**2 * be**2 + be**4)
    assert _k_membrane_single_term(A, a, b, "sine", False) == pytest.approx(k_mov, rel=1e-12)
    assert _k_membrane_single_term(A, a, b, "sine", True) == pytest.approx(k_imm, rel=1e-12)
    if a == b:
        assert k_mov == pytest.approx(math.pi**4 * E * h / (32 * a**2), rel=1e-12)
        assert k_imm == pytest.approx(
            math.pi**4 * E * h * (3 - nu) / (32 * a**2 * (1 - nu)), rel=1e-12
        )


@pytest.mark.parametrize("mode", ["sine", "clamped"])
@pytest.mark.parametrize("immovable", [False, True])
@pytest.mark.parametrize("layup", [_QI16, [0, 90] * 4, [0, 0, 45, -45, 0, 0]])
def test_membrane_closed_form_matches_independent_ritz(mode, immovable, layup):
    """The stress-function closed form equals a direct minimisation of the
    membrane energy over the in-plane displacements, for orthotropic
    laminates and both assumed modes."""
    lam = Laminate(_M, layup, 0.152)
    A, _, _ = lam.abd_matrices()
    A = A.copy()
    A[0, 2] = A[2, 0] = A[1, 2] = A[2, 1] = 0.0  # the closed form ignores A16, A26
    for a, b in ((150.0, 100.0), (60.0, 60.0)):
        k_cf = _k_membrane_single_term(A, a, b, mode, immovable)
        k_ritz = _ritz_membrane(A, a, b, mode, immovable)
        assert k_cf == pytest.approx(k_ritz, rel=1e-8)


def test_k_membrane_per_boundary():
    """simply_supported: sine mode, movable edges; clamped: clamped mode,
    immovable edges; free: no membrane stiffness."""
    lam = Laminate(_M, _QI16, 0.152)
    A, _, _ = lam.abd_matrices()
    ss = _k_membrane(lam, PanelGeometry(150, 100, "simply_supported"))
    cl = _k_membrane(lam, PanelGeometry(150, 100, "clamped"))
    assert ss == pytest.approx(_k_membrane_single_term(A, 150, 100, "sine", False))
    assert cl == pytest.approx(_k_membrane_single_term(A, 150, 100, "clamped", True))
    assert _k_membrane(lam, PanelGeometry(150, 100, "free")) == 0.0
    assert 0.0 < ss < cl


def test_k_membrane_scales_with_thickness_and_inverse_square_of_size():
    """``k_m`` is proportional to ``A`` (so to ``h``) and to ``1 / L^2``."""
    lam1 = Laminate(_M, _QI16, 0.152)
    lam2 = Laminate(_M, _QI16, 0.304)
    pan = PanelGeometry(150, 100)
    assert _k_membrane(lam2, pan) / _k_membrane(lam1, pan) == pytest.approx(2.0, rel=1e-9)
    big = PanelGeometry(300, 200)
    assert _k_membrane(lam1, pan) / _k_membrane(lam1, big) == pytest.approx(4.0, rel=1e-9)


# ---------- the cubic ----------


@pytest.mark.parametrize("k_m", [1e-12, 1e-3, 1.0, 1e3, 1e6])
def test_plate_deflection_solves_the_cubic(k_m):
    k_b, P = 500.0, 3000.0
    w = plate_deflection(k_b, k_m, P)
    assert k_b * w + k_m * w**3 == pytest.approx(P, rel=1e-12)
    assert 0.0 < w <= P / k_b


def test_plate_deflection_limits():
    k_b, P = 500.0, 3000.0
    assert plate_deflection(k_b, 0.0, P) == P / k_b
    assert plate_deflection(k_b, 1e-15, P) == pytest.approx(P / k_b, rel=1e-9)
    # membrane-dominated: w -> (P / k_m)^(1/3)
    assert plate_deflection(1e-6, 10.0, P) == pytest.approx((P / 10.0) ** (1 / 3), rel=1e-6)
    assert plate_deflection(k_b, 1.0, 0.0) == 0.0
    with pytest.raises(ValueError):
        plate_deflection(0.0, 1.0, P)
    with pytest.raises(ValueError):
        plate_deflection(k_b, -1.0, P)


# ---------- contact energy ----------


@pytest.mark.parametrize(
    "shape, fraction", [("hemispherical", 2.0 / 5.0), ("flat", 1.0 / 2.0), ("conical", 1.0 / 3.0)]
)
def test_contact_energy_is_the_integral_of_the_contact_law(shape, fraction):
    """``int_0^delta P d delta`` for ``P = K delta^n`` is ``P delta / (n + 1)``:
    2/5 for Hertz (not the secant's 1/2), 1/2 for a flat punch and 1/3 for
    a cone (not the tangent's 1/4)."""
    lam = Laminate(_M, _QI16, 0.152)
    imp = ImpactorGeometry(16.0, shape=shape)
    P = 2500.0
    delta = _contact_indentation(lam, imp, P)
    assert _contact_energy(lam, imp, P) == pytest.approx(fraction * P * delta, rel=1e-12)
    ds = np.linspace(0.0, delta, 20001)
    Ps = np.array([_P_of_delta(lam, imp, d) for d in ds])
    assert _contact_energy(lam, imp, P) == pytest.approx(np.trapezoid(Ps, ds), rel=1e-6)


def _P_of_delta(lam, imp, d):
    K_c, n = olsson._contact_law(lam, imp)
    return K_c * d**n


def test_hertz_contact_law():
    """``P = (4/3) E* sqrt(R) delta^1.5`` (Johnson 1985, §4.2)."""
    lam = Laminate(_M, _QI16, 0.152)
    imp = ImpactorGeometry(16.0)
    E_eff = olsson._contact_E_eff(lam)
    delta = _contact_indentation(lam, imp, 2000.0)
    assert (4 / 3) * E_eff * math.sqrt(8.0) * delta**1.5 == pytest.approx(2000.0, rel=1e-12)


# ---------- onset energy ----------


def _onset_parts(lam, pan, imp):
    Pc = threshold_load(lam, pan, imp)
    k_b = _k_bending(lam, pan, pan.Lx_mm / 2, pan.Ly_mm / 2)
    return Pc, k_b, _contact_energy(lam, imp, Pc)


def test_onset_energy_is_the_stored_energy_at_the_threshold_load():
    lam = Laminate(_M, _QI8, 0.152)
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    Pc, k_b, e_c = _onset_parts(lam, pan, imp)
    k_m = _k_membrane(lam, pan)
    w = plate_deflection(k_b, k_m, Pc)
    expected = (0.5 * k_b * w**2 + 0.25 * k_m * w**4 + e_c) * 1e-3
    assert onset_energy(lam, pan, imp) == pytest.approx(expected, rel=1e-12)
    # and it equals the work of the plate law integrated numerically
    ws = np.linspace(0.0, w, 20001)
    work = np.trapezoid(k_b * ws + k_m * ws**3, ws)
    assert (work + e_c) * 1e-3 == pytest.approx(expected, rel=1e-6)


def test_onset_energy_small_deflection_limit_is_the_linear_result(monkeypatch):
    """With ``k_m = 0`` the plate part is ``Pc^2 / (2 k_b)``."""
    lam = Laminate(_M, _QI8, 0.152)
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    Pc, k_b, e_c = _onset_parts(lam, pan, imp)
    monkeypatch.setattr(olsson, "_k_membrane", lambda lam, pan: 0.0)
    assert onset_energy(lam, pan, imp) == pytest.approx((Pc**2 / (2 * k_b) + e_c) * 1e-3)


def test_membrane_term_negligible_when_deflection_is_small():
    """A thick laminate on a small window deflects a fraction of its
    thickness at Pc: the membrane term changes the plate energy < 3%."""
    lam = Laminate(_M, [0, 45, -45, 90] * 6 + [90, -45, 45, 0] * 6, 0.152)  # 48 plies, 7.3 mm
    pan = PanelGeometry(100, 100)
    imp = ImpactorGeometry()
    Pc, k_b, e_c = _onset_parts(lam, pan, imp)
    w_lin = Pc / k_b
    assert w_lin / lam.thickness_mm < 0.3
    e_plate = onset_energy(lam, pan, imp) * 1e3 - e_c
    assert e_plate == pytest.approx(Pc**2 / (2 * k_b), rel=0.03)
    assert e_plate < Pc**2 / (2 * k_b)


def test_membrane_stiffening_lowers_the_onset_of_thin_large_panels():
    """An 8-ply plate deflects many thicknesses at Pc on the linear law
    (14x on 150 x 100 mm, 56x on 300 x 200 mm); membrane stretching
    carries the load at a few thicknesses, so the onset energy is far below
    ``Pc^2 / (2 k_b)`` and grows only slowly with the panel size."""
    lam = Laminate(_M, _QI8, 0.152)
    imp = ImpactorGeometry()
    onsets = []
    for a, b in ((150, 100), (300, 200)):
        pan = PanelGeometry(a, b)
        Pc, k_b, e_c = _onset_parts(lam, pan, imp)
        w = plate_deflection(k_b, _k_membrane(lam, pan), Pc)
        assert 1.5 < w / lam.thickness_mm < 6.0
        e = onset_energy(lam, pan, imp)
        assert e < 0.25 * (Pc**2 / (2 * k_b) + e_c) * 1e-3
        onsets.append(e)
    assert onsets[0] == pytest.approx(0.84, abs=0.01)
    assert onsets[1] == pytest.approx(1.24, abs=0.01)


def test_free_panel_has_no_membrane_term():
    lam = Laminate(_M, _QI16, 0.152)
    pan = PanelGeometry(150, 100, "free")
    imp = ImpactorGeometry()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Pc, k_b, e_c = _onset_parts(lam, pan, imp)
        e = onset_energy(lam, pan, imp)
    assert e == pytest.approx((Pc**2 / (2 * k_b) + e_c) * 1e-3, rel=1e-12)
