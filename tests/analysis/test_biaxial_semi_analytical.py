"""Biaxial loading in the semi_analytical tier (issue #115, PR 2).

The film of the delaminated sublaminate takes the laminate's biaxial
strain: its compression along x is the x-row of ``A_sub a* (s, s beta, 0) /
h`` and along y the y-row of ``A_sub a* (0, s beta, 0) / h`` (Option A of
the design: the Poisson-mismatch y-stress and the film shear are left out,
as at the uniaxial load). That biaxial prestress enters the clamped-ellipse
onset Ritz, the post-buckled 2D blister and the whole-panel Ritz. Checked
here:

- ``beta == 0`` runs the unchanged code, bit for bit, in every function;
- closed forms: the isotropic clamped circle under equibiaxial compression
  (``j_{1,1}^2 D / R^2``) and its exact axisymmetric blister (Evans &
  Hutchinson 1984), the isotropic simply supported square under
  equibiaxial compression (``2 pi^2 D / a^2``, half the uniaxial value),
  the biaxial Navier formula for a cross-ply panel, and an isotropic film
  whose stress is the laminate's exactly;
- lateral tension stabilises the film, lateral compression destabilises
  it, and a film in tension along both axes does not buckle;
- a 90 deg rotation of the load and the ellipse gives the same growth;
- continuity at ``beta -> 0``;
- the pruned search finds the full minimum;
- pinned knockdowns, and tension with lateral compression (the growth
  channel of D10 in the design);
- ``semi_analytical <= empirical`` (also in
  ``tests/test_cross_tier_consistency.py``).

Nothing here is validated against measurements: there are no biaxial
residual-strength data.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import replace

import numpy as np
import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, PanelGeometry
from bvidfe.analysis import semi_analytical as sa
from bvidfe.analysis.semi_analytical import (
    _clamped_ellipse_load_factor,
    _film_compression_per_unit_load,
    _PostbuckledFilm,
    _ssss_ritz_load_factor,
    _sublaminate_film,
    mixed_mode_toughness,
    panel_buckling_load,
    semi_analytical_cai,
    semi_analytical_tai,
    sublaminate_buckling_load,
    sublaminate_growth_stress,
    weakest_sublaminate_growth,
)
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from tests.analysis.test_blister_energy_release import _exact_circular_blister

_IM7 = MATERIAL_LIBRARY["IM7/8552"]
QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
CP8 = [0, 90, 0, 90, 90, 0, 90, 0]
_J11 = 3.831705970207512  # first zero of J1

# An in-plane isotropic ply: the film stress is the laminate stress, so
# the film load of every load ray is known exactly.
_E, _NU = 70000.0, 0.3
_ISO = replace(
    _IM7,
    name="isotropic",
    E11=_E,
    E22=_E,
    nu12=_NU,
    G12=_E / (2 * (1 + _NU)),
    G13=_E / (2 * (1 + _NU)),
    G23=_E / (2 * (1 + _NU)),
)
_T = 0.25  # ply thickness (mm)
_D_PLY = _E * _T**3 / (12 * (1 - _NU**2))


def _iso_lam(n_plies: int = 4) -> Laminate:
    return Laminate(_ISO, [0] * n_plies, _T)


def _cfg(layup, energy, loading, beta, tier="semi_analytical"):
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=list(layup),
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        impact=ImpactEvent(energy_J=energy),
        tier=tier,
        loading=loading,
        transverse_stress_ratio=beta,
    )


def _run(cfg):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return BvidAnalysis(cfg).run()


_DAMAGE_CACHE: dict = {}


def _damage(layup, energy) -> DamageState:
    key = (tuple(layup), energy)
    if key not in _DAMAGE_CACHE:
        _DAMAGE_CACHE[key] = _run(_cfg(layup, energy, "compression", 0.0, "empirical")).damage
    return _DAMAGE_CACHE[key]


def _largest(damage: DamageState) -> dict:
    out: dict = {}
    for e in damage.delaminations:
        if e.interface_index not in out or e.area_mm2 > out[e.interface_index].area_mm2:
            out[e.interface_index] = e
    return out


# --- beta = 0 runs the unchanged code --------------------------------------


@pytest.fixture
def no_biaxial_code(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("biaxial code called at the uniaxial load")

    for name in (
        "_biaxial_growth_terms",
        "_film_compression_per_unit_load",
        "_biaxial_film_load",
        "_ssss_ritz_load_factor",
        "_weakest_biaxial_growth",
        "_semi_analytical_tai_result",
    ):
        monkeypatch.setattr(sa, name, boom)


@pytest.mark.parametrize("zero", [0.0, -0.0])
def test_zero_ratio_is_the_uniaxial_code_in_every_function(zero, no_biaxial_code):
    lam = Laminate(_IM7, QI8, 0.152)
    damage = _damage(QI8, 10.0)
    e = _largest(damage)[5]
    kw = dict(loading="compression", transverse_stress_ratio=zero)
    assert sublaminate_buckling_load(lam, e, **kw) == sublaminate_buckling_load(lam, e)
    assert sublaminate_growth_stress(lam, e, **kw) == sublaminate_growth_stress(lam, e)
    assert weakest_sublaminate_growth(lam, damage, **kw) == weakest_sublaminate_growth(lam, damage)
    assert panel_buckling_load(lam, 150.0, 100.0, **kw) == panel_buckling_load(lam, 150.0, 100.0)
    assert semi_analytical_cai(
        lam, damage, 601.3, transverse_stress_ratio=zero
    ) == semi_analytical_cai(lam, damage, 601.3)
    assert semi_analytical_tai(
        lam, damage, 753.2, 100.0, transverse_stress_ratio=zero
    ) == semi_analytical_tai(lam, damage, 753.2, 100.0)
    assert semi_analytical_cai(lam, damage, 601.3).growth is None


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_zero_ratio_run_does_not_touch_the_biaxial_code(loading, no_biaxial_code):
    r = _run(_cfg(QI8, 10.0, loading, -0.0))
    assert r.biaxial_interaction_factor is None
    assert not any("Biaxial" in n for n in r.notes)


def test_the_new_growth_fields_default_to_the_uniaxial_load():
    lam = Laminate(_IM7, QI8, 0.152)
    g = sublaminate_growth_stress(lam, _largest(_damage(QI8, 10.0))[5])
    assert (g.film_stress_axis, g.film_stress_ratio) == ("x", 0.0)


def test_arguments_are_validated():
    lam = Laminate(_IM7, QI8, 0.152)
    e = _largest(_damage(QI8, 10.0))[5]
    with pytest.raises(ValueError, match="loading"):
        sublaminate_buckling_load(lam, e, loading="biaxial")
    with pytest.raises(ValueError, match="transverse_stress_ratio"):
        sublaminate_growth_stress(lam, e, transverse_stress_ratio=float("nan"))
    with pytest.raises(ValueError, match="boundary"):
        semi_analytical_tai(
            lam, _damage(QI8, 10.0), 753.2, 100.0, boundary="pinned", transverse_stress_ratio=-0.5
        )


# --- the film load (Option A) ----------------------------------------------


@pytest.mark.parametrize("loading, s", [("compression", -1.0), ("tension", 1.0)])
@pytest.mark.parametrize("beta", [-1.0, -0.3, 0.6, 1.0])
def test_film_with_the_laminates_layup_takes_the_laminate_stress(loading, s, beta):
    """For a film with the laminate's own in-plane stiffness per unit
    thickness (an isotropic ply, or [0/90] in a [0/90]2s laminate) ``A_sub a*
    = h_sub I``: the film stress is the laminate stress, and Option A leaves
    nothing out."""
    for lam, iface in ((_iso_lam(), 0), (Laminate(_IM7, CP8, 0.152), 1)):
        c_x, c_y = _film_compression_per_unit_load(lam, iface, loading, beta)
        assert c_x == pytest.approx(-s, abs=1e-12)
        assert c_y == pytest.approx(-s * beta, abs=1e-12)


def test_lateral_tension_compresses_a_zero_degree_film_more_along_x():
    """Poisson's effect: under lateral tension the laminate contracts along
    x, and a 0 deg surface ply, stiff along x, takes more axial
    compression. This is why lateral tension can lower the growth stress of
    a quasi-isotropic laminate (docs/physics_models.md, Biaxial loading)."""
    lam = Laminate(_IM7, QI8, 0.152)
    c0, _ = _film_compression_per_unit_load(lam, 6, "compression", 1e-300)
    c_x, c_y = _film_compression_per_unit_load(lam, 6, "compression", -1.0)
    assert c_x > 1.2 * c0 and c_y < 0


# --- closed forms -----------------------------------------------------------


def test_isotropic_clamped_circle_under_equibiaxial_compression():
    """Through the whole laminate path: a one-ply isotropic film over a
    circle of radius R buckles at ``j_{1,1}^2 D / R^2`` under equibiaxial
    compression (Timoshenko & Gere 1961), and at the uniaxial ``27.22 D /
    R^2`` at beta -> 0."""
    lam = _iso_lam()
    for R in (3.0, 12.0):
        e = DelaminationEllipse(0, (0.0, 0.0), R, R, 20.0)
        equi = sublaminate_buckling_load(lam, e, transverse_stress_ratio=1.0)
        assert equi * R**2 / _D_PLY == pytest.approx(_J11**2, rel=1e-6)
        uni = sublaminate_buckling_load(lam, e, transverse_stress_ratio=1e-12)
        assert uni == pytest.approx(sublaminate_buckling_load(lam, e), rel=1e-9)
        assert uni * R**2 / _D_PLY == pytest.approx(27.22, abs=0.01)


@pytest.mark.parametrize("theta", [0.0, 30.0, 90.0])
def test_isotropic_ellipse_under_equibiaxial_compression(theta):
    """Hydrostatic compression of an isotropic film: the onset does not
    depend on the ellipse's orientation, scales as 1/size^2, and lies
    between the bounds of the Rayleigh quotient, ``1 / (1/N_x + 1/N_y) <=
    N_equi <= min(N_x, N_y)`` with ``N_x``, ``N_y`` the uniaxial loads
    along the two axes."""
    lam = _iso_lam()
    e = DelaminationEllipse(0, (0.0, 0.0), 12.0, 5.0, theta)
    N = sublaminate_buckling_load(lam, e, transverse_stress_ratio=1.0)
    ref = sublaminate_buckling_load(
        lam, DelaminationEllipse(0, (0.0, 0.0), 12.0, 5.0, 0.0), transverse_stress_ratio=1.0
    )
    assert N == pytest.approx(ref, rel=1e-8)
    big = DelaminationEllipse(0, (0.0, 0.0), 24.0, 10.0, theta)
    assert sublaminate_buckling_load(lam, big, transverse_stress_ratio=1.0) == pytest.approx(
        N / 4, rel=1e-8
    )
    D_star, _, a, b, (n_xx, n_yy, n_xy) = _sublaminate_film(lam, e)
    N_x = _clamped_ellipse_load_factor(D_star, a, b, n_xx, n_yy, n_xy)
    N_y = _clamped_ellipse_load_factor(D_star, a, b, n_yy, n_xx, -n_xy)
    assert 1 / (1 / N_x + 1 / N_y) < N < min(N_x, N_y)


def test_isotropic_circle_grows_at_the_exact_blister_stress():
    """Equibiaxial compression of a one-ply isotropic film over a circle is
    the axisymmetric circular blister: at the growth stress the exact
    energy release rate (Evans & Hutchinson 1984, solved as a boundary value
    problem in test_blister_energy_release) is the toughness."""
    lam = _iso_lam()
    R = 8.0
    g = sublaminate_growth_stress(
        lam, DelaminationEllipse(0, (0.0, 0.0), R, R, 0.0), transverse_stress_ratio=1.0
    )
    sigma_c = _J11**2 * _D_PLY / R**2 / _T
    assert g.film_buckling_stress_MPa == pytest.approx(sigma_c, rel=1e-6)
    assert g.growth_stress_MPa == pytest.approx(g.film_growth_stress_MPa, rel=1e-12)  # k = 1
    ratio = g.film_growth_stress_MPa / sigma_c
    assert 3.0 < ratio < 20.0
    (G_over_G0,) = _exact_circular_blister(_NU, [ratio])
    G = G_over_G0 * (1 - _NU) * _T * g.film_growth_stress_MPa**2 / _E
    assert G == pytest.approx(mixed_mode_toughness(_ISO), rel=2e-3)


def test_isotropic_square_panel_under_equibiaxial_compression():
    """Simply supported isotropic square: ``N_cr = 2 pi^2 D / a^2`` under
    equibiaxial compression, half the uniaxial ``4 pi^2 D / a^2``
    (Timoshenko & Gere 1961 §9.2)."""
    lam = _iso_lam(8)
    D = _E * lam.thickness_mm**3 / (12 * (1 - _NU**2))
    a = 120.0
    uni = panel_buckling_load(lam, a, a)
    assert uni == pytest.approx(4 * math.pi**2 * D / a**2, rel=1e-8)
    equi = panel_buckling_load(lam, a, a, transverse_stress_ratio=1.0)
    assert equi == pytest.approx(2 * math.pi**2 * D / a**2, rel=1e-8)
    assert panel_buckling_load(lam, a, a, transverse_stress_ratio=1e-12) == pytest.approx(
        uni, rel=1e-9
    )
    assert panel_buckling_load(lam, a, a, "clamped", transverse_stress_ratio=1.0) == (
        pytest.approx(1.9 * equi, rel=1e-12)
    )


def _navier_biaxial(D, a, b, n_x, n_y, n_max=10):
    """Navier SSSS biaxial (D16 = D26 = 0): the lowest load factor over the
    positive denominators, and its (m, n)."""
    best = (math.inf, None)
    for m in range(1, n_max + 1):
        for n in range(1, n_max + 1):
            den = n_x * m**2 / a**2 + n_y * n**2 / b**2
            if den <= 0:
                continue
            num = math.pi**2 * (
                D[0, 0] * m**4 / a**4
                + 2 * (D[0, 1] + 2 * D[2, 2]) * m**2 * n**2 / (a * b) ** 2
                + D[1, 1] * n**4 / b**4
            )
            best = min(best, (num / den, (m, n)), key=lambda t: t[0])
    return best


@pytest.mark.parametrize("a, b", [(150.0, 100.0), (100.0, 150.0), (300.0, 80.0)])
@pytest.mark.parametrize("n_x, n_y", [(1, 1), (1, 0.5), (1, -0.5), (1, -1), (-1, 0.7), (0.3, 1)])
def test_panel_ritz_is_the_biaxial_navier_value_for_a_cross_ply(a, b, n_x, n_y):
    lam = Laminate(_IM7, CP8, 0.152)
    D = lam.abd_matrices()[2]
    assert abs(D[0, 2]) < 1e-9 and abs(D[1, 2]) < 1e-9
    navier, (m, n) = _navier_biaxial(D, a, b, n_x, n_y)
    assert m < 10 and n < 10  # the closed form's range is not what limits it
    assert _ssss_ritz_load_factor(D, a, b, n_x, n_y) == pytest.approx(navier, rel=1e-7)


def test_panel_buckling_along_the_ray():
    lam = Laminate(_IM7, QI8, 0.152)
    values = [
        panel_buckling_load(lam, 150.0, 100.0, transverse_stress_ratio=b)
        for b in (-1.0, -0.5, 0.0, 0.5, 1.0)
    ]
    assert all(v1 > v2 for v1, v2 in zip(values, values[1:]))  # lateral compression lowers it
    for beta in (0.0, 0.5, 1.0):
        assert panel_buckling_load(
            lam, 150.0, 100.0, loading="tension", transverse_stress_ratio=beta
        ) == float("inf")
    # Axial tension with lateral compression buckles the panel across.
    assert math.isfinite(
        panel_buckling_load(lam, 150.0, 100.0, loading="tension", transverse_stress_ratio=-1.0)
    )


# --- stabilisation, net tension ---------------------------------------------


def test_lateral_tension_raises_the_film_onset_and_growth_load():
    """At a fixed axial film compression, lateral film tension removes work
    from the Rayleigh quotient, so the onset rises (strictly here); the
    post-buckled film also needs more axial load to grow, and lateral
    compression less. A 2-ply [45/0] film over a 15 x 6 mm ellipse at 20
    deg."""
    lam = Laminate(_IM7, QI8, 0.152)
    e = DelaminationEllipse(5, (0.0, 0.0), 15.0, 6.0, 20.0)
    D_star, A, a, b, (n_xx, n_yy, n_xy) = _sublaminate_film(lam, e)
    y = (n_yy, n_xx, -n_xy)
    G_c = mixed_mode_toughness(_IM7)
    onset, growth = [], []
    for t in (0.5, 0.25, 0.0, -0.1, -0.25, -0.5):  # lateral compression to tension
        n = (n_xx + t * y[0], n_yy + t * y[1], n_xy + t * y[2])
        onset.append(_clamped_ellipse_load_factor(D_star, a, b, *n))
        growth.append(_PostbuckledFilm(D_star, A, a, b, *n).growth_load(G_c)[0])
    assert all(v1 > v0 for v0, v1 in zip(onset, onset[1:]))
    assert all(v1 > v0 for v0, v1 in zip(growth, growth[1:]))


def test_isotropic_laminate_growth_rises_with_lateral_tension():
    """Without Poisson's effect on the film's axial stress (an isotropic
    laminate), the growth stress rises monotonically from biaxial
    compression to lateral tension."""
    lam = _iso_lam()
    e = DelaminationEllipse(0, (0.0, 0.0), 6.0, 3.0, 30.0)
    values = [
        sublaminate_growth_stress(lam, e, transverse_stress_ratio=b).growth_stress_MPa
        for b in (1.0, 0.5, 0.25, 1e-9, -0.25, -0.5, -1.0)
    ]
    assert all(v1 > v0 for v0, v1 in zip(values, values[1:]))
    assert values[3] == pytest.approx(sublaminate_growth_stress(lam, e).growth_stress_MPa, rel=1e-6)


@pytest.mark.parametrize("beta", [0.0, 1e-9, 0.5, 1.0])
def test_a_film_in_tension_along_both_axes_does_not_buckle(beta):
    lam = Laminate(_IM7, QI8, 0.152)
    damage = _damage(QI8, 10.0)
    for e in _largest(damage).values():
        kw = dict(loading="tension", transverse_stress_ratio=beta)
        assert sublaminate_buckling_load(lam, e, **kw) == math.inf
        g = sublaminate_growth_stress(lam, e, **kw)
        assert g.growth_stress_MPa == math.inf and g.growth_axis is None
    g = weakest_sublaminate_growth(lam, damage, loading="tension", transverse_stress_ratio=beta)
    assert g.growth_stress_MPa == math.inf
    D_star, _, a, b, _ = _sublaminate_film(lam, _largest(damage)[5])
    assert _clamped_ellipse_load_factor(D_star, a, b, -1.0, -0.3, 0.2) == math.inf


@pytest.mark.parametrize("beta", [0.5, 1.0])
def test_biaxial_tension_is_the_hole_model(beta):
    s = _run(_cfg(CP8, 10.0, "tension", beta))
    e = _run(_cfg(CP8, 10.0, "tension", beta, tier="empirical"))
    assert s.residual_strength_MPa == e.residual_strength_MPa
    assert s.critical_sublaminate is None and s.buckling_eigenvalues is None
    assert "no delamination grows" in s.notes[0]


# --- direction, rotation ----------------------------------------------------


@pytest.mark.parametrize("theta", [0.0, 25.0])
def test_quarter_turn_of_the_load_and_the_ellipse_gives_the_same_growth(theta):
    """Isotropic film: axial tension with lateral compression (tension, beta
    = -1) on an ellipse at theta is axial compression with lateral tension
    (compression, beta = -1) on the ellipse at theta + 90 deg, turned by
    90 deg. The film stress axis turns with it; the growth axis does not."""
    lam = _iso_lam()
    e = DelaminationEllipse(0, (0.0, 0.0), 6.0, 3.0, theta)
    turned = replace(e, orientation_deg=theta + 90.0)
    g_t = sublaminate_growth_stress(lam, e, loading="tension", transverse_stress_ratio=-1.0)
    g_c = sublaminate_growth_stress(
        lam, turned, loading="compression", transverse_stress_ratio=-1.0
    )
    assert g_t.film_stress_axis == "y" and g_c.film_stress_axis == "x"
    assert g_t.film_stress_ratio == pytest.approx(-1.0, abs=1e-12)
    assert g_c.film_stress_ratio == pytest.approx(-1.0, abs=1e-12)
    assert g_t.growth_axis == g_c.growth_axis
    assert g_t.buckling_load_N_per_mm == pytest.approx(g_c.buckling_load_N_per_mm, rel=1e-8)
    assert g_t.growth_stress_MPa == pytest.approx(g_c.growth_stress_MPa, rel=1e-6)


def test_tension_with_lateral_compression_grows_the_cross_ply_delamination():
    """D10: under axial tension a laterally compressed film buckles along y
    and grows; on the 8-ply cross-ply at 10 J it governs over the hole
    model."""
    lam = Laminate(_IM7, CP8, 0.152)
    g = weakest_sublaminate_growth(
        lam, _damage(CP8, 10.0), loading="tension", transverse_stress_ratio=-1.0
    )
    assert g.film_stress_axis == "y" and g.film_stress_ratio < 0
    assert g.interface_index == 4 and g.growth_axis == "major"
    assert g.growth_stress_MPa == pytest.approx(578.4, abs=0.2)
    r = _run(_cfg(CP8, 10.0, "tension", -1.0))
    assert r.residual_strength_MPa == pytest.approx(g.growth_stress_MPa, rel=1e-9)
    assert r.critical_sublaminate == 4
    assert r.buckling_eigenvalues == [pytest.approx(g.buckling_load_N_per_mm, rel=1e-12)]
    note = r.notes[0]
    assert "Biaxial sigma_y/sigma_x = -1.00: experimental, not validated" in note
    assert "interface 4 grows" in note and "more compressed along y" in note
    assert "whole panel buckles" in note
    assert any(n.startswith("TAI: equivalent hole") for n in r.notes)


# --- continuity, pruning ----------------------------------------------------


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_results_are_continuous_at_zero_ratio(loading):
    base = _run(_cfg(QI8, 10.0, loading, 0.0))
    for beta in (1e-9, -1e-9):
        r = _run(_cfg(QI8, 10.0, loading, beta))
        assert r.knockdown == pytest.approx(base.knockdown, rel=1e-6)
        assert r.biaxial_interaction_factor == pytest.approx(1.0, rel=1e-6)
    lam = Laminate(_IM7, QI8, 0.152)
    damage = _damage(QI8, 10.0)
    g0 = weakest_sublaminate_growth(lam, damage)
    for beta in (1e-9, -1e-9):
        g = weakest_sublaminate_growth(lam, damage, transverse_stress_ratio=beta)
        assert g.interface_index == g0.interface_index and g.growth_axis == g0.growth_axis
        assert g.growth_stress_MPa == pytest.approx(g0.growth_stress_MPa, rel=1e-6)
        assert g.buckling_load_N_per_mm == pytest.approx(g0.buckling_load_N_per_mm, rel=1e-6)


def _check_pruned_is_full_minimum(layup, energy, loading, beta):
    lam = Laminate(_IM7, layup, 0.152)
    damage = _damage(layup, energy)
    pruned = weakest_sublaminate_growth(lam, damage, loading=loading, transverse_stress_ratio=beta)
    full = {
        i: sublaminate_growth_stress(lam, e, loading=loading, transverse_stress_ratio=beta)
        for i, e in _largest(damage).items()
    }
    # The skip test is a lower bound: no delamination grows before its
    # sublaminate buckles.
    for i, e in _largest(damage).items():
        terms = sa._growth_terms(lam, e, "simply_supported", loading, beta)
        assert full[i].growth_stress_MPa >= sa._onset_stress(terms)
    best = min(full.values(), key=lambda g: (g.growth_stress_MPa, g.interface_index))
    if math.isinf(best.growth_stress_MPa):
        assert math.isinf(pruned.growth_stress_MPa)
        return
    # Skipped interfaces only: the others are solved in full, not capped
    # (under a mixed load a capped solve can jump to another post-buckled
    # branch; 8-ply QI, 20 J, compression, beta = -1 found one), so the
    # result is the full minimum exactly.
    assert pruned == best


@pytest.mark.parametrize(
    "layup, energy, loading, beta",
    [(QI8, 10.0, "compression", 1.0), (CP8, 10.0, "tension", -1.0)],
)
def test_pruned_search_finds_the_full_minimum(layup, energy, loading, beta):
    _check_pruned_is_full_minimum(layup, energy, loading, beta)


@pytest.mark.slow
@pytest.mark.parametrize("layup", [QI8, CP8, [0, 45, -45, 90] * 2 + [90, -45, 45, 0] * 2])
@pytest.mark.parametrize("energy", [5.0, 10.0, 20.0])
@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", [-1.0, -0.5, -0.05, 0.5, 1.0])
def test_pruned_search_finds_the_full_minimum_grid(layup, energy, loading, beta):
    _check_pruned_is_full_minimum(layup, energy, loading, beta)


# --- pinned knockdowns ------------------------------------------------------

# 150 x 100 mm IM7/8552 panels, 0.152 mm plies, 16 mm hemisphere, 5.5 kg.
# Unvalidated: pinned so that a change is seen, not because they are right.
_PINNED = [
    # layup, J, loading, beta, knockdown, governs
    (QI8, 5.0, "compression", -1.0, 0.5459, "bound"),
    (QI8, 5.0, "compression", 0.5, 0.4803, "growth"),
    (QI8, 5.0, "compression", 1.0, 0.3839, "growth"),
    (QI8, 10.0, "compression", 1.0, 0.3937, "growth"),
    (CP8, 5.0, "compression", 1.0, 0.4853, "growth"),
    (CP8, 5.0, "tension", -1.0, 0.8397, "growth"),
    (CP8, 10.0, "tension", -1.0, 0.7035, "growth"),
    (CP8, 10.0, "tension", -0.5, 0.8313, "bound"),
    (QI8, 10.0, "tension", -1.0, 0.8409, "bound"),
]


@pytest.mark.parametrize("layup, energy, loading, beta, knockdown, governs", _PINNED)
def test_pinned_knockdowns(layup, energy, loading, beta, knockdown, governs):
    s = _run(_cfg(layup, energy, loading, beta))
    e = _run(_cfg(layup, energy, loading, beta, tier="empirical"))
    assert s.pristine_strength_MPa == e.pristine_strength_MPa
    assert s.knockdown == pytest.approx(knockdown, abs=5e-4)
    if governs == "growth":
        assert s.knockdown < e.knockdown - 1e-3
    else:
        assert s.knockdown == e.knockdown
    assert s.residual_transverse_stress_MPa == pytest.approx(
        (1 if loading == "tension" else -1) * beta * s.residual_strength_MPa, rel=1e-12
    )
    uniaxial = _run(_cfg(layup, energy, loading, 0.0))
    assert s.biaxial_interaction_factor == pytest.approx(
        s.residual_strength_MPa / uniaxial.residual_strength_MPa, rel=1e-12
    )
    assert "semi_analytical puts the film's lateral stress" in s.notes[0]


def test_biaxial_compression_lowers_the_growth_stress():
    """Lateral compression destabilises the film: on the 8-ply QI at 5 J
    the growth stress falls from beta = 0 to beta = 1, and growth then
    governs over the homothetic Soutis bound."""
    lam = Laminate(_IM7, QI8, 0.152)
    damage = _damage(QI8, 5.0)
    values = [
        weakest_sublaminate_growth(lam, damage, transverse_stress_ratio=b).growth_stress_MPa
        for b in (1e-12, 0.5, 1.0)
    ]
    assert all(v1 < v0 for v0, v1 in zip(values, values[1:]))
    assert np.isclose(values[0], weakest_sublaminate_growth(lam, damage).growth_stress_MPa)
