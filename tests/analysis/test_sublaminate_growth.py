"""Compression strength from buckling-driven delamination growth.

The delaminated sublaminate buckling at all does not fail the laminate: a
thin sublaminate over a 30 mm delamination buckles at a few MPa, so taking
buckling onset as the compressive strength gave CAI knockdowns of 0.002-0.04
against measured values of 0.3-0.6. The laminate fails when the buckled
sublaminate releases enough energy to grow the delamination, G = G_c, with
G_c the Benzeggagh-Kenane mix of G_Ic and G_IIc at the thin-film phase angle
(52.1 deg, eta = 2). G is that of the post-buckled elliptical film, a 2D
blister (``test_blister_energy_release.py``); until the 2D blister it was the
straight-sided strip's (Chai, Babcock & Knauss 1981; Hutchinson & Suo 1992),

    G = h / (2 E_f) * (sigma - sigma_c) * (sigma + 3 sigma_c),

with sigma the film stress, sigma_c the film buckling stress, h the
sublaminate thickness and E_f the film modulus, kept as
``_strip_release_rate`` for reference.
"""

from __future__ import annotations

import math

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.analysis.semi_analytical import (
    SublaminateGrowth,
    _strip_growth_film_stress,
    _strip_release_rate,
    mixed_mode_toughness,
    semi_analytical_cai,
    sublaminate_buckling_load,
    sublaminate_growth_stress,
    weakest_sublaminate_growth,
)
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.soutis_openhole import soutis_cai
from bvidfe.impact.mapping import ImpactEvent

_M = MATERIAL_LIBRARY["IM7/8552"]
_QI16 = [0, 45, -45, 90] * 4
_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]


def _film_terms(layup, interface_index, ply_t=0.152):
    """Film thickness and modulus for the (upper, thinner) sublaminate,
    computed independently of the implementation: E_f is the film's x-stress
    per unit laminate x-strain with the laminate's Poisson contraction."""
    lam = Laminate(_M, layup, ply_t)
    sub = layup[: interface_index + 1]
    h = ply_t * len(sub)
    A_sub, _, _ = Laminate(_M, sub, ply_t).abd_matrices()
    Ex, _, _, nu = lam.effective_engineering_constants()
    E_f = (A_sub[0, 0] - nu * A_sub[0, 1]) / h
    return lam, h, E_f, Ex


def test_mixed_mode_toughness_im7_8552():
    expected = _M.G_Ic + (_M.G_IIc - _M.G_Ic) * math.sin(math.radians(52.1)) ** 4
    assert mixed_mode_toughness(_M) == pytest.approx(expected, rel=1e-12)
    assert mixed_mode_toughness(_M) == pytest.approx(0.478, abs=1e-3)


@pytest.mark.parametrize("name", sorted(MATERIAL_LIBRARY))
def test_mixed_mode_toughness_between_mode_i_and_mode_ii(name):
    m = MATERIAL_LIBRARY[name]
    assert m.G_Ic < mixed_mode_toughness(m) < m.G_IIc


def test_strip_reference_releases_exactly_the_mixed_mode_toughness():
    """The straight-sided strip (the growth law before the 2D blister)."""
    _, h, E_f, _ = _film_terms(_QI16, 3)
    sigma_c = 120.0
    sigma = _strip_growth_film_stress(sigma_c, E_f, mixed_mode_toughness(_M), h)
    G = h / (2.0 * E_f) * (sigma - sigma_c) * (sigma + 3.0 * sigma_c)
    assert G == pytest.approx(mixed_mode_toughness(_M), rel=1e-12)
    assert _strip_release_rate(sigma, sigma_c, E_f, h) == pytest.approx(G, rel=1e-12)
    assert _strip_release_rate(sigma_c, sigma_c, E_f, h) == 0.0
    assert sigma > sigma_c  # growth needs a buckled film


def test_2d_growth_stress_is_above_the_strips_on_this_sublaminate():
    """The 2D blister releases less energy than the strip at the same
    buckling stress, so it grows at a higher film stress: 417.3 against
    272.3 MPa (1.53x) for this 4-ply film over a near-circular 20 x 15 mm
    ellipse, which buckles at 106.8 MPa and grows across the load."""
    ell = DelaminationEllipse(3, (0.0, 0.0), 20.0, 15.0, 0.0)
    lam, h, E_f, _ = _film_terms(_QI16, 3)
    g = sublaminate_growth_stress(lam, ell)
    strip = _strip_growth_film_stress(g.film_buckling_stress_MPa, E_f, mixed_mode_toughness(_M), h)
    assert g.film_growth_stress_MPa > g.film_buckling_stress_MPa
    assert g.film_growth_stress_MPa == pytest.approx(417.32, abs=0.01)
    assert g.film_growth_stress_MPa / strip == pytest.approx(1.532, abs=1e-3)
    assert g.growth_axis == "minor"


def test_film_buckling_stress_is_sublaminate_buckling_load_over_thickness():
    ell = DelaminationEllipse(3, (0.0, 0.0), 20.0, 15.0, 0.0)
    lam, h, _, _ = _film_terms(_QI16, 3)
    g = sublaminate_growth_stress(lam, ell, boundary="clamped")
    N_cr = sublaminate_buckling_load(lam, ell, boundary="clamped")
    assert g.buckling_load_N_per_mm == pytest.approx(N_cr, rel=1e-12)
    assert g.film_buckling_stress_MPa == pytest.approx(N_cr / h, rel=1e-12)


def test_growth_stress_converts_film_stress_to_laminate_stress():
    """Film and laminate share the far-field strain, so the laminate stress
    at growth is Ex * (film stress / E_f)."""
    ell = DelaminationEllipse(3, (0.0, 0.0), 20.0, 15.0, 0.0)
    lam, _, E_f, Ex = _film_terms(_QI16, 3)
    g = sublaminate_growth_stress(lam, ell)
    assert g.growth_stress_MPa == pytest.approx(Ex * g.film_growth_stress_MPa / E_f, rel=1e-9)


def test_large_delamination_reaches_a_steady_state_growth_stress():
    """When the film buckles at a negligible stress, the growth stress no
    longer depends on the delamination's size. For the strip it is the
    stress whose stored film energy equals G_c, sqrt(2 E_f G_c / h); the 2D
    blister's (continued past 100x onset with the strip's shape) is 1.08x
    that here."""
    lam, h, E_f, _ = _film_terms(_QI16, 3)
    steady = math.sqrt(2.0 * E_f * mixed_mode_toughness(_M) / h)
    assert _strip_growth_film_stress(1e-6 * steady, E_f, mixed_mode_toughness(_M), h) == (
        pytest.approx(steady, rel=1e-5)
    )
    g = sublaminate_growth_stress(lam, DelaminationEllipse(3, (0.0, 0.0), 4000.0, 3000.0, 0.0))
    g2 = sublaminate_growth_stress(lam, DelaminationEllipse(3, (0.0, 0.0), 8000.0, 6000.0, 0.0))
    assert g.film_buckling_stress_MPa < 1e-3 * steady
    assert g2.film_growth_stress_MPa == pytest.approx(g.film_growth_stress_MPa, rel=1e-3)
    assert g.film_growth_stress_MPa / steady == pytest.approx(1.08, abs=0.01)


def test_weakest_growth_is_the_minimum_over_interfaces():
    lam = Laminate(_M, _QI16, 0.152)
    small_at_3 = DelaminationEllipse(3, (0.0, 0.0), 5.0, 4.0, 0.0)
    large_at_3 = DelaminationEllipse(3, (0.0, 0.0), 20.0, 15.0, 0.0)
    at_7 = DelaminationEllipse(7, (0.0, 0.0), 25.0, 20.0, 0.0)
    at_14 = DelaminationEllipse(14, (0.0, 0.0), 30.0, 10.0, 0.0)
    damage = DamageState([small_at_3, large_at_3, at_7, at_14], dent_depth_mm=0.3)
    weakest = weakest_sublaminate_growth(lam, damage)
    assert isinstance(weakest, SublaminateGrowth)
    # The largest ellipse at each interface is the one that grows. The
    # interfaces that cannot be the weakest are not solved, and the others
    # stop at the lowest stress found so far, so the result agrees with the
    # full minimum to the root finder's tolerance.
    per_interface = [sublaminate_growth_stress(lam, e) for e in (large_at_3, at_7, at_14)]
    expected = min(per_interface, key=lambda g: g.growth_stress_MPa)
    assert weakest.interface_index == expected.interface_index
    assert weakest.growth_stress_MPa == pytest.approx(expected.growth_stress_MPa, rel=1e-6)
    assert weakest.growth_axis == expected.growth_axis


def test_weakest_growth_is_none_without_delaminations():
    lam = Laminate(_M, _QI16, 0.152)
    assert weakest_sublaminate_growth(lam, DamageState([], dent_depth_mm=0.0)) is None


def test_semi_analytical_cai_is_min_of_soutis_and_growth():
    lam = Laminate(_M, _QI16, 0.152)
    damage = DamageState(
        [
            DelaminationEllipse(3, (0.0, 0.0), 20.0, 15.0, 0.0),
            DelaminationEllipse(7, (0.0, 0.0), 25.0, 20.0, 0.0),
        ],
        dent_depth_mm=0.3,
    )
    result = semi_analytical_cai(lam, damage, 500.0, 15000.0)
    weakest = weakest_sublaminate_growth(lam, damage)
    soutis = soutis_cai(_M, damage.projected_damage_area_mm2, 15000.0, 500.0)
    assert result.residual_strength_MPa == pytest.approx(
        min(soutis, weakest.growth_stress_MPa), rel=1e-12
    )
    assert result.critical_interface_index == weakest.interface_index
    assert result.critical_buckling_load_N == pytest.approx(weakest.buckling_load_N_per_mm)


@pytest.mark.parametrize("layup", [_QI8, _QI16], ids=["8-ply", "16-ply"])
def test_semi_analytical_cai_knockdown_is_not_collapsed_by_buckling_onset(layup):
    """Buckling onset gave knockdowns of 0.002-0.04 here; delamination growth
    gives values in the measured CAI range."""
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=layup,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        loading="compression",
        tier="semi_analytical",
        impact=ImpactEvent(15.0, ImpactorGeometry(diameter_mm=16.0), mass_kg=5.5),
    )
    r = BvidAnalysis(cfg).run()
    assert 0.15 < r.knockdown < 0.8, r.knockdown


def _fe3d_cfg(damage: DamageState, Lx: float = 200.0, Ly: float = 150.0) -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(Lx, Ly),
        loading="compression",
        tier="fe3d",
        damage=damage,
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0),
    )


def test_fe3d_compression_knockdown_excludes_whole_panel_buckling():
    """An 8-ply 200x150 mm panel buckles as a whole at ~12.2 MPa. With panel
    buckling in both the damaged and undamaged strengths, any delamination
    that grows above 12.2 MPa was invisible (knockdown 1.0). Panel buckling
    is a test-fixture effect (CAI fixtures carry anti-buckling guides), so
    it is reported in the notes but kept out of the knockdown."""
    damage = DamageState([DelaminationEllipse(3, (100.0, 75.0), 15.0, 10.0, 0.0)], 0.0)
    r = BvidAnalysis(_fe3d_cfg(damage)).run()
    assert r.knockdown < 0.9, r.notes
    assert any("panel buckling" in n.lower() and "excluded" in n.lower() for n in r.notes), r.notes
    assert r.buckling_eigenvalues is not None and r.buckling_eigenvalues[0] > 0
    assert r.warnings == []


def test_fe3d_compression_residual_uses_growth_interface():
    damage = DamageState(
        [
            DelaminationEllipse(1, (100.0, 75.0), 12.0, 8.0, 0.0),
            DelaminationEllipse(3, (100.0, 75.0), 15.0, 10.0, 0.0),
        ],
        0.0,
    )
    r = BvidAnalysis(_fe3d_cfg(damage)).run()
    lam = Laminate(_M, _QI8, 0.152)
    weakest = weakest_sublaminate_growth(lam, damage)
    assert r.critical_sublaminate == weakest.interface_index
    assert any(f"interface {weakest.interface_index}" in n for n in r.notes), r.notes
