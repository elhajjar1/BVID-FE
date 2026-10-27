"""Fibre-damage open-hole TAI model (``fibre_damage_tai``).

The equivalent hole is ``tai_notch_fraction`` x the delamination footprint
width across the load, with the Whitney-Nuismer point-stress strength and
Tan's (1988) orthotropic finite-width correction. Checked here: the
analytical limits (no damage -> pristine, a vanishing hole -> knockdown 1,
the finite-width factor -> 1 for a narrow hole and -> 0 as the hole fills
the panel, the isotropic closed form at Kt_inf = 3), pinned values, the
energy trend on IM7/8552 panels, and the TAI note on the results.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY, OrthotropicMaterial
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.soutis_openhole import (
    fibre_damage_tai,
    finite_width_factor,
    lekhnitskii_kt_infinity,
    tai_hole_radius_mm,
    transverse_damage_width_mm,
)

IM7 = MATERIAL_LIBRARY["IM7/8552"]
QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
QI16 = [0, 45, -45, 90] * 2 + [90, -45, 45, 0] * 2
CP8 = [0, 90, 0, 90, 90, 0, 90, 0]


def _damage(width_mm: float) -> DamageState:
    """One delamination whose footprint is ``width_mm`` wide across the load."""
    return DamageState([DelaminationEllipse(0, (75.0, 50.0), width_mm / 2, 1.0, 90.0)])


def test_no_damage_returns_pristine():
    assert fibre_damage_tai(IM7, DamageState([]), 900.0, 3.0, 100.0) == 900.0


@pytest.mark.parametrize("Kt_inf", [3.0, 5.18])
def test_vanishing_hole_gives_unit_knockdown(Kt_inf):
    sigma = fibre_damage_tai(IM7, _damage(1e-4), 900.0, Kt_inf, 100.0)
    assert sigma == pytest.approx(900.0, rel=1e-9)
    assert sigma <= 900.0


def test_point_stress_value_pinned():
    """40 mm wide footprint -> 2 mm hole (xi = 0.5): 2 / (2 + 1/4 + 3/16),
    times the isotropic finite-width factor at 2R/W = 0.02."""
    d = _damage(40.0)
    assert tai_hole_radius_mm(IM7, d, 100.0) == pytest.approx(1.0)
    kd_inf = 2.0 / (2.0 + 0.25 + 3.0 / 16.0)
    fwc = 3 * 0.98 / (2 + 0.98**3)
    assert fibre_damage_tai(IM7, d, 1.0, 3.0, 100.0) == pytest.approx(kd_inf * fwc, rel=1e-12)
    assert fibre_damage_tai(IM7, d, 1.0, 3.0, None) == pytest.approx(kd_inf, rel=1e-12)
    assert fibre_damage_tai(IM7, d, 1000.0, 3.0, 100.0) == pytest.approx(820.2, abs=0.05)


def test_finite_width_factor_limits():
    for Kt in (2.0, 3.0, 5.18, 7.8):
        assert finite_width_factor(0.0, Kt) == 1.0
        assert finite_width_factor(1.0, Kt) == 0.0
        assert finite_width_factor(1e-6, Kt) == pytest.approx(1.0, abs=1e-9)
        assert finite_width_factor(1.0 - 1e-9, Kt) == pytest.approx(0.0, abs=1e-6)
        series = [finite_width_factor(x / 20, Kt) for x in range(21)]
        assert all(a > b for a, b in zip(series, series[1:])), series


@pytest.mark.parametrize("x", [0.05, 0.2, 0.5, 0.8])
def test_finite_width_factor_isotropic_closed_form(x):
    """At Kt_inf = 3 Tan's factor is 1 / {[2 + (1 - x)^3] / [3 (1 - x)]}."""
    assert finite_width_factor(x, 3.0) == pytest.approx(3 * (1 - x) / (2 + (1 - x) ** 3))


def test_finite_width_factor_orthotropic_term():
    """Tan (1988): a higher Kt_inf lowers K_T / K_T_inf, i.e. raises the factor."""
    assert finite_width_factor(0.5, 5.18) > finite_width_factor(0.5, 3.0)
    assert finite_width_factor(0.5, 5.18) == pytest.approx(0.70588 + 0.02093 * 2.18, abs=2e-4)


def test_narrow_panel_lowers_strength():
    m = dataclasses.replace(IM7, tai_notch_fraction=0.5)
    d = _damage(40.0)  # 20 mm hole
    wide = fibre_damage_tai(m, d, 1.0, 3.0, 1e6)
    narrow = fibre_damage_tai(m, d, 1.0, 3.0, 50.0)
    assert fibre_damage_tai(m, d, 1.0, 3.0, None) == pytest.approx(wide, rel=1e-9)
    assert narrow == pytest.approx(wide * finite_width_factor(20.0 / 50.0, 3.0))
    assert narrow < wide


def test_transverse_width_rotated_offset_and_capped():
    e0 = DelaminationEllipse(0, (75.0, 50.0), 20.0, 10.0, 0.0)  # major along the load
    e90 = DelaminationEllipse(1, (75.0, 50.0), 20.0, 10.0, 90.0)
    e45 = DelaminationEllipse(2, (75.0, 50.0), 20.0, 10.0, 45.0)
    assert transverse_damage_width_mm(DamageState([e0])) == pytest.approx(20.0)
    assert transverse_damage_width_mm(DamageState([e90])) == pytest.approx(40.0)
    assert transverse_damage_width_mm(DamageState([e45])) == pytest.approx(2 * math.sqrt(250.0))
    shifted = DelaminationEllipse(1, (75.0, 60.0), 20.0, 10.0, 0.0)
    assert transverse_damage_width_mm(DamageState([e0, shifted])) == pytest.approx(30.0)
    assert transverse_damage_width_mm(DamageState([e90]), panel_width_mm=25.0) == 25.0
    assert transverse_damage_width_mm(DamageState([])) == 0.0


def test_fibre_break_core_sets_a_minimum_hole():
    d = dataclasses.replace(_damage(40.0), fiber_break_radius_mm=3.0)
    assert tai_hole_radius_mm(IM7, d, 100.0) == 3.0
    assert fibre_damage_tai(IM7, d, 1.0, 3.0, 100.0) < fibre_damage_tai(
        IM7, _damage(40.0), 1.0, 3.0, 100.0
    )


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5])
def test_tai_notch_fraction_validated(bad):
    with pytest.raises(ValueError, match="tai_notch_fraction"):
        dataclasses.replace(IM7, tai_notch_fraction=bad)


def test_presets_carry_the_fitted_fraction():
    assert OrthotropicMaterial.__dataclass_fields__["tai_notch_fraction"].default == 0.05
    assert all(m.tai_notch_fraction == 0.05 for m in MATERIAL_LIBRARY.values())


def _run(layup, energy_J, loading, tier="empirical"):
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=layup,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        impact=ImpactEvent(energy_J=energy_J),
        loading=loading,
        tier=tier,
    )
    return BvidAnalysis(cfg).run()


@pytest.mark.parametrize("layup", [QI8, QI16, CP8], ids=["QI8", "QI16", "CP8"])
def test_tai_monotonic_in_energy_and_above_cai(layup):
    """The knockdown falls with impact energy (5-40 J) and, unlike the
    whole-footprint hole, stays above the CAI knockdown once damaged."""
    energies = (5.0, 10.0, 15.0, 20.0, 30.0, 40.0)
    tai = [_run(layup, E, "tension").knockdown for E in energies]
    assert all(a >= b for a, b in zip(tai, tai[1:])), tai
    assert tai[-1] < tai[0]
    for E, kd_t in zip(energies, tai):
        kd_c = _run(layup, E, "compression").knockdown
        if kd_c < 1.0:
            assert kd_t > kd_c, (E, kd_t, kd_c)


def test_tai_note_on_closed_form_tiers_only():
    for tier in ("empirical", "semi_analytical"):
        notes = _run(QI8, 15.0, "tension", tier).notes
        assert any(n.startswith("TAI: equivalent hole") and "not validated" in n for n in notes)
    assert not any(n.startswith("TAI:") for n in _run(QI8, 15.0, "compression").notes)


def test_bvid_analysis_matches_the_function():
    r = _run(CP8, 15.0, "tension")
    lam = Laminate(IM7, CP8, 0.152)
    expected = fibre_damage_tai(
        IM7, r.damage, r.pristine_strength_MPa, lekhnitskii_kt_infinity(lam), 100.0
    )
    assert r.residual_strength_MPa == pytest.approx(expected, rel=1e-12)
