import math

import pytest

from bvidfe.analysis.semi_analytical import (
    SemiAnalyticalResult,
    _sublaminate_D_matrix,
    find_critical_interface,
    semi_analytical_cai,
    semi_analytical_tai,
    sublaminate_buckling_load,
)
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse


def test_buckling_load_positive_for_realistic_case():
    m = MATERIAL_LIBRARY["IM7/8552"]
    layup = [0, 45, -45, 90] * 4  # 16 plies
    lam = Laminate(m, layup, 0.152)
    # Delamination at interface 3 (between ply 3 and ply 4)
    ell = DelaminationEllipse(3, (75, 50), major_mm=20, minor_mm=10, orientation_deg=0)
    N = sublaminate_buckling_load(lam, ell)
    assert N > 0


def test_buckling_load_decreases_with_ellipse_size():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ell_small = DelaminationEllipse(3, (0, 0), 10, 6, 0)
    ell_large = DelaminationEllipse(3, (0, 0), 40, 20, 0)
    N_small = sublaminate_buckling_load(lam, ell_small)
    N_large = sublaminate_buckling_load(lam, ell_large)
    assert N_small > N_large  # smaller ellipse => harder to buckle


def test_find_critical_interface_picks_largest_ellipse_near_surface():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState(
        [
            DelaminationEllipse(1, (0, 0), 10, 5, 0),  # near impact face
            DelaminationEllipse(7, (0, 0), 20, 10, 0),  # mid-plane, larger
            DelaminationEllipse(14, (0, 0), 8, 4, 0),  # near back face
        ],
        dent_depth_mm=0.3,
    )
    idx = find_critical_interface(ds, lam)
    # Interface 7 has largest area, but let's check scoring:
    # interface 7 z-upper = 7 * 0.152 = 1.064 (distance from top)
    # interface 7 z-lower = (16-8) * 0.152 = 1.216 (distance from bottom)
    # max |z| = 1.216
    # score_7 = pi*20*10 * 1.216 ≈ 764
    # score_1 = pi*10*5 * (16-2)*0.152 = pi*50 * 2.128 ≈ 334 (distance from bottom larger)
    # score_14 = pi*8*4 * (16-15)*0.152 = small
    # Expect 7 to win
    assert idx == 7


def test_find_critical_interface_returns_none_for_empty():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState([], dent_depth_mm=0.0)
    idx = find_critical_interface(ds, lam)
    assert idx is None


def test_sublaminate_above_interface_uses_correct_plies():
    """Critical sublaminate includes plies above the interface (from top surface to the delamination)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    layup = [0, 45, -45, 90, 90, -45, 45, 0]  # 8 plies
    lam = Laminate(m, layup, 0.152)
    # Delamination at interface 1 (between ply 1 and ply 2)
    # Upper sublaminate = plies 0, 1 (2 plies); thinner sublaminate buckles first.
    ell = DelaminationEllipse(1, (0, 0), 20, 10, 0)
    N = sublaminate_buckling_load(lam, ell)
    assert N > 0


def test_semi_analytical_cai_pristine_returns_pristine():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState([], dent_depth_mm=0.0)
    result = semi_analytical_cai(lam, ds, 500.0, 15000.0)
    assert isinstance(result, SemiAnalyticalResult)
    assert result.residual_strength_MPa == 500.0
    assert result.critical_interface_index is None
    assert result.critical_buckling_load_N is None


def test_semi_analytical_cai_returns_less_than_pristine():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState(
        [DelaminationEllipse(3, (0, 0), 30, 20, 0)],
        dent_depth_mm=0.5,
    )
    result = semi_analytical_cai(lam, ds, 500.0, 15000.0)
    assert isinstance(result, SemiAnalyticalResult)
    assert 0 < result.residual_strength_MPa < 500.0
    assert result.critical_interface_index == 3
    assert result.critical_buckling_load_N is not None
    assert result.critical_buckling_load_N > 0


def test_semi_analytical_cai_positive_for_large_delamination():
    """A large delamination still leaves a positive residual (min of Soutis
    and delamination growth; see test_sublaminate_growth.py)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds_large = DamageState(
        [DelaminationEllipse(3, (0, 0), 60, 40, 0)],
        dent_depth_mm=0.5,
    )
    result_large = semi_analytical_cai(lam, ds_large, 500.0, 15000.0)
    assert result_large.residual_strength_MPa > 0


def test_semi_analytical_result_is_frozen_dataclass():
    """``SemiAnalyticalResult`` is an immutable, named-attribute container —
    callers must not be able to mutate its fields after construction, and
    the structured access must replace positional-tuple unpacking."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState(
        [DelaminationEllipse(3, (0, 0), 30, 20, 0)],
        dent_depth_mm=0.5,
    )
    result = semi_analytical_cai(lam, ds, 500.0, 15000.0)
    # Frozen dataclass: assignment to a field raises.
    with pytest.raises((AttributeError, Exception)):
        result.residual_strength_MPa = 0.0  # type: ignore[misc]
    # Named-attribute access (the whole point of this dataclass).
    assert hasattr(result, "residual_strength_MPa")
    assert hasattr(result, "critical_interface_index")
    assert hasattr(result, "critical_buckling_load_N")


def test_semi_analytical_tai_pristine_returns_pristine():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState([], dent_depth_mm=0.0)
    sigma_tai = semi_analytical_tai(lam, ds, 800.0)
    assert sigma_tai == 800.0


def test_semi_analytical_tai_monotonic_in_dpa():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds_small = DamageState([DelaminationEllipse(3, (0, 0), 10, 5, 0)], 0.3)
    ds_large = DamageState([DelaminationEllipse(3, (0, 0), 30, 20, 0)], 0.3)
    t1 = semi_analytical_tai(lam, ds_small, 800.0)
    t2 = semi_analytical_tai(lam, ds_large, 800.0)
    assert t1 > t2


def test_buckling_load_matches_closed_form_ssss_square():
    """Regression test for Issue #29.

    The closed-form SSSS Navier eigenvalue for an orthotropic rectangle is

        N_cr(m,n) = (pi^2 / a^2)
                    * [D11*m^4 + 2*(D12 + 2*D66)*(m*a/b)^2 * n^2
                       + D22 * (a*n/b)^4] / m^2

    where (a, b) are the FULL plate dimensions (the sine basis
    sin(m*pi*x/a) * sin(n*pi*y/b) vanishes at x=0, x=a, y=0, y=b).
    For an ellipse with semi-axes (major_mm, minor_mm), the enclosing
    rectangle is (2*major_mm) x (2*minor_mm), so a = 2*major_mm.

    Prior to the fix, the implementation passed a = major_mm directly,
    inflating N_cr by ~4x (N_cr ∝ 1/a^2). This test pins the corrected
    convention against a hand-computed value for a square sublaminate
    with a balanced [45, -45]s layup (D11 = D22 and m=n=1 governs).
    """
    m = MATERIAL_LIBRARY["IM7/8552"]
    # Balanced ±45 sublaminate -> D11 == D22 by symmetry; D16 == D26 != 0
    # but the closed form uses only D11, D22, D12, D66.
    layup = [45.0, -45.0, -45.0, 45.0]
    ply_t = 0.152
    lam = Laminate(m, layup, ply_t)

    # Square ellipse so the enclosing rectangle is square (a == b),
    # which makes m=n=1 the analytic minimizer for any D with D11 == D22.
    semi = 25.0  # mm
    ell = DelaminationEllipse(
        interface_index=1,
        centroid_mm=(0.0, 0.0),
        major_mm=semi,
        minor_mm=semi,
        orientation_deg=0.0,
    )

    # The sublaminate selection picks the thinner stack. With interface=1
    # and 4 plies total: upper = plies 0..1 (2 plies), lower = plies 2..3
    # (2 plies); ties go to upper. So sub_layup = [45, -45].
    sub_layup = [45.0, -45.0]
    D = _sublaminate_D_matrix(m, sub_layup, ply_t)
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]

    # Full rectangle side lengths (a, b) per the corrected convention.
    a = 2.0 * semi
    b = 2.0 * semi
    aspect = a / b  # == 1

    # Hand-evaluate the closed form at m=n=1.
    pi2 = math.pi * math.pi
    N_cr_11 = (
        (pi2 / a**2)
        * (
            D11 * 1**4
            + 2.0 * (D12 + 2.0 * D66) * (1 * aspect) ** 2 * 1**2
            + D22 * (aspect * 1) ** 4
        )
        / 1**2
    )

    # Sanity: m=n=1 really is the minimizer here (a == b, D11 == D22 by ±45
    # symmetry). Spot-check (2,1) and (1,2) to be safe.
    def _N_mn(mm, nn):
        return (
            (pi2 / a**2)
            * (
                D11 * mm**4
                + 2.0 * (D12 + 2.0 * D66) * (mm * aspect) ** 2 * nn**2
                + D22 * (aspect * nn) ** 4
            )
            / mm**2
        )

    assert _N_mn(2, 1) > N_cr_11
    assert _N_mn(1, 2) > N_cr_11
    assert _N_mn(2, 2) > N_cr_11

    N_actual = sublaminate_buckling_load(lam, ell, boundary="simply_supported")
    assert N_actual == pytest.approx(N_cr_11, rel=1e-10)

    # And confirm the boundary factor is applied.
    N_clamped = sublaminate_buckling_load(lam, ell, boundary="clamped")
    assert N_clamped == pytest.approx(N_cr_11 * 1.9, rel=1e-10)


def test_buckling_load_uses_full_plate_dimensions_not_semi_axes():
    """Regression test for Issue #29 — verify the 4x-overprediction bug
    does not regress.

    If the implementation reverted to using ellipse semi-axes directly as
    (a, b), N_cr would be inflated by exactly (2*r / r)^2 = 4 (since the
    eigenvalue is homogeneous of degree -2 in length when a/b is held fixed).
    Build two cases differing only in scale and check the corrected
    scaling, then independently check the magnitude against the
    "semi-axis-as-a" miscalculation.
    """
    mat = MATERIAL_LIBRARY["IM7/8552"]
    ply_t = 0.152
    lam = Laminate(mat, [0.0, 90.0, 90.0, 0.0], ply_t)

    semi = 20.0  # mm
    ell = DelaminationEllipse(1, (0.0, 0.0), semi, semi, 0.0)

    # Sub layup is the thinner of [0, 90] vs [90, 0]; tie -> upper -> [0, 90].
    sub_layup = [0.0, 90.0]
    D = _sublaminate_D_matrix(mat, sub_layup, ply_t)
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]

    # Closed form with the CORRECT full dimensions (a = b = 2*semi).
    a_full = 2.0 * semi
    pi2 = math.pi * math.pi
    N_correct = (pi2 / a_full**2) * (D11 + 2.0 * (D12 + 2.0 * D66) + D22)

    # Closed form with the BUGGY semi-axis-as-side-length (a = semi).
    a_bug = semi
    N_buggy = (pi2 / a_bug**2) * (D11 + 2.0 * (D12 + 2.0 * D66) + D22)

    # The bug would inflate by exactly 4x.
    assert N_buggy == pytest.approx(4.0 * N_correct, rel=1e-12)

    # Verify D11 == D22 for this symmetric cross-ply (so m=n=1 governs at a=b)
    assert D11 == pytest.approx(D22, rel=1e-12)

    N_actual = sublaminate_buckling_load(lam, ell, boundary="simply_supported")
    assert N_actual == pytest.approx(N_correct, rel=1e-10)
    # And explicitly rule out the buggy value.
    assert N_actual < 0.5 * N_buggy


# --- #165: the sublaminate rectangle is the rotated ellipse's bounding box ---

_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]


def _lam8():
    return Laminate(MATERIAL_LIBRARY["IM7/8552"], _QI8, 0.152)


def test_buckling_load_rotated_90_equals_swapped_axes():
    """An ellipse elongated across the load (orientation 90 deg) must buckle
    like the unrotated ellipse with its semi-axes swapped; before #165
    ``orientation_deg`` was ignored and both gave the along-load value."""
    lam = _lam8()
    across = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 90.0))
    swapped = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 8.0, 20.0, 0.0))
    along = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 0.0))
    assert across == pytest.approx(swapped, rel=1e-12)
    assert across != pytest.approx(along, rel=1e-3)
    # Fibre directions are axial: 180 deg is 0 deg and -90 deg is 90 deg.
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 180.0)
    ) == pytest.approx(along, rel=1e-12)
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, -90.0)
    ) == pytest.approx(across, rel=1e-12)


def test_buckling_load_45_deg_is_the_square_bounding_box():
    """At 45 deg the bounding box is square with half-side
    hypot(major, minor) / sqrt(2): larger in both directions than either
    axis-aligned box, so the enclosing-rectangle N_cr is below both the
    along-load and the across-load values (it does not lie between them)."""
    lam = _lam8()
    tilted = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 45.0))
    r = math.hypot(20.0, 8.0) / math.sqrt(2.0)
    square = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), r, r, 0.0))
    assert tilted == pytest.approx(square, rel=1e-12)
    along = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 0.0))
    across = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 90.0))
    assert tilted < min(along, across)


def test_buckling_load_pins_issue_165_probe():
    """The 20 x 8 mm ellipse at interface 2 of the 8-ply QI laminate from
    #165: 78.4 N/mm along the load (unchanged) and 39.5 N/mm across it
    (was 78.4, the orientation being ignored)."""
    lam = _lam8()
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 0.0)
    ) == pytest.approx(78.37, abs=0.05)
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 90.0)
    ) == pytest.approx(39.50, abs=0.05)


def test_buckling_rectangle_uses_the_tai_width_extents():
    """The sublaminate rectangle and the TAI damage width read the same
    extents off the ellipse: buckling a rotated ellipse equals buckling the
    axis-aligned ellipse of its half-extents, and twice its y half-extent is
    the TAI width."""
    from bvidfe.failure.soutis_openhole import transverse_damage_width_mm

    lam = _lam8()
    for theta in (0.0, 30.0, 45.0, 90.0, -45.0):
        e = DelaminationEllipse(3, (75, 50), 18.0, 7.0, theta)
        ex, ey = e.half_extents_mm()
        aligned = DelaminationEllipse(3, (75, 50), ex, ey, 0.0)
        assert sublaminate_buckling_load(lam, e) == pytest.approx(
            sublaminate_buckling_load(lam, aligned), rel=1e-12
        )
        assert transverse_damage_width_mm(DamageState([e])) == pytest.approx(2.0 * ey, rel=1e-12)


def test_template_damage_growth_sees_the_orientation():
    """End-to-end regression pins for #165 on 150 x 100 mm IM7/8552 panels
    (16 mm, 5.5 kg impactor). The peanut template elongates the back-face
    delaminations along the +/-45 and 90 deg plies, so they are wider
    across the load than ``2 * major_mm`` along it: on the 8-ply QI panel at
    5 J the growth stress falls from 300.5 to 281.8 MPa (knockdown 0.500 ->
    0.469, growth governs); at 10 J the film buckling stress at the
    governing interface 4 falls from 210 to 78 MPa but Soutis still governs
    (0.448). The [0/90]2s panel keeps 0.424 / 0.389 at 5 / 10 J."""
    from bvidfe import AnalysisConfig, BvidAnalysis
    from bvidfe.analysis.semi_analytical import weakest_sublaminate_growth
    from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
    from bvidfe.impact.mapping import ImpactEvent

    def run(layup, energy):
        cfg = AnalysisConfig(
            material="IM7/8552",
            layup_deg=layup,
            ply_thickness_mm=0.152,
            panel=PanelGeometry(150.0, 100.0),
            loading="compression",
            tier="semi_analytical",
            impact=ImpactEvent(energy, ImpactorGeometry(diameter_mm=16.0), mass_kg=5.5),
        )
        return BvidAnalysis(cfg).run()

    r5 = run(_QI8, 5.0)
    assert r5.knockdown == pytest.approx(0.469, abs=5e-4)
    assert r5.residual_strength_MPa == pytest.approx(281.8, abs=0.1)
    r10 = run(_QI8, 10.0)
    assert r10.knockdown == pytest.approx(0.448, abs=5e-4)
    g10 = weakest_sublaminate_growth(_lam8(), r10.damage)
    assert g10.interface_index == 4
    assert g10.film_buckling_stress_MPa == pytest.approx(78.4, abs=0.1)
    assert g10.growth_stress_MPa > r10.residual_strength_MPa  # Soutis governs
    xp = [0, 90, 0, 90, 90, 0, 90, 0]
    assert run(xp, 5.0).knockdown == pytest.approx(0.424, abs=5e-4)
    assert run(xp, 10.0).knockdown == pytest.approx(0.389, abs=5e-4)
