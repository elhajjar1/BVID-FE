import pytest

from bvidfe.damage.state import DamageState
from bvidfe.impact.shape_templates import (
    _aspect_ratio,
    _orientation_deg,
    _ply_mismatch_deg,
    distribute_damage,
)


def test_one_ellipse_per_interface():
    layup = [0, 45, -45, 90, 0, 90, -45, 45, 0]  # 9 plies => 8 interfaces
    ellipses = distribute_damage(
        layup_deg=layup,
        target_dpa_mm2=1000.0,
        dent_depth_mm=0.4,
        fiber_break_radius_mm=0.0,
    )
    assert len({e.interface_index for e in ellipses}) == 8


def test_dpa_conservation_within_tolerance():
    layup = [0, 45, -45, 90, 90, -45, 45, 0]  # 8 plies => 7 interfaces
    target = 800.0
    ellipses = distribute_damage(
        layup_deg=layup,
        target_dpa_mm2=target,
        dent_depth_mm=0.3,
        fiber_break_radius_mm=0.0,
    )
    ds = DamageState(ellipses, 0.3, 0.0)
    assert abs(ds.projected_damage_area_mm2 - target) / target < 0.01


def test_aspect_ratio_grows_with_ply_angle_mismatch():
    ellipses_aligned = distribute_damage([0, 0, 0], 400.0, 0.3, 0.0)
    ellipses_cross = distribute_damage([0, 90, 0], 400.0, 0.3, 0.0)

    def _ar(es):
        return max(e.major_mm / e.minor_mm for e in es)

    assert _ar(ellipses_cross) > _ar(ellipses_aligned)


def test_empty_when_single_ply():
    assert distribute_damage([0], 100.0, 0.3, 0.0) == []


def test_empty_when_nonpositive_target():
    assert distribute_damage([0, 90, 0], 0.0, 0.3, 0.0) == []
    assert distribute_damage([0, 90, 0], -5.0, 0.3, 0.0) == []


def test_back_face_growth_axiom():
    """Issue #16: delaminations must grow toward the back face (the
    experimentally observed BVID signature encoded by _relative_size).

    Constant-angle layup -> aspect_ratio == 1 at every interface, so any
    size difference is purely the back-face growth ramp. A regression that
    flipped or zeroed the ramp slope would still pass the existing
    DPA-conservation / aspect-ratio tests but fail this one.
    """
    ellipses = distribute_damage(
        [0] * 6,
        target_dpa_mm2=600.0,
        dent_depth_mm=0.3,
        fiber_break_radius_mm=0.0,
        centroid_mm=(75.0, 50.0),
    )
    ellipses.sort(key=lambda e: e.interface_index)
    # Model ratio is (0.3 + 0.7*1) / (0.3 + 0.7*(1/n)) ~= 2.3; assert >= 2x.
    assert ellipses[-1].major_mm > 2.0 * ellipses[0].major_mm
    assert ellipses[-1].minor_mm > 2.0 * ellipses[0].minor_mm


# --- Fibre directions are axial: theta and theta + 180 are the same ply ----


def _footprint(layup, area=1500.0):
    ells = distribute_damage(layup, area, 0.3, 0.0, (0.0, 0.0))
    return [
        (e.interface_index, e.major_mm, e.minor_mm, e.orientation_deg)
        for e in sorted(ells, key=lambda e: e.interface_index)
    ]


_QI16 = [45, 0, -45, 90] * 2
_QI16 = _QI16 + _QI16[::-1]


@pytest.mark.parametrize(
    "rewrite",
    [
        {-45: 135},  # -45 written as 135
        {90: -90},  # 90 written as -90
        {45: -135, 0: 180},
    ],
)
def test_equivalent_layup_notations_give_identical_footprints(rewrite):
    other = [rewrite.get(a, a) for a in _QI16]
    assert other != _QI16
    a, b = _footprint(_QI16), _footprint(other)
    assert len(a) == len(b)
    for ea, eb in zip(a, b):
        assert ea[0] == eb[0] and ea[3] == eb[3]
        assert ea[1:3] == pytest.approx(eb[1:3], rel=1e-12)


@pytest.mark.parametrize(
    "a, b, expected",
    [(-45, 90, 45), (135, 90, 45), (45, -45, 90), (0, 90, 90), (90, -90, 0), (30, -60, 90)],
)
def test_ply_mismatch_is_the_acute_angle(a, b, expected):
    assert float(_ply_mismatch_deg(a, b)) == pytest.approx(expected)
    assert _aspect_ratio(b - a) == pytest.approx(1.0 + 0.025 * expected)


def test_minus45_90_interface_aspect_ratio_is_not_clipped():
    """-45/90 is a 45 deg mismatch (AR 2.125), not 135 deg (AR clipped to 4)."""
    (e,) = distribute_damage([-45, 90], 500.0, 0.3, 0.0)
    assert e.major_mm / e.minor_mm == pytest.approx(2.125)


@pytest.mark.parametrize(
    "impact_side, back_side, expected",
    [(45, 0, 0.0), (0, 45, 45.0), (-45, 90, 90.0), (90, 135, -45.0), (0, -90, 90.0)],
)
def test_orientation_is_the_back_face_side_ply(impact_side, back_side, expected):
    """Ply 0 is the impact face, so interface i elongates along ply i + 1."""
    assert _orientation_deg(impact_side, back_side) == expected
    (e,) = distribute_damage([impact_side, back_side], 500.0, 0.3, 0.0)
    assert e.orientation_deg == expected
    assert -90.0 < e.orientation_deg <= 90.0


def test_same_direction_interfaces_carry_no_delamination():
    layup = [0, 45, -45, 90, 90, -45, 45, 0]  # 90/90 at interface 3
    ells = distribute_damage(layup, 800.0, 0.3, 0.0)
    assert sorted(e.interface_index for e in ells) == [0, 1, 2, 4, 5, 6]
    # 90 and -90 are the same direction too.
    ells = distribute_damage([0, 90, -90, 0], 400.0, 0.3, 0.0)
    assert sorted(e.interface_index for e in ells) == [0, 2]
    assert DamageState(ells).projected_damage_area_mm2 == pytest.approx(400.0, rel=0.01)


def test_unidirectional_layup_keeps_every_interface():
    ells = distribute_damage([0, 180, 0, 0], 300.0, 0.3, 0.0)
    assert sorted(e.interface_index for e in ells) == [0, 1, 2]
    assert all(e.major_mm == pytest.approx(e.minor_mm) for e in ells)


def test_layup_notation_does_not_change_the_analysis():
    """End to end: writing -45 as 135 leaves the semi-analytical CAI unchanged."""
    from bvidfe.analysis import AnalysisConfig, BvidAnalysis
    from bvidfe.core.geometry import PanelGeometry
    from bvidfe.impact.mapping import ImpactEvent

    def run(layup):
        cfg = AnalysisConfig(
            material="IM7/8552",
            layup_deg=layup,
            ply_thickness_mm=0.152,
            panel=PanelGeometry(150.0, 100.0),
            loading="compression",
            tier="semi_analytical",
            impact=ImpactEvent(energy_J=20.0),
        )
        return BvidAnalysis(cfg).run()

    a = run(_QI16)
    b = run([135 if x == -45 else x for x in _QI16])
    assert b.dpa_mm2 == pytest.approx(a.dpa_mm2, rel=1e-12)
    assert b.knockdown == pytest.approx(a.knockdown, rel=1e-12)
