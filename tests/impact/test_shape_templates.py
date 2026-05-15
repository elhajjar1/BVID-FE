from bvidfe.damage.state import DamageState
from bvidfe.impact.shape_templates import distribute_damage


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


def test_back_face_delaminations_grow_relative_to_impact_face():
    """Peanut-template back-face growth axiom, pinned directly.

    A constant-angle layup makes aspect_ratio == 1 at every interface, so
    all size variation comes purely from the _relative_size back-face ramp
    (0.3 near the impact face -> 1.0 near the back face). The existing
    count / DPA-conservation / aspect-ratio tests all still pass if that
    slope is flipped or zeroed (DPA conservation auto-rescales; aspect
    ratios are independent). This asserts the physics claim itself:
    delaminations near the back face are markedly larger than near the
    impact face.
    """
    ellipses = distribute_damage(
        layup_deg=[0] * 6,  # 6 plies => 5 interfaces, AR == 1 everywhere
        target_dpa_mm2=600.0,
        dent_depth_mm=0.3,
        fiber_break_radius_mm=0.0,
        centroid_mm=(75.0, 50.0),
    )
    ellipses.sort(key=lambda e: e.interface_index)
    # Model ratio is rel(back)/rel(impact) = 1.0 / 0.44 ~= 2.3; assert >= 2x.
    assert ellipses[-1].major_mm > 2.0 * ellipses[0].major_mm
    assert ellipses[-1].minor_mm > 2.0 * ellipses[0].minor_mm


def test_empty_when_single_ply():
    assert distribute_damage([0], 100.0, 0.3, 0.0) == []


def test_empty_when_nonpositive_target():
    assert distribute_damage([0, 90, 0], 0.0, 0.3, 0.0) == []
    assert distribute_damage([0, 90, 0], -5.0, 0.3, 0.0) == []
