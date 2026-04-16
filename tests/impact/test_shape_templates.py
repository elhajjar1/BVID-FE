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


def test_empty_when_single_ply():
    assert distribute_damage([0], 100.0, 0.3, 0.0) == []


def test_empty_when_nonpositive_target():
    assert distribute_damage([0, 90, 0], 0.0, 0.3, 0.0) == []
    assert distribute_damage([0, 90, 0], -5.0, 0.3, 0.0) == []
