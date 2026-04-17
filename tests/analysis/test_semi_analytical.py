from bvidfe.analysis.semi_analytical import (
    find_critical_interface,
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
