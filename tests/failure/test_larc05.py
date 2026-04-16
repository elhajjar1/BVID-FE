from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.larc05 import larc05_index


def test_uniaxial_Xt_tension_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert abs(larc05_index(m, [m.Xt, 0, 0, 0, 0, 0]) - 1.0) < 1e-6


def test_uniaxial_Xc_compression_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert abs(larc05_index(m, [-m.Xc, 0, 0, 0, 0, 0]) - 1.0) < 1e-6


def test_uniaxial_Yt_transverse_tension_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert abs(larc05_index(m, [0, m.Yt, 0, 0, 0, 0]) - 1.0) < 1e-6


def test_uniaxial_Yc_transverse_compression_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert abs(larc05_index(m, [0, -m.Yc, 0, 0, 0, 0]) - 1.0) < 1e-6


def test_pure_in_plane_shear_S12():
    m = MATERIAL_LIBRARY["IM7/8552"]
    # S12 shear alone should give sqrt(1.0) = 1 in matrix modes
    # our modes: tau_12 contributes to both mt and mc equally
    idx = larc05_index(m, [0, 0, 0, 0, 0, m.S12])
    assert abs(idx - 1.0) < 1e-6


def test_zero_stress_gives_zero():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert larc05_index(m, [0, 0, 0, 0, 0, 0]) == 0.0


def test_returns_max_over_modes():
    m = MATERIAL_LIBRARY["IM7/8552"]
    # A stress state combining fiber tension (dominant) and some shear
    idx = larc05_index(m, [m.Xt * 0.9, 0, 0, 0, 0, m.S12 * 0.1])
    # Fiber tension mode f_ft = 0.81; shear-only in matrix mode = 0.01. Max = 0.81.
    assert abs(idx - 0.81) < 0.01
