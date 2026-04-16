from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.tsai_wu import tsai_wu_index, tsai_wu_strength_uniaxial


def test_uniaxial_Xt_gives_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    s = [m.Xt, 0, 0, 0, 0, 0]
    assert abs(tsai_wu_index(m, s) - 1.0) < 0.05


def test_uniaxial_negative_Xc_gives_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    s = [-m.Xc, 0, 0, 0, 0, 0]
    assert abs(tsai_wu_index(m, s) - 1.0) < 0.05


def test_zero_stress_gives_zero_index():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert tsai_wu_index(m, [0] * 6) == 0.0


def test_strength_uniaxial_tension_matches_Xt_within_10pct():
    m = MATERIAL_LIBRARY["IM7/8552"]
    sig = tsai_wu_strength_uniaxial(m, direction=1, sign=+1)
    assert 0.9 * m.Xt <= sig <= 1.1 * m.Xt


def test_strength_uniaxial_compression_matches_Xc_within_10pct():
    m = MATERIAL_LIBRARY["IM7/8552"]
    sig = tsai_wu_strength_uniaxial(m, direction=1, sign=-1)
    assert 0.9 * m.Xc <= sig <= 1.1 * m.Xc
