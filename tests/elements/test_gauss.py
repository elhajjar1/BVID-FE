import pytest

from bvidfe.elements.gauss import gauss_points_1d, gauss_points_hex


def test_gauss_1d_n1():
    pts, wts = gauss_points_1d(1)
    assert pts.shape == (1,)
    assert wts.shape == (1,)
    assert abs(pts[0]) < 1e-12
    assert abs(wts[0] - 2.0) < 1e-12


def test_gauss_1d_n2_integrates_x2():
    pts, wts = gauss_points_1d(2)
    # integral of x^2 from -1 to 1 = 2/3
    val = sum(w * p * p for p, w in zip(pts, wts))
    assert abs(val - 2.0 / 3.0) < 1e-12


def test_gauss_1d_n3_integrates_x4():
    pts, wts = gauss_points_1d(3)
    # integral of x^4 from -1 to 1 = 2/5
    val = sum(w * p**4 for p, w in zip(pts, wts))
    assert abs(val - 2.0 / 5.0) < 1e-12


def test_gauss_1d_invalid_order_raises():
    with pytest.raises(ValueError):
        gauss_points_1d(0)
    with pytest.raises(ValueError):
        gauss_points_1d(4)


def test_gauss_hex_order2_has_8_points():
    gp, wt = gauss_points_hex(order=2)
    assert gp.shape == (8, 3)
    assert wt.shape == (8,)
    # Unit cube volume integral: sum(weights) = 2^3 = 8
    assert abs(wt.sum() - 8.0) < 1e-12


def test_gauss_hex_order3_has_27_points():
    gp, wt = gauss_points_hex(order=3)
    assert gp.shape == (27, 3)
    assert abs(wt.sum() - 8.0) < 1e-12
