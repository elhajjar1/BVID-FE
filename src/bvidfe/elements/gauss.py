"""Gauss-Legendre quadrature points and weights for 1D and 3D tensor-product
integration over the reference cube [-1, 1]^3."""

from __future__ import annotations

import numpy as np

# Tabulated 1D Gauss-Legendre points/weights for order n = 1, 2, 3.
_GAUSS_1D = {
    1: (np.array([0.0]), np.array([2.0])),
    2: (
        np.array([-1.0 / np.sqrt(3.0), +1.0 / np.sqrt(3.0)]),
        np.array([1.0, 1.0]),
    ),
    3: (
        np.array([-np.sqrt(3.0 / 5.0), 0.0, +np.sqrt(3.0 / 5.0)]),
        np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]),
    ),
}


def gauss_points_1d(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (points, weights) for n-point Gauss-Legendre on [-1, 1]. Supports n in {1, 2, 3}."""
    if n not in _GAUSS_1D:
        raise ValueError(f"gauss_points_1d supports n in {{1,2,3}} (got {n})")
    pts, wts = _GAUSS_1D[n]
    return pts.copy(), wts.copy()


def gauss_points_hex(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (xi_eta_zeta, weights) for the tensor product of 1D Gauss points over the 3D cube.

    Returns
    -------
    points : (n_gp, 3) array of reference coordinates
    weights : (n_gp,) array of weights (product of 1D weights)
    """
    p1, w1 = gauss_points_1d(order)
    pts = np.array([[xi, eta, zeta] for xi in p1 for eta in p1 for zeta in p1])
    wts = np.array([wi * wj * wk for wi in w1 for wj in w1 for wk in w1])
    return pts, wts
