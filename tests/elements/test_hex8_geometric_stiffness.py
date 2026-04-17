"""Tests for Hex8Element.geometric_stiffness_matrix."""

import numpy as np
import pytest

from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.elements.hex8 import Hex8Element


def _unit_cube_nodes():
    return np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )


def test_geometric_stiffness_shape_and_symmetry():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    sigma = np.array([[-1.0, 0, 0], [0, 0, 0], [0, 0, 0]])
    Kg = elem.geometric_stiffness_matrix(sigma)
    assert Kg.shape == (24, 24)
    assert np.allclose(Kg, Kg.T, atol=1e-6)


def test_geometric_stiffness_zero_for_zero_stress():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    Kg = elem.geometric_stiffness_matrix(np.zeros((3, 3)))
    assert np.allclose(Kg, 0.0, atol=1e-12)


def test_geometric_stiffness_linear_in_stress():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    sigma_a = np.array([[-1.0, 0, 0], [0, 0, 0], [0, 0, 0]])
    sigma_b = 2.5 * sigma_a
    Kg_a = elem.geometric_stiffness_matrix(sigma_a)
    Kg_b = elem.geometric_stiffness_matrix(sigma_b)
    assert np.allclose(Kg_b, 2.5 * Kg_a, atol=1e-6)


def test_geometric_stiffness_wrong_shape_raises():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    with pytest.raises(ValueError, match="sigma_bar must be"):
        elem.geometric_stiffness_matrix(np.zeros((6,)))
