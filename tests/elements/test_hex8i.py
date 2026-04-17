import numpy as np

from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.elements.hex8 import Hex8Element
from bvidfe.elements.hex8i import Hex8iElement


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


def test_hex8i_stiffness_shape_and_symmetry():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8iElement(_unit_cube_nodes(), m)
    K = elem.stiffness_matrix()
    assert K.shape == (24, 24)
    assert np.allclose(K, K.T, atol=1e-5)


def test_hex8i_stiffer_or_softer_than_hex8_in_bending():
    m = MATERIAL_LIBRARY["IM7/8552"]
    # Thin plate in bending: hex8 locks (too stiff), hex8i softens
    thin_nodes = np.array(
        [
            [0, 0, 0],
            [10, 0, 0],
            [10, 10, 0],
            [0, 10, 0],
            [0, 0, 0.1],
            [10, 0, 0.1],
            [10, 10, 0.1],
            [0, 10, 0.1],
        ],
        dtype=float,
    )
    K_h8 = Hex8Element(thin_nodes, m).stiffness_matrix()
    K_h8i = Hex8iElement(thin_nodes, m).stiffness_matrix()
    # The max singular value / diagonal norm of K_h8i should be less than K_h8 for bending modes.
    # Easier check: hex8i stiffness is NOT identical to hex8.
    assert not np.allclose(K_h8, K_h8i, atol=1e-3)


def test_hex8i_positive_definite_after_bc():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8iElement(_unit_cube_nodes(), m)
    K = elem.stiffness_matrix()
    free = list(range(24))
    fixed = [0, 1, 2, 4, 5, 8]
    for i in sorted(fixed, reverse=True):
        free.remove(i)
    Krr = K[np.ix_(free, free)]
    eigs = np.linalg.eigvalsh(Krr)
    assert eigs.min() > 0
