import numpy as np

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


def test_shape_functions_partition_of_unity():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    for xi, eta, zeta in [(-0.5, 0.3, 0.1), (0.7, -0.2, 0.8), (0, 0, 0)]:
        N = elem.shape_functions(xi, eta, zeta)
        assert N.shape == (8,)
        assert abs(N.sum() - 1.0) < 1e-12


def test_shape_function_nodes_values():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    # At node 0 corner (-1,-1,-1), N_0 = 1, others = 0
    N = elem.shape_functions(-1, -1, -1)
    assert abs(N[0] - 1.0) < 1e-12
    assert np.allclose(N[1:], 0.0, atol=1e-12)


def test_jacobian_unit_cube_is_half():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    J = elem.jacobian(0, 0, 0)
    assert J.shape == (3, 3)
    # Unit cube (side = 1 mm), natural length = 2, so dx/dxi = 0.5
    assert abs(np.linalg.det(J) - 0.125) < 1e-12  # det(J) = (0.5)^3


def test_B_matrix_shape():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    B, detJ = elem.B_matrix(0, 0, 0)
    assert B.shape == (6, 24)
    assert detJ > 0


def test_stiffness_matrix_shape_and_symmetry():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    K = elem.stiffness_matrix()
    assert K.shape == (24, 24)
    assert np.allclose(K, K.T, atol=1e-6)


def test_stiffness_matrix_positive_definite_after_bc():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m)
    K = elem.stiffness_matrix()
    # Kill 6 DOF to remove rigid body modes (fix node 0 fully + node 1 yz + node 3 z)
    free = list(range(24))
    fixed = [0, 1, 2, 4, 5, 8]
    for i in sorted(fixed, reverse=True):
        free.remove(i)
    Krr = K[np.ix_(free, free)]
    eigs = np.linalg.eigvalsh(Krr)
    assert eigs.min() > 0


def test_ply_rotation_produces_stiffness_in_laminate_frame():
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem0 = Hex8Element(_unit_cube_nodes(), m, ply_angle_deg=0.0)
    elem90 = Hex8Element(_unit_cube_nodes(), m, ply_angle_deg=90.0)
    # K differs between 0 and 90 because the fiber direction rotates
    K0 = elem0.stiffness_matrix()
    K90 = elem90.stiffness_matrix()
    assert not np.allclose(K0, K90, atol=1.0)
