import numpy as np
import scipy.sparse as sp

from bvidfe.solver.boundary import (
    BoundaryCondition,
    apply_dirichlet_penalty,
    compression_bcs,
    tension_bcs,
)


def test_apply_dirichlet_penalty_enforces_prescribed_value():
    # Small 3-DOF system: u = [u0, u1, u2]
    K = sp.csc_matrix(np.array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]], dtype=float))
    F = np.zeros(3)
    bcs = [BoundaryCondition(dof=0, value=1.5), BoundaryCondition(dof=2, value=-0.5)]
    Kmod, Fmod = apply_dirichlet_penalty(K, F, bcs, penalty=1e10)
    u = sp.linalg.spsolve(Kmod, Fmod)
    assert abs(u[0] - 1.5) < 1e-6
    assert abs(u[2] - (-0.5)) < 1e-6


def test_compression_bcs_returns_xmin_and_xmax_bcs():
    # 8 nodes of a unit cube
    node_coords = np.array(
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
    bcs = compression_bcs(node_coords, applied_strain=-0.01)
    # 4 nodes on x_min get u_x=0, 4 nodes on x_max get u_x=-0.01*Lx=-0.01
    xmin_bcs = [b for b in bcs if abs(b.value) < 1e-12]
    xmax_bcs = [b for b in bcs if abs(b.value - (-0.01)) < 1e-6]
    assert len(xmin_bcs) >= 4  # includes symmetry y and z
    assert len(xmax_bcs) == 4


def test_tension_bcs_positive_displacement():
    node_coords = np.array(
        [
            [0, 0, 0],
            [2, 0, 0],
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [2, 0, 1],
            [2, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    bcs = tension_bcs(node_coords, applied_strain=0.005)
    xmax = [b for b in bcs if abs(b.value - (0.005 * 2)) < 1e-6]
    assert len(xmax) == 4
