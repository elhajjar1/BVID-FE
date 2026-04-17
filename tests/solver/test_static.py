"""Tests for the linear static FE solver."""

import numpy as np

from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.elements.hex8 import Hex8Element
from bvidfe.solver.boundary import compression_bcs
from bvidfe.solver.static import solve_linear_static


def _single_cube_model():
    """Create a simple single-element cube for testing."""
    m = MATERIAL_LIBRARY["IM7/8552"]
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
    elem = Hex8Element(node_coords, m)
    # global DOF map: node i -> DOFs 3*i, 3*i+1, 3*i+2
    dof_map = np.arange(24)
    return node_coords, [elem], [dof_map], m


def test_solve_compression_returns_prescribed_displacement_at_xmax():
    """Test that prescribed displacement is achieved at x_max."""
    node_coords, elems, dofs, m = _single_cube_model()
    n_dof = 3 * node_coords.shape[0]
    bcs = compression_bcs(node_coords, applied_strain=-0.01)
    u = solve_linear_static(elems, dofs, n_dof, bcs)
    assert u.shape == (n_dof,)
    # Node 1 is at x=1 (x_max); u_x should be ≈ -0.01
    assert abs(u[3 * 1 + 0] - (-0.01)) < 1e-3


def test_solve_preserves_rigid_y_displacement_at_xmin():
    """Test that clamped boundary conditions are preserved at x_min."""
    node_coords, elems, dofs, m = _single_cube_model()
    n_dof = 3 * node_coords.shape[0]
    bcs = compression_bcs(node_coords, applied_strain=-0.01)
    u = solve_linear_static(elems, dofs, n_dof, bcs)
    # Node 0 at x=y=z=0: u_x=u_y=u_z=0
    for k in range(3):
        assert abs(u[0 + k]) < 1e-6


def test_solve_axial_strain_matches_prescribed_for_single_cube():
    """Axial strain in a single cube under uniform compression should match the prescribed value."""
    node_coords, elems, dofs, m = _single_cube_model()
    n_dof = 3 * node_coords.shape[0]
    applied_strain = -0.01
    bcs = compression_bcs(node_coords, applied_strain=applied_strain)
    u = solve_linear_static(elems, dofs, n_dof, bcs)
    # Node 2 at x=1: u_x ≈ applied_strain * Lx
    assert abs(u[3 * 2 + 0] - applied_strain) < 1e-3
