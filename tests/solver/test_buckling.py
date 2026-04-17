import numpy as np
import scipy.sparse as sp

from bvidfe.solver.buckling import linear_buckling


def test_linear_buckling_2x2_diagonal_eigenvalues():
    """Handcrafted 2x2 diagonal generalized eigenproblem with known eigenvalues."""
    K = sp.csc_matrix(np.diag([2.0, 4.0]))
    Kg = sp.csc_matrix(np.diag([1.0, 1.0]))
    eigs, modes = linear_buckling(K, Kg, n_modes=2)
    eigs = sorted(eigs)
    assert abs(eigs[0] - 2.0) < 1e-6
    assert abs(eigs[1] - 4.0) < 1e-6


def test_linear_buckling_returns_requested_number_of_modes():
    # Small 5x5 SPD system
    rng = np.random.default_rng(42)
    A = rng.random((5, 5))
    K_dense = A @ A.T + 5 * np.eye(5)
    Kg_dense = np.eye(5)
    K = sp.csc_matrix(K_dense)
    Kg = sp.csc_matrix(Kg_dense)
    eigs, modes = linear_buckling(K, Kg, n_modes=3)
    assert len(eigs) == 3
    assert modes.shape == (5, 3)


def test_linear_buckling_modes_are_normalized_or_finite():
    """Mode shapes are finite and non-trivial."""
    K = sp.csc_matrix(np.diag([2.0, 4.0, 6.0]))
    Kg = sp.csc_matrix(np.eye(3))
    eigs, modes = linear_buckling(K, Kg, n_modes=2)
    assert np.all(np.isfinite(modes))
    # Each column mode should have at least one nonzero entry
    col_norms = np.linalg.norm(modes, axis=0)
    assert np.all(col_norms > 1e-6)
