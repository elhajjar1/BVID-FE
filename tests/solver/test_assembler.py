import numpy as np
import scipy.sparse as sp

from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.elements.hex8 import Hex8Element
from bvidfe.solver.assembler import assemble_global_stiffness


def _single_element_system():
    m = MATERIAL_LIBRARY["IM7/8552"]
    nodes = np.array(
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
    elem = Hex8Element(nodes, m)
    dof_map = np.arange(24)  # 24 DOFs for 8 nodes
    return [elem], [dof_map], 24


def _two_element_system():
    """Two unit hex elements sharing the face x=1."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    # 12 nodes total; first cube occupies x in [0,1], second in [1,2]
    n_nodes = 12
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
            [2, 0, 0],
            [2, 1, 0],
            [2, 0, 1],
            [2, 1, 1],
        ],
        dtype=float,
    )
    # Element 1 nodes indices
    e1_nodes = [0, 1, 2, 3, 4, 5, 6, 7]
    # Element 2 nodes indices (shares nodes 1,2,5,6)
    e2_nodes = [1, 8, 9, 2, 5, 10, 11, 6]
    elem1 = Hex8Element(node_coords[e1_nodes], m)
    elem2 = Hex8Element(node_coords[e2_nodes], m)

    # Global DOF index arrays (24 each)
    def _dof(ni):
        return np.array([3 * n + k for n in ni for k in range(3)])

    return [elem1, elem2], [_dof(e1_nodes), _dof(e2_nodes)], 3 * n_nodes


def test_assemble_single_element_shape():
    elems, dofs, n_dof = _single_element_system()
    K = assemble_global_stiffness(elems, dofs, n_dof)
    assert sp.issparse(K)
    assert K.shape == (24, 24)
    assert K.nnz > 0


def test_assemble_symmetric():
    elems, dofs, n_dof = _two_element_system()
    K = assemble_global_stiffness(elems, dofs, n_dof)
    diff = K - K.T
    assert np.allclose(diff.data, 0.0, atol=1e-6)


def test_assemble_six_rigid_body_modes():
    elems, dofs, n_dof = _two_element_system()
    K = assemble_global_stiffness(elems, dofs, n_dof)
    # Convert to dense for eigenvalue analysis (small system)
    Kd = K.toarray()
    eigs = np.linalg.eigvalsh(Kd)
    # 6 rigid-body zero eigenvalues expected (within numerical tolerance)
    near_zero = eigs[np.abs(eigs) < 1e-3]
    assert len(near_zero) == 6


def test_assemble_summation_of_overlapping_dofs():
    """Shared face DOFs should sum contributions from both elements."""
    elems, dofs, n_dof = _two_element_system()
    K = assemble_global_stiffness(elems, dofs, n_dof)
    # Node 1 is shared; K[3,3] should include contributions from both elements
    # We can't easily check the exact value, but K should be positive-definite on the constrained system
    Kd = K.toarray()
    free = list(range(n_dof))
    # Fix 6 DOFs to remove rigid-body modes
    fixed = [0, 1, 2, 4, 5, 11]  # arbitrary but enough
    for f in sorted(fixed, reverse=True):
        free.remove(f)
    Krr = Kd[np.ix_(free, free)]
    eigs = np.linalg.eigvalsh(Krr)
    assert eigs.min() > 0


def test_assembly_sums_overlapping_dof_contributions():
    """Shared-DOF summation invariant, pinned exactly.

    Where two elements share a node, that node's global stiffness block
    must equal the *sum* of each element's local block — not an overwrite
    or a deduplicated (row, col) pair. A regression flipping the
    assembler's accumulate to an assignment would silently corrupt every
    multi-element solve; the symmetry / PD / rigid-body tests above would
    not catch it.
    """
    elems, dofs, n_dof = _two_element_system()
    Ke1 = elems[0].stiffness_matrix()
    Ke2 = elems[1].stiffness_matrix()
    K = assemble_global_stiffness(elems, dofs, n_dof).toarray()

    # _two_element_system: e1 nodes [0,1,2,3,4,5,6,7], e2 nodes [1,8,9,2,5,10,11,6].
    # Node 1 is shared: local index 1 in element 1, local index 0 in element 2.
    e1_local, e2_local = 1, 0
    g = [3 * 1 + k for k in range(3)]  # node 1's global DOFs
    l1 = [3 * e1_local + k for k in range(3)]
    l2 = [3 * e2_local + k for k in range(3)]
    expected = Ke1[np.ix_(l1, l1)] + Ke2[np.ix_(l2, l2)]
    np.testing.assert_allclose(K[np.ix_(g, g)], expected, atol=1e-9, rtol=0)

    # Sanity: an unshared node (node 0, only in element 1) must equal
    # exactly element 1's local block — no spurious doubling.
    g0 = [0, 1, 2]
    np.testing.assert_allclose(
        K[np.ix_(g0, g0)], Ke1[np.ix_(g0, g0)], atol=1e-9, rtol=0
    )
