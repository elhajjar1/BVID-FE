"""Sparse global stiffness assembly from element contributions."""

from __future__ import annotations

from typing import Sequence

import numpy as np
import scipy.sparse as sp


def assemble_global_stiffness(
    elements: Sequence,
    element_dof_maps: Sequence[np.ndarray],
    n_dof: int,
) -> sp.csc_matrix:
    """Assemble the global stiffness matrix in CSC format.

    Parameters
    ----------
    elements : sequence of elements each exposing `.stiffness_matrix() -> (24, 24)`
    element_dof_maps : sequence of int arrays, each length 24, giving the global DOF
        index for each of the 24 element DOFs.
    n_dof : total number of global DOFs.
    """
    if len(elements) != len(element_dof_maps):
        raise ValueError(
            f"elements (n={len(elements)}) and element_dof_maps (n={len(element_dof_maps)}) length mismatch"
        )

    # Vectorized COO assembly: one numpy meshgrid per element instead of
    # a 24x24 Python loop. ~10x faster on large meshes (the inner Python
    # overhead was the dominant cost on 4k+ element meshes).
    n_elem = len(elements)
    rows_chunks: list[np.ndarray] = []
    cols_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []

    for elem, dof_map in zip(elements, element_dof_maps):
        Ke = elem.stiffness_matrix()
        if Ke.shape != (24, 24):
            raise ValueError(f"element stiffness must be (24,24), got {Ke.shape}")
        dof_arr = np.asarray(dof_map, dtype=np.int64)
        if dof_arr.shape != (24,):
            raise ValueError(f"dof_map must have 24 entries, got {len(dof_map)}")
        # Outer product of dof indices → 24x24 row/col index arrays
        rows_chunks.append(np.broadcast_to(dof_arr[:, None], (24, 24)).ravel())
        cols_chunks.append(np.broadcast_to(dof_arr[None, :], (24, 24)).ravel())
        data_chunks.append(Ke.ravel())

    if n_elem == 0:
        return sp.csc_matrix((n_dof, n_dof))

    rows = np.concatenate(rows_chunks)
    cols = np.concatenate(cols_chunks)
    data = np.concatenate(data_chunks)

    K = sp.coo_matrix((data, (rows, cols)), shape=(n_dof, n_dof))
    return K.tocsc()
