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

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    for elem, dof_map in zip(elements, element_dof_maps):
        Ke = elem.stiffness_matrix()
        if Ke.shape != (24, 24):
            raise ValueError(f"element stiffness must be (24,24), got {Ke.shape}")
        if len(dof_map) != 24:
            raise ValueError(f"dof_map must have 24 entries, got {len(dof_map)}")
        for i in range(24):
            gi = int(dof_map[i])
            for j in range(24):
                gj = int(dof_map[j])
                rows.append(gi)
                cols.append(gj)
                data.append(Ke[i, j])

    K = sp.coo_matrix((data, (rows, cols)), shape=(n_dof, n_dof))
    return K.tocsc()
