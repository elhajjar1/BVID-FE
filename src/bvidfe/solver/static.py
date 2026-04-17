"""Linear static finite-element solver (direct sparse solve)."""

from __future__ import annotations

from typing import Sequence

import numpy as np
import scipy.sparse.linalg as spla

from bvidfe.solver.assembler import assemble_global_stiffness
from bvidfe.solver.boundary import BoundaryCondition, apply_dirichlet_penalty


def solve_linear_static(
    elements: Sequence,
    element_dof_maps: Sequence[np.ndarray],
    n_dof: int,
    bcs: Sequence[BoundaryCondition],
    external_force: np.ndarray | None = None,
    penalty: float = 1.0e10,
) -> np.ndarray:
    """Assemble, apply BCs via penalty, and solve K u = F.

    Parameters
    ----------
    elements : sequence of FE elements exposing `.stiffness_matrix() -> (24,24)`.
    element_dof_maps : global DOF indices per element.
    n_dof : total number of global DOFs.
    bcs : Dirichlet boundary conditions.
    external_force : optional applied nodal force vector (length n_dof, default zero).
    penalty : penalty scalar for Dirichlet enforcement.

    Returns
    -------
    np.ndarray
        Displacement vector of length n_dof.
    """
    K = assemble_global_stiffness(elements, element_dof_maps, n_dof)
    F = (
        np.zeros(n_dof)
        if external_force is None
        else np.asarray(external_force, dtype=float).copy()
    )
    K_mod, F_mod = apply_dirichlet_penalty(K, F, bcs, penalty=penalty)
    u = spla.spsolve(K_mod.tocsc(), F_mod)
    return np.asarray(u)
