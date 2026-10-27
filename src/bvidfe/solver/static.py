"""Linear static finite-element solver.

Dirichlet BCs are imposed by eliminating the constrained DOFs: with ``c``
the constrained and ``f`` the free DOFs, the solve is

    K_ff u_f = F_f - K_fc u_c

and ``u_c`` is copied into the result exactly. ``K_ff`` is symmetric
positive definite, so it is solved by conjugate gradients preconditioned
with an incomplete LU factorisation (``spilu``). When CG does not reach the
tolerance, GMRES retries with the same factor, and a direct sparse LU of
``K_ff`` is the last resort, so every returned solution meets the
tolerance.

Measured on an 8-ply 150 x 100 mm IM7/8552 panel at a 2.5 mm mesh
(67.5k DOF), the previous path (1e10 penalty on the diagonal, then
``spsolve`` with SuperLU's default COLAMD ordering) took ~110 s and ~2.9 GB
per solve. ILU-CG on the reduced system takes ~10 s in about 50
iterations, and agrees with the penalty solve to ~1e-7, which is the
penalty method's own constraint error.
"""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from bvidfe.solver.assembler import assemble_global_stiffness
from bvidfe.solver.boundary import BoundaryCondition, apply_dirichlet_penalty

_log = logging.getLogger("bvidfe.fe3d.solver")

# Relative residual ||rhs - K_ff u_f|| / ||rhs|| the iterative solve must
# reach. Tight enough that the fe3d knockdown does not move at the precision
# it is reported to (the reduced direct solve reaches ~1e-13).
DEFAULT_RTOL = 1.0e-10

# ILU preconditioner settings, from the 2.5 mm reduced system above
# (63k free DOFs, 4.5M nonzeros):
#   drop_tol  fill_factor  nnz(L+U)  factor  CG iterations
#   1e-3      10           11.9M      8.5 s  43
#   1e-4      10           21.7M     14.0 s  13
#   1e-5      10           33.1M     20.2 s  no convergence in 300 (GMRES: 230)
#   1e-5      30           34.9M     23.0 s  7
# SuperLU caps the fill column by column at fill_factor x the column of A
# and drops harder where the cap binds, which can leave an unstable factor
# CG cannot use: it is why drop_tol=1e-5 failed at fill_factor=10, and
# likely why drop_tol=1e-4 (fill_factor=10) failed at a 1.25 mm mesh (265k
# DOFs; not rerun). drop_tol=1e-3 is the cheapest overall and gives the same
# factor at fill_factor=5; the generous fill_factor keeps the cap from
# binding on finer meshes.
ILU_DROP_TOL = 1.0e-3
ILU_FILL_FACTOR = 30.0

# Iteration caps. A freshly factorised ILU converges in tens of CG
# iterations, so hitting CG_MAXITER means the preconditioner is bad; GMRES
# (restarted every GMRES_RESTART iterations) then gets GMRES_MAXITER
# iterations with the same factor before the direct solve. One iteration
# costs roughly 1/250 of an ILU factorisation at 2.5 mm.
CG_MAXITER = 200
GMRES_RESTART = 50
GMRES_MAXITER = 200

# Below this many free DOFs the direct solve is cheaper than building an
# ILU, and exact.
DIRECT_MAX_FREE_DOF = 5000

_METHODS = ("auto", "iterative", "direct", "penalty")


def dirichlet_dofs_and_values(
    bcs: Sequence[BoundaryCondition], n_dof: int
) -> tuple[np.ndarray, np.ndarray]:
    """Unique constrained DOFs (sorted) and their prescribed values.

    Repeated BCs on a DOF are merged when their values agree (to 1e-12
    relative); conflicting values raise ``ValueError``, since elimination
    cannot satisfy both (the penalty method silently averaged them).
    """
    dofs = np.fromiter((bc.dof for bc in bcs), dtype=np.int64, count=len(bcs))
    values = np.fromiter((bc.value for bc in bcs), dtype=float, count=len(bcs))
    if dofs.size == 0:
        return dofs, values
    if dofs.min() < 0 or dofs.max() >= n_dof:
        raise ValueError(f"boundary condition DOF out of range [0, {n_dof})")
    if not np.all(np.isfinite(values)):
        raise ValueError("boundary condition values must be finite")
    unique, first, inverse = np.unique(dofs, return_index=True, return_inverse=True)
    unique_values = values[first]
    spread = np.abs(values - unique_values[inverse])
    tol = 1e-12 * np.maximum(1.0, np.abs(values))
    if np.any(spread > tol):
        bad = int(dofs[np.argmax(spread - tol)])
        conflicting = sorted(set(values[dofs == bad].tolist()))
        raise ValueError(f"conflicting boundary conditions on DOF {bad}: {conflicting}")
    return unique, unique_values


def _ilu(Kff: sp.csc_matrix):
    return spla.spilu(
        Kff, drop_tol=ILU_DROP_TOL, fill_factor=ILU_FILL_FACTOR, permc_spec="MMD_AT_PLUS_A"
    )


def _direct(Kff: sp.csc_matrix, rhs: np.ndarray) -> np.ndarray:
    # Minimum-degree ordering on K + K^T suits the symmetric FE stiffness;
    # SuperLU's default COLAMD ordering fills in several times more.
    lu = spla.splu(Kff, permc_spec="MMD_AT_PLUS_A", options={"SymmetricMode": True})
    return lu.solve(rhs)


def _krylov(
    Kff: sp.csc_matrix, rhs: np.ndarray, ilu, rtol: float, maxiter: int, solver: str = "cg"
) -> tuple[np.ndarray | None, int]:
    """CG (or GMRES) preconditioned by ``ilu``; returns ``(x, iterations)``,
    with ``x`` ``None`` unless the true residual meets ``rtol`` (10x slack
    over the solver's own recurrence estimate)."""
    n_iter = 0

    def _count(_arg) -> None:
        nonlocal n_iter
        n_iter += 1

    M = spla.LinearOperator(Kff.shape, matvec=ilu.solve, dtype=float)
    if solver == "cg":
        x, info = spla.cg(Kff, rhs, rtol=rtol, atol=0.0, maxiter=maxiter, M=M, callback=_count)
    else:
        # GMRES does not need the preconditioner to be SPD, so it rescues a
        # factor CG stalls on. Its maxiter counts restart cycles.
        restart = min(GMRES_RESTART, maxiter)
        x, info = spla.gmres(
            Kff,
            rhs,
            rtol=rtol,
            atol=0.0,
            restart=restart,
            maxiter=max(1, maxiter // restart),
            M=M,
            callback=_count,
            callback_type="pr_norm",
        )
    if info != 0 or not np.all(np.isfinite(x)):
        return None, n_iter
    residual = np.linalg.norm(rhs - Kff @ x)
    if residual > 10.0 * rtol * np.linalg.norm(rhs):
        return None, n_iter
    return x, n_iter


def solve_spd_reduced(
    Kff: sp.spmatrix,
    rhs: np.ndarray,
    *,
    method: str = "auto",
    rtol: float = DEFAULT_RTOL,
    maxiter: int = CG_MAXITER,
) -> np.ndarray:
    """Solve the reduced SPD system ``Kff x = rhs``.

    ``method="iterative"`` runs ILU-preconditioned CG; when CG does not
    converge to ``rtol`` within ``maxiter`` iterations, GMRES retries with
    the same ILU (``GMRES_MAXITER`` iterations), and when that fails too
    (or the ILU cannot be built) a direct sparse LU solves the system.
    ``"direct"`` goes straight to the LU; ``"auto"`` picks direct below
    ``DIRECT_MAX_FREE_DOF`` unknowns and iterative above.
    """
    Kff = sp.csc_matrix(Kff)
    rhs = np.asarray(rhs, dtype=float)
    n = Kff.shape[0]
    if n == 0:
        return np.zeros(0)
    if method == "auto":
        method = "direct" if n <= DIRECT_MAX_FREE_DOF else "iterative"
    if method == "direct":
        return _direct(Kff, rhs)
    if method != "iterative":
        raise ValueError(f"unknown method {method!r}")
    if not np.any(rhs):
        return np.zeros(n)

    try:
        ilu = _ilu(Kff)
    except RuntimeError as exc:  # e.g. an exactly singular pivot
        _log.info("ILU factorisation failed (%s); using the direct solver", exc)
        return _direct(Kff, rhs)
    x, n_iter = _krylov(Kff, rhs, ilu, rtol, maxiter)
    if x is None:
        _log.info("ILU-CG did not reach rtol=%.0e in %d iterations; trying GMRES", rtol, n_iter)
        x, n_iter = _krylov(Kff, rhs, ilu, rtol, min(maxiter, GMRES_MAXITER), solver="gmres")
        if x is None:
            _log.info(
                "ILU-GMRES did not reach rtol=%.0e in %d iterations; using the direct solver",
                rtol,
                n_iter,
            )
            del ilu
            return _direct(Kff, rhs)
        _log.debug("ILU-GMRES: %d iterations, %d DOFs", n_iter, n)
        return x
    _log.debug("ILU-CG: %d iterations, %d DOFs", n_iter, n)
    return x


def solve_linear_static(
    elements: Sequence,
    element_dof_maps: Sequence[np.ndarray],
    n_dof: int,
    bcs: Sequence[BoundaryCondition],
    external_force: np.ndarray | None = None,
    penalty: float = 1.0e10,
    *,
    method: str = "auto",
    rtol: float = DEFAULT_RTOL,
    maxiter: int = CG_MAXITER,
) -> np.ndarray:
    """Assemble K, impose the Dirichlet BCs, and solve K u = F.

    Parameters
    ----------
    elements : sequence of FE elements exposing `.stiffness_matrix() -> (24,24)`.
    element_dof_maps : global DOF indices per element.
    n_dof : total number of global DOFs.
    bcs : Dirichlet boundary conditions. Repeated entries on a DOF must
        agree (see ``dirichlet_dofs_and_values``).
    external_force : optional applied nodal force vector (length n_dof, default zero).
    penalty : penalty scalar, used only by ``method="penalty"``.
    method : ``"auto"`` (default) / ``"iterative"`` / ``"direct"`` eliminate
        the constrained DOFs and solve the reduced system with
        ``solve_spd_reduced``. ``"penalty"`` is the previous path (penalty
        on the diagonal, ``spsolve``), kept for comparison.
    rtol, maxiter : relative-residual tolerance and CG iteration cap for the
        iterative solve (GMRES gets ``min(maxiter, GMRES_MAXITER)``).

    Returns
    -------
    np.ndarray
        Displacement vector of length n_dof; constrained DOFs hold their
        prescribed values exactly (except with ``method="penalty"``).
    """
    if method not in _METHODS:
        raise ValueError(f"unknown method {method!r}; valid options are {list(_METHODS)}")
    K = assemble_global_stiffness(elements, element_dof_maps, n_dof)
    F = (
        np.zeros(n_dof)
        if external_force is None
        else np.asarray(external_force, dtype=float).copy()
    )
    if method == "penalty":
        K_mod, F_mod = apply_dirichlet_penalty(K, F, bcs, penalty=penalty)
        return np.asarray(spla.spsolve(K_mod.tocsc(), F_mod))

    dofs, values = dirichlet_dofs_and_values(bcs, n_dof)
    u = np.zeros(n_dof)
    u[dofs] = values
    is_free = np.ones(n_dof, dtype=bool)
    is_free[dofs] = False
    free = np.flatnonzero(is_free)

    K_rows = K.tocsr()[free]
    del K
    rhs = F[free] - K_rows[:, dofs] @ values
    Kff = K_rows[:, free].tocsc()
    del K_rows
    u[free] = solve_spd_reduced(Kff, rhs, method=method, rtol=rtol, maxiter=maxiter)
    return u
