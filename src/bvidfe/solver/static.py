"""Linear static finite-element solver.

Dirichlet BCs are imposed by eliminating the constrained DOFs: with ``c``
the constrained and ``f`` the free DOFs, the solve is

    K_ff u_f = F_f - K_fc u_c

and ``u_c`` is copied into the result exactly. ``K_ff`` is symmetric
positive definite. It is solved by GMRES preconditioned with an incomplete
LU factorisation (``spilu``), and a direct sparse LU of ``K_ff`` is the
fallback, so every returned solution meets the tolerance.

Why GMRES and not conjugate gradients. ``spilu`` is an incomplete LU
(SuperLU's ILUTP), not an incomplete Cholesky: it drops entries of ``L``
and of ``U`` separately, so ``L U`` is not symmetric although ``K_ff`` is
(to 6e-17), and neither is the preconditioner ``M = (L U)^-1``:
``|v . M u - u . M v| / |v . M u|`` is 0.3-0.4 for random ``u``, ``v`` on
the 8-ply 5 mm systems, with no row pivoting involved. CG needs a
symmetric positive definite preconditioner; GMRES does not. Before the
mid-surface ``u_z`` restraint (``uniaxial_x_bcs``) the reduced system was
well enough conditioned for ILU-CG to converge regardless (undamaged 8-ply
150 x 100 mm panel at 5 mm: Jacobi-scaled condition number 7.5e4, 101
iterations); with it the condition number is 1.3e6 (1.7e6 damaged) and
ILU-CG did not converge in 200 iterations on any 5 mm system measured (8
to 40 plies, damaged or not, compression and tension; on the 8-ply panel
it stalls at a relative residual of 1e-6 to 1e-5), while GMRES with the
same factor reaches 1e-10 in 26-42 iterations. CG with the symmetric part of the same preconditioner,
``(M + M^T) / 2``, converges in as many iterations as GMRES, which
confirms the asymmetry as the cause, but needs two ILU solves an
iteration; Jacobi and 3 x 3 nodal block-Jacobi CG need ~2,650 iterations
at 5 mm (3 s against 1 s for ILU-GMRES); a direct LU takes 3 s. GMRES, the
fallback after 200 wasted CG iterations before, is now the only
Krylov method.

Measured on an 8-ply 150 x 100 mm panel at a 2.5 mm mesh (64k free
DOFs), the original path (1e10 penalty on the diagonal, then ``spsolve``
with SuperLU's default COLAMD ordering) took ~110 s and ~2.9 GB per solve;
ILU-GMRES takes ~7 s in about 44 iterations (5 s of it the factorisation).
"""

from __future__ import annotations

import logging
import time
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

# ILU preconditioner settings. ILU-GMRES on the damaged 8-ply 150 x 100 mm
# CAI system at 2.5 mm (64k free DOFs, 4.5M nonzeros, one thread,
# GMRES_RESTART = 100):
#   drop_tol  nnz(L+U)  factor  GMRES iterations  total
#   1e-2       5.8M     3.4 s   no convergence in 200 (residual 4e-10)
#   3e-3       9.0M     4.1 s   70                6.5 s
#   1e-3      12.6M     5.0 s   44                6.8 s
#   3e-4      17.2M     6.2 s   28                8.0 s
#   1e-4      22.2M     8.0 s   19               10.0 s
# 3e-3 is within the timing noise of 1e-3 (5 mm 8-ply: 1.0 s both;
# 24-ply NCAMP coupon at 5 mm: 6.4 s against 7.4 s) at 1.5-2x the
# iterations and next to the 1e-2 that does not converge, so 1e-3 is kept.
# SuperLU caps the fill column by column at fill_factor x the column of A
# and drops harder where the cap binds, which can leave an unstable factor
# (on the 2.5 mm panel with the earlier BCs, drop_tol=1e-5 at
# fill_factor=10 left one CG could not use in 300 iterations and GMRES
# needed 230 for, against 7 CG iterations at fill_factor=30); the generous
# fill_factor keeps the cap from binding on finer meshes (1e-3 gives the
# same factor at 10 and 30 here).
ILU_DROP_TOL = 1.0e-3
ILU_FILL_FACTOR = 30.0

# GMRES iteration cap (inner iterations, over restart cycles of
# GMRES_RESTART). A freshly factorised ILU converges in 26-58 iterations
# (5 mm to 2 mm meshes, 8 to 40 plies), so hitting the cap means the
# preconditioner is bad and the direct solve takes over. Restarting every
# 50 iterations cost 30-40% more iterations where more than 50 are needed
# (8-ply at 2 mm: 74 against 56; 16-ply at 2.5 mm: 80 against 58); the
# 100-vector Krylov basis costs 800 bytes per free DOF (120 MB at the
# FE3D_MAX_DOF cap), against ~2.4 kB (~200 nonzeros) for the ILU factor.
GMRES_RESTART = 100
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


def _ilu_gmres(
    Kff: sp.csc_matrix, rhs: np.ndarray, ilu, rtol: float, maxiter: int
) -> tuple[np.ndarray | None, int]:
    """GMRES preconditioned by ``ilu``, restarted every ``GMRES_RESTART``
    iterations; returns ``(x, iterations)``, with ``x`` ``None`` unless the
    true residual meets ``rtol`` (10x slack over GMRES's own estimate)."""
    n_iter = 0

    def _count(_arg) -> None:
        nonlocal n_iter
        n_iter += 1

    M = spla.LinearOperator(Kff.shape, matvec=ilu.solve, dtype=float)
    restart = min(GMRES_RESTART, maxiter)
    # A breakdown (a preconditioner that maps the residual to zero; with
    # CG, r . M r = 0, seen with the ILU on macOS) makes the solver divide
    # by zero. The non-finite result is caught below and the caller moves
    # on to the direct solve, so numpy's RuntimeWarning for it is silenced
    # rather than passed to the user.
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        x, info = spla.gmres(
            Kff,
            rhs,
            rtol=rtol,
            atol=0.0,
            restart=restart,
            maxiter=max(1, maxiter // restart),  # counts restart cycles
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
    maxiter: int = GMRES_MAXITER,
) -> np.ndarray:
    """Solve the reduced SPD system ``Kff x = rhs``.

    ``method="iterative"`` runs ILU-preconditioned GMRES (``maxiter``
    iterations, restarted every ``GMRES_RESTART``); when it does not
    converge to ``rtol`` (or the ILU cannot be built) a direct sparse LU
    solves the system. ``"direct"`` goes straight to the LU; ``"auto"``
    picks direct below ``DIRECT_MAX_FREE_DOF`` unknowns and iterative
    above.

    Logs an INFO line on entry and one on exit (the solver that produced
    the answer, its iterations and the wall time) on ``bvidfe.fe3d.solver``.
    """
    Kff = sp.csc_matrix(Kff)
    rhs = np.asarray(rhs, dtype=float)
    n = Kff.shape[0]
    if method == "auto":
        method = "direct" if n <= DIRECT_MAX_FREE_DOF else "iterative"
    t0 = time.perf_counter()
    _log.info("solve start: method=%s, %d free DOFs, %d nonzeros", method, n, Kff.nnz)
    x, how = _solve_spd_reduced(Kff, rhs, method, rtol, maxiter)
    _log.info("solve done: %s (%.2fs)", how, time.perf_counter() - t0)
    return x


def _solve_spd_reduced(
    Kff: sp.csc_matrix, rhs: np.ndarray, method: str, rtol: float, maxiter: int
) -> tuple[np.ndarray, str]:
    """``solve_spd_reduced`` with ``"auto"`` resolved: the solution and a
    description of the solver that produced it, for the log."""
    n = Kff.shape[0]
    if n == 0:
        return np.zeros(0), "no free DOFs"
    if method == "direct":
        return _direct(Kff, rhs), "direct LU"
    if method != "iterative":
        raise ValueError(f"unknown method {method!r}")
    if not np.any(rhs):
        return np.zeros(n), "zero load, nothing to solve"

    try:
        ilu = _ilu(Kff)
    except RuntimeError as exc:  # e.g. an exactly singular pivot
        _log.info("ILU factorisation failed (%s); using the direct solver", exc)
        return _direct(Kff, rhs), "direct LU (ILU factorisation failed)"
    x, n_iter = _ilu_gmres(Kff, rhs, ilu, rtol, maxiter)
    if x is not None:
        return x, f"ILU-GMRES, {n_iter} iterations"
    _log.info(
        "ILU-GMRES did not reach rtol=%.0e in %d iterations; using the direct solver",
        rtol,
        n_iter,
    )
    del ilu
    return _direct(Kff, rhs), f"direct LU (after {n_iter} ILU-GMRES iterations)"


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
    maxiter: int = GMRES_MAXITER,
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
    rtol, maxiter : relative-residual tolerance and GMRES iteration cap for
        the iterative solve.

    Returns
    -------
    np.ndarray
        Displacement vector of length n_dof; constrained DOFs hold their
        prescribed values exactly (except with ``method="penalty"``).
    """
    if method not in _METHODS:
        raise ValueError(f"unknown method {method!r}; valid options are {list(_METHODS)}")
    t0 = time.perf_counter()
    _log.info("assembly start: %d elements, %d DOFs", len(elements), n_dof)
    K = assemble_global_stiffness(elements, element_dof_maps, n_dof)
    _log.info("assembly done: %d nonzeros (%.2fs)", K.nnz, time.perf_counter() - t0)
    F = (
        np.zeros(n_dof)
        if external_force is None
        else np.asarray(external_force, dtype=float).copy()
    )
    t0 = time.perf_counter()
    _log.info(
        "BCs start: %d conditions, %s",
        len(bcs),
        "penalty" if method == "penalty" else "elimination",
    )
    if method == "penalty":
        K_mod, F_mod = apply_dirichlet_penalty(K, F, bcs, penalty=penalty)
        _log.info(
            "BCs done: %.0e penalty on the diagonal (%.2fs)", penalty, time.perf_counter() - t0
        )
        t0 = time.perf_counter()
        _log.info("solve start: method=penalty, %d DOFs", n_dof)
        u = np.asarray(spla.spsolve(K_mod.tocsc(), F_mod))
        _log.info("solve done: spsolve (%.2fs)", time.perf_counter() - t0)
        return u

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
    _log.info(
        "BCs done: %d constrained DOFs eliminated, %d free (%.2fs)",
        dofs.size,
        free.size,
        time.perf_counter() - t0,
    )
    u[free] = solve_spd_reduced(Kff, rhs, method=method, rtol=rtol, maxiter=maxiter)
    return u
