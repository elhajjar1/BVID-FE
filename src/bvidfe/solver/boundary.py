"""Displacement-controlled boundary conditions.

``solve_linear_static`` imposes these by eliminating the constrained DOFs;
``apply_dirichlet_penalty`` is the earlier penalty path, still used by its
``method="penalty"`` option.

Standard patterns:
- uniaxial_x_bcs: prescribed end shortening (or extension) between rigid,
  frictionless platens on the x faces, ``u_z = 0`` on the mid-surface node
  layer and ``u_y = 0`` at one node. The sign of ``applied_strain`` selects
  compression (negative) or tension (positive) — the BC topology is
  identical either way.
- apply_dirichlet_penalty: multiply K[i,i] by penalty and F[i] = penalty * value.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import numpy as np
import scipy.sparse as sp


@dataclass
class BoundaryCondition:
    """A single prescribed-displacement boundary condition on DOF `dof`, value `value`."""

    dof: int
    value: float


def apply_dirichlet_penalty(
    K: sp.spmatrix,
    F: np.ndarray,
    bcs: Sequence[BoundaryCondition],
    penalty: float = 1.0e10,
) -> tuple[sp.csc_matrix, np.ndarray]:
    """Apply Dirichlet BCs to (K, F) via the penalty method.

    Returns a new (K_mod, F_mod) without mutating the inputs.

    The diagonal is updated value-only: every DOF is self-coupled in an
    assembled FE stiffness matrix, so ``K[d, d]`` already exists and adding
    ``penalty`` via ``K + sp.diags(p)`` introduces no new sparsity (avoids
    the SparseEfficiencyWarning storm of scalar ``K_csr[d, d] = ...``
    assignment). Numerically identical to the previous per-DOF loop.
    """
    n_dof = K.shape[0]
    F_out = F.copy()
    bc_dofs = np.fromiter((bc.dof for bc in bcs), dtype=int, count=len(bcs))
    bc_values = np.fromiter((bc.value for bc in bcs), dtype=float, count=len(bcs))

    p = np.zeros(n_dof)
    # np.add.at accumulates duplicates, matching the loop's repeated
    # K_csr[d, d] += penalty / F_out[d] += penalty * value for repeated DOFs.
    np.add.at(p, bc_dofs, penalty)
    np.add.at(F_out, bc_dofs, penalty * bc_values)

    K_mod = (K + sp.diags(p, format="csc")).tocsc()
    return K_mod, F_out


def _nodes_on_plane(
    node_coords: np.ndarray, axis: int, coord: float, tol: float = 1e-9
) -> np.ndarray:
    """Return indices of nodes with coordinate along axis close to `coord`."""
    return np.where(np.abs(node_coords[:, axis] - coord) < tol)[0]


def mid_surface_z(node_coords: np.ndarray) -> float:
    """z of the node layer nearest the laminate mid-surface: the mid-surface
    itself for an even number of element layers, the lower of the two
    layers half an element from it for an odd number."""
    z_layers = np.unique(node_coords[:, 2])
    distance = np.abs(z_layers - 0.5 * (z_layers[0] + z_layers[-1]))
    tol = 1e-9 * max(1.0, float(z_layers[-1] - z_layers[0]))
    return float(z_layers[np.flatnonzero(distance <= distance.min() + tol)[0]])


def uniaxial_x_bcs(
    node_coords: np.ndarray,
    applied_strain: float,
    boundary: str = "simply_supported",
) -> List[BoundaryCondition]:
    """Build BCs for a uniaxial test along x (compression or tension).

    The sign of ``applied_strain`` is the only difference between the
    compression and tension load cases — the constrained-DOF topology is
    identical, so a single builder serves both:

    - every ``x_min`` face node: ``u_x = 0``; every ``x_max`` face node:
      ``u_x = applied_strain * Lx`` (negative = compression, positive =
      tension). These are rigid, frictionless platens: ``u_y`` and ``u_z``
      stay free on the loaded faces.
    - every node of the mid-surface layer (``mid_surface_z``): ``u_z = 0``.
    - one node, on the ``x_min`` face at the mid-surface nearest
      ``y = Ly / 2``: ``u_y = 0``.

    No other DOF is constrained, and no DOF twice. With the platens (which
    fix ``u_x`` and so the rotations about y and z) the last two remove the
    remaining rigid-body modes: z translation and rotation about x (the
    mid-surface layer), y translation (the single node, whose reaction is
    zero).

    Why the mid-surface layer. A uniform stress state (each ply's free
    Poisson expansion in y and through the thickness) satisfies every one of
    these restraints, so an undamaged unidirectional coupon fails exactly at
    its ply strength and an undamaged laminate at its CLT first-ply failure
    load, apart from the physical free-edge stresses at ``y = 0`` and
    ``y = Ly``. When the laminate and its damage are symmetric about the
    mid-plane (a symmetric layup: the damage zone softens every ply), the
    exact solution has ``u_z = 0`` on the mid-plane, so the restraint is a
    symmetry condition and carries no load. An unsymmetric layup is held
    flat instead of bending, and with an odd number of element layers the
    restrained layer lies half an element off the mid-plane: the solve is
    a membrane (pre-buckling) state. Bending is left to the closed-form
    channels: delamination growth by sublaminate buckling and Rayleigh-Ritz
    panel buckling, which use ``boundary``.

    ``boundary`` is accepted for the callers' signature but does not change
    the constraints: in the membrane solve of a flat panel the out-of-plane
    edge supports carry no load. Supporting only the mid-surface lines of
    the four edges (a simple support) and leaving the panel free to bend
    changed the knockdowns of the calibration coupons and the 8-ply panels
    by under 1%, and the reduced system lost the conditioning the
    ILU-preconditioned Krylov solve needs (a 40-ply coupon's reference
    solve at 5 mm: 48 s by ILU-GMRES, 183 s by the direct fallback; the
    whole analyses took 2-4x as long). Its plate bending would also rest
    on fully integrated C3D8 bricks a few tenths of a millimetre thick and
    several millimetres wide, which lock in bending.

    The previous builder held ``u_z = 0`` on every node of the ``z_min``
    face, of the two loaded faces (and of the ``y`` faces for
    ``"clamped"``), and ``u_y = 0`` on the whole ``y_min`` face. On the
    loaded and ``y`` faces that blocked the plies' through-thickness
    Poisson expansion (``nu23`` = 0.3-0.7 for the presets) and so raised
    the stress next to the edges: a ``[90]_4`` coupon failed 3-4% below
    ``Yc / E22`` in compression, and the undamaged 8-ply QI panel 12.6%
    below its CLT first-ply failure load in tension. The ``boundary``
    dependence it gave the first-ply failure strain (clamped < simply
    supported < free) came from those faces alone, since the ``z_min``
    face already held the panel flat; a face of ``u_z = 0`` does not fix
    the edge rotation either (the face rotates about its line, ``u_y``
    varying through the thickness).
    """
    del boundary  # see the docstring: the membrane solve does not use it
    Lx = node_coords[:, 0].max() - node_coords[:, 0].min()
    xmin_nodes = _nodes_on_plane(node_coords, 0, node_coords[:, 0].min())
    xmax_nodes = _nodes_on_plane(node_coords, 0, node_coords[:, 0].max())
    mid_nodes = _nodes_on_plane(node_coords, 2, mid_surface_z(node_coords))

    y = node_coords[:, 1]
    x_min_mid = np.intersect1d(xmin_nodes, mid_nodes)
    y_node = int(x_min_mid[np.argmin(np.abs(y[x_min_mid] - 0.5 * (y.min() + y.max())))])

    bcs: List[BoundaryCondition] = []
    for n in xmin_nodes:
        bcs.append(BoundaryCondition(dof=3 * int(n) + 0, value=0.0))
    for n in xmax_nodes:
        bcs.append(BoundaryCondition(dof=3 * int(n) + 0, value=applied_strain * Lx))
    bcs.append(BoundaryCondition(dof=3 * y_node + 1, value=0.0))
    for n in mid_nodes:
        bcs.append(BoundaryCondition(dof=3 * int(n) + 2, value=0.0))
    return bcs
