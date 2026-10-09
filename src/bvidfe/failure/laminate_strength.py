"""Unnotched laminate strength under uniaxial in-plane load.

Ply-discount last-ply failure on classical lamination theory (Daniel & Ishai,
*Engineering Mechanics of Composite Materials*, 2006, sec. 7.6). The
laminate is loaded by a uniaxial average stress sigma_x. Each ply's
material-frame stresses are checked with the maximum-stress criterion:

* matrix failure (sigma_2 against Yt/Yc, |tau_12| against S12) cracks the
  ply: its E22 and G12 drop to ``PLY_DISCOUNT_FACTOR`` of their values and
  the load it shed is redistributed at the same applied stress;
* fibre failure (sigma_1 against Xt/Xc) in any ply fails the laminate;
* a laminate whose plies have all cracked has nothing left to stiffen it
  and fails at that load (a [90]n laminate fails at Yt / Yc).

Plies share the laminate strain, not its stress, so the 0 deg plies carry
most of an axial load and break first. The thickness-weighted average of the
ply strengths this replaces credited the off-axis plies with strength they
never reach: it was 1.17-2.85x the NCAMP unnotched strengths of AS4/8552 and
IM7/8552 laminates, which this model matches within 0.83-1.18.

Membrane response only: bending-extension coupling of unsymmetric layups is
neglected, as elsewhere in the closed-form tiers.
"""

from __future__ import annotations

import math

import numpy as np

from bvidfe._types import LoadingMode
from bvidfe.core.laminate import Laminate, _reduced_stiffness, _transform_reduced_stiffness

#: Fraction of E22 and G12 a ply keeps once its matrix has failed.
PLY_DISCOUNT_FACTOR = 0.01


def _ply_stresses_per_unit_load(lam: Laminate, cracked: list[bool], sign: float) -> np.ndarray:
    """(n_plies, 3) material-frame (sigma_1, sigma_2, tau_12) per 1 MPa of
    laminate average stress sigma_x = sign * 1."""
    m = lam.material
    Q_intact = _reduced_stiffness(m)
    Q_cracked = Q_intact.copy()
    Q_cracked[1, :2] *= PLY_DISCOUNT_FACTOR
    Q_cracked[0, 1] *= PLY_DISCOUNT_FACTOR
    Q_cracked[2, 2] *= PLY_DISCOUNT_FACTOR

    Qs = [Q_cracked if c else Q_intact for c in cracked]
    A = sum(
        _transform_reduced_stiffness(Q, math.radians(theta)) * t
        for Q, theta, t in zip(Qs, lam.layup_deg, lam.ply_thicknesses_mm)
    )
    eps_x, eps_y, gamma_xy = np.linalg.solve(A, np.array([sign * lam.thickness_mm, 0.0, 0.0]))

    out = np.empty((lam.n_plies, 3))
    for k, (Q, theta) in enumerate(zip(Qs, lam.layup_deg)):
        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        eps_1 = c * c * eps_x + s * s * eps_y + c * s * gamma_xy
        eps_2 = s * s * eps_x + c * c * eps_y - c * s * gamma_xy
        gamma_12 = 2 * c * s * (eps_y - eps_x) + (c * c - s * s) * gamma_xy
        out[k] = Q @ np.array([eps_1, eps_2, gamma_12])
    return out


def _load_to(strength_pos: float, strength_neg: float, stress: np.ndarray) -> np.ndarray:
    """Load factor at which ``stress`` (per unit load) reaches its strength."""
    with np.errstate(divide="ignore"):
        return np.where(
            stress > 0,
            strength_pos / stress,
            np.where(stress < 0, strength_neg / -stress, np.inf),
        )


def unnotched_strength(lam: Laminate, loading: LoadingMode) -> float:
    """Laminate strength (MPa, positive) under uniaxial x tension or compression.

    See the module docstring for the model.
    """
    m = lam.material
    sign = 1.0 if loading == "tension" else -1.0
    cracked = [False] * lam.n_plies
    load = 0.0
    while True:
        s = _ply_stresses_per_unit_load(lam, cracked, sign)
        fibre = _load_to(m.Xt, m.Xc, s[:, 0])
        matrix = np.minimum(_load_to(m.Yt, m.Yc, s[:, 1]), _load_to(m.S12, m.S12, s[:, 2]))
        matrix[cracked] = np.inf
        if fibre.min() <= matrix.min():
            return max(float(fibre.min()), load)
        # Crack every ply that fails at this load, then re-solve at the same
        # applied stress: the shed load may crack or break further plies.
        load = max(float(matrix.min()), load)
        for k in np.flatnonzero(matrix <= load * (1.0 + 1e-9)):
            cracked[k] = True
        if all(cracked):
            return load
