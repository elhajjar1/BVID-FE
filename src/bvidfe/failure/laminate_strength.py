"""Unnotched laminate strength under a proportional in-plane load.

Ply-discount last-ply failure on classical lamination theory (Daniel & Ishai,
*Engineering Mechanics of Composite Materials*, 2006, sec. 7.6). The
laminate is loaded by a uniaxial average stress sigma_x or, for biaxial
loading (``unnotched_strength(..., transverse_stress_ratio)``,
``unnotched_strength_along``), by a proportional average stress
``lambda * (sigma_x, sigma_y, tau_xy)``. Each ply's
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
IM7/8552 laminates, which this model matches within 0.83-1.14.

Membrane response only: bending-extension coupling of unsymmetric layups is
neglected, as elsewhere in the closed-form tiers.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from bvidfe._types import LoadingMode
from bvidfe.core.laminate import Laminate, _reduced_stiffness, _transform_reduced_stiffness

#: Fraction of E22 and G12 a ply keeps once its matrix has failed.
PLY_DISCOUNT_FACTOR = 0.01


def _ply_stresses_per_unit_direction(
    lam: Laminate, cracked: list[bool], direction: np.ndarray
) -> np.ndarray:
    """(n_plies, 3) material-frame (sigma_1, sigma_2, tau_12) per unit
    multiplier of the laminate average stress ``(sigma_x, sigma_y, tau_xy)
    = direction``.

    The right-hand side is ``h * direction``. For a uniaxial direction
    ``(sign, 0.0, 0.0)`` that is the same IEEE vector ``(sign * h, 0.0,
    0.0)`` the uniaxial model has always solved with, so its results are
    unchanged bit for bit."""
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
    eps_x, eps_y, gamma_xy = np.linalg.solve(A, lam.thickness_mm * direction)

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


def unnotched_strength(
    lam: Laminate, loading: LoadingMode, transverse_stress_ratio: float = 0.0
) -> float:
    """Laminate strength (MPa, positive): ``|sigma_x|`` at last-ply failure
    along the proportional in-plane load ``(sigma_x, sigma_y, tau_xy) =
    lambda * (s, s * beta, 0)``, ``s = +1`` for ``"tension"`` and ``-1``
    for ``"compression"``, ``beta = transverse_stress_ratio`` (signed
    ``sigma_y / sigma_x``).

    ``beta == 0`` (the default) is uniaxial x tension or compression, the
    model as it has always been. See the module docstring for the model and
    :func:`unnotched_strength_along` for a general direction.

    >>> from bvidfe import MATERIAL_LIBRARY
    >>> from bvidfe.core.laminate import Laminate
    >>> qi = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    >>> round(unnotched_strength(qi, "compression"), 1)
    601.3
    >>> round(unnotched_strength(qi, "compression", 1.0), 1)  # equibiaxial
    848.3
    """
    sign = 1.0 if loading == "tension" else -1.0
    if transverse_stress_ratio == 0:
        direction = np.array([sign, 0.0, 0.0])
    else:
        direction = np.array([sign, sign * float(transverse_stress_ratio), 0.0])
    return _ply_discount_multiplier(lam, direction)


def unnotched_strength_along(lam: Laminate, direction: Sequence[float]) -> float:
    """Ply-discount last-ply failure multiplier ``lambda0`` (positive) of the
    laminate along a general proportional in-plane direction ``direction =
    (sigma_x, sigma_y, tau_xy)``: the pristine laminate fails at the average
    stress ``lambda0 * direction``.

    The model is the module docstring's, with the load ``h * direction`` on
    the right-hand side of the CLT solve instead of ``h * (+-1, 0, 0)``, so
    ``unnotched_strength_along(lam, (+-1, 0, 0))`` is ``unnotched_strength``
    exactly. Under biaxial tension a laminate can end with every ply cracked
    before any fibre breaks; it then fails at that load (the uniaxial
    ``[90]n`` rule), which is conservative where a fibre network would
    carry on (see "Biaxial loading" in docs/physics_models.md).

    Raises ``ValueError`` for a zero or non-finite direction.

    >>> from bvidfe import MATERIAL_LIBRARY
    >>> from bvidfe.core.laminate import Laminate
    >>> qi = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    >>> round(unnotched_strength_along(qi, (1.0, -1.0, 0.0)), 1)  # pure shear at 45 deg
    411.6
    """
    d = np.asarray(direction, dtype=float)
    if d.shape != (3,) or not np.all(np.isfinite(d)) or not np.any(d):
        raise ValueError(
            "direction must be three finite numbers (sigma_x, sigma_y, tau_xy), "
            f"not all zero (got {direction!r})"
        )
    return _ply_discount_multiplier(lam, d)


def _ply_discount_multiplier(lam: Laminate, direction: np.ndarray) -> float:
    """Last-ply failure multiplier along ``direction`` (the shared loop of
    :func:`unnotched_strength` and :func:`unnotched_strength_along`)."""
    m = lam.material
    cracked = [False] * lam.n_plies
    load = 0.0
    while True:
        s = _ply_stresses_per_unit_direction(lam, cracked, direction)
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
