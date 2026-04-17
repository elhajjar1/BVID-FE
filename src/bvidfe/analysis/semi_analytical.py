"""Semi-analytical tier: sublaminate Rayleigh-Ritz buckling + critical interface scoring.

The ellipse is approximated as its enclosing simply-supported rectangle
(2a x 2b in the panel frame). For orthotropic simply-supported rectangles
under uniaxial compression, the sine basis is exact and the eigenvalue is
closed-form — we minimize over integer modes (m, n) in [1..5] x [1..5].

Sublaminate selection: the plies above the delaminated interface form the
thinner buckling sublaminate (closer to the impact face for interfaces in
the upper half of the laminate). We always use the smaller of the two
sublaminates because it buckles first.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from bvidfe.core.laminate import Laminate
from bvidfe.core.material import OrthotropicMaterial
from bvidfe.damage.state import DamageState, DelaminationEllipse


def _sublaminate_D_matrix(
    material: OrthotropicMaterial,
    sub_layup_deg: list[float],
    ply_thickness_mm: float,
) -> np.ndarray:
    """CLT D matrix (3, 3) for a sublaminate. z origin at sublaminate midplane."""
    sub_lam = Laminate(material, sub_layup_deg, ply_thickness_mm)
    _, _, D = sub_lam.abd_matrices()
    return D


def sublaminate_buckling_load(lam: Laminate, ellipse: DelaminationEllipse) -> float:
    """Critical buckling force per unit width (N/mm) for the sublaminate above
    the given ellipse's interface, modeled as a simply-supported rectangle
    with semi-axes = ellipse major/minor.
    """
    i = ellipse.interface_index
    full_layup = lam.layup_deg

    # Choose the thinner sublaminate between "above" (plies 0..i) and "below" (plies i+1..)
    upper_layup = full_layup[: i + 1]
    lower_layup = full_layup[i + 1 :]
    sub_layup = upper_layup if len(upper_layup) <= len(lower_layup) else lower_layup
    if len(sub_layup) == 0:
        return float("inf")

    D = _sublaminate_D_matrix(lam.material, sub_layup, lam.ply_thickness_mm)
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]

    # Rectangle dimensions (panel frame). Ellipse semi-axes = a, b.
    a = ellipse.major_mm
    b = ellipse.minor_mm
    if a <= 0 or b <= 0:
        return float("inf")

    # Minimum over (m, n) in 1..5 for uniaxial compression N0_x:
    # N_cr(m,n) = (pi^2 / a^2) * [D11*m^4 + 2*(D12+2*D66)*(m*a/b)^2*n^2 + D22*(a*n/b)^4] / m^2
    pi2 = math.pi * math.pi
    best = float("inf")
    for m_mode in range(1, 6):
        for n_mode in range(1, 6):
            num = (
                D11 * m_mode**4
                + 2.0 * (D12 + 2.0 * D66) * (m_mode * a / b) ** 2 * n_mode**2
                + D22 * (a * n_mode / b) ** 4
            )
            N_mn = (pi2 / a**2) * num / m_mode**2
            if N_mn < best:
                best = N_mn
    return best


def find_critical_interface(damage: DamageState, lam: Laminate) -> Optional[int]:
    """Return the interface index that would fail first under compression.

    Scoring: max_area_i * max(|z_upper_i|, |z_lower_i|), where z is distance
    from interface to the top/bottom laminate surface. Largest wins.
    """
    if not damage.delaminations:
        return None
    n_plies = len(lam.layup_deg)
    h = lam.ply_thickness_mm
    per_iface_max_area: dict[int, float] = {}
    for e in damage.delaminations:
        per_iface_max_area[e.interface_index] = max(
            per_iface_max_area.get(e.interface_index, 0.0), e.area_mm2
        )

    best_idx: Optional[int] = None
    best_score = -1.0
    for idx, area in per_iface_max_area.items():
        z_upper = (idx + 1) * h  # distance from top of laminate (plies 0..idx above)
        z_lower = (n_plies - idx - 1) * h
        score = area * max(z_upper, z_lower)
        if score > best_score:
            best_score = score
            best_idx = idx
    return best_idx
