"""Olsson quasi-static impact threshold load and onset energy.

Closed-form threshold load from plate bending + delamination fracture balance:

    Pc = pi * sqrt(8 * G_IIc * D_eff / 9)
    E_onset = Pc^2 / (2 * k_cb)
    k_cb = 1 / (1/k_bending + 1/k_contact_linearized_at_Pc)

Reference: Olsson (2001), Composites Part A, 32(9); Olsson (2010), IJSS 47(21).
"""

from __future__ import annotations

import math

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate

NAVIER_N: int = 11  # Navier series truncation (N x N modes)


def _k_bending_ssss(
    lam: Laminate,
    pan: PanelGeometry,
    x0: float,
    y0: float,
    n_modes: int = NAVIER_N,
) -> float:
    """Navier series point-load stiffness of a simply-supported rectangular
    orthotropic plate at (x0, y0). Returns N/mm.

    w(x0,y0)/P = (4 / (a*b)) * sum_{m,n} sin^2(m*pi*x0/a) sin^2(n*pi*y0/b) / D_mn
    D_mn = D11 * (m*pi/a)^4 + 2*(D12 + 2*D66) * (m*pi/a)^2 * (n*pi/b)^2
           + D22 * (n*pi/b)^4
    k_bending = 1 / (w/P)
    """
    _, _, D = lam.abd_matrices()
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]
    a, b = pan.Lx_mm, pan.Ly_mm
    w_over_P = 0.0
    for m in range(1, n_modes + 1):
        for n in range(1, n_modes + 1):
            sin_mx = math.sin(m * math.pi * x0 / a)
            sin_ny = math.sin(n * math.pi * y0 / b)
            Dmn = (
                D11 * (m * math.pi / a) ** 4
                + 2 * (D12 + 2 * D66) * (m * math.pi / a) ** 2 * (n * math.pi / b) ** 2
                + D22 * (n * math.pi / b) ** 4
            )
            w_over_P += (sin_mx * sin_ny) ** 2 / Dmn
    w_over_P *= 4.0 / (a * b)
    return 1.0 / w_over_P


def _k_contact_hertz_linearized(lam: Laminate, imp: ImpactorGeometry, P: float) -> float:
    """Linear-equivalent Hertzian contact stiffness at load P (N), returns N/mm.

    Nonlinear Hertz: P = k_nl * delta^(3/2), k_nl = (4/3) * sqrt(R) * E_eff.
    Linearized secant stiffness at P:
        delta = (P / k_nl)^(2/3);   k_lin = P / delta = k_nl^(2/3) * P^(1/3).
    """
    R = imp.diameter_mm / 2.0
    E_steel = 200e3  # MPa = N/mm^2
    nu_steel = 0.3
    E_plate = lam.material.E22
    nu_plate = 0.3
    inv_E = (1 - nu_steel**2) / E_steel + (1 - nu_plate**2) / E_plate
    E_eff = 1.0 / inv_E
    k_nl = (4.0 / 3.0) * math.sqrt(R) * E_eff
    return (k_nl ** (2.0 / 3.0)) * (P ** (1.0 / 3.0))


def threshold_load(lam: Laminate, pan: PanelGeometry, imp: ImpactorGeometry) -> float:
    """Olsson damage-threshold load Pc (N). Uses geometric-mean flexural rigidity D_eff."""
    D_eff = lam.flexural_rigidity_Deff()  # N*mm
    G_IIc = lam.material.G_IIc  # N/mm
    return math.pi * math.sqrt(8 * G_IIc * D_eff / 9.0)


def onset_energy(
    lam: Laminate,
    pan: PanelGeometry,
    imp: ImpactorGeometry,
    location_xy_mm: tuple[float, float] | None = None,
) -> float:
    """Impact energy (J) at which BVID damage onsets."""
    if location_xy_mm is None:
        location_xy_mm = (pan.Lx_mm / 2, pan.Ly_mm / 2)
    x0, y0 = location_xy_mm
    Pc = threshold_load(lam, pan, imp)  # N
    k_b = _k_bending_ssss(lam, pan, x0, y0, n_modes=NAVIER_N)  # N/mm
    k_c = _k_contact_hertz_linearized(lam, imp, P=Pc)  # N/mm (linearized at Pc)
    k_cb = 1.0 / (1.0 / k_b + 1.0 / k_c)  # N/mm
    E_mJ = Pc**2 / (2.0 * k_cb)  # N*mm = mJ
    return E_mJ * 1e-3  # J
