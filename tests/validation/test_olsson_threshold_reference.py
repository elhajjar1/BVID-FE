"""Olsson threshold-load + onset-energy validation.

The quasi-static (large-mass) delamination threshold for one delamination
(Davies & Robinson 1992; Suemasu & Majima 1996; Olsson 2001) is
   Pc = pi * sqrt(32 * D* * G_IIc / 3)
   D* = sqrt(D11 * D22 * (A + 1) / 2),  A = (D12 + 2 * D66) / sqrt(D11 * D22)
(Gonzalez et al. 2011, Eq. 8; Salvetti et al. 2016, Eq. 4). These tests
pin the closed form and its shape (toughness scaling, thickness scaling,
panel size invariance), and check the absolute value against a measured
threshold load (Gonzalez et al. 2011, AS4/8552 quasi-isotropic).

The shape relationships follow directly from the closed form and must hold
for every BVID-FE material:
  * Pc scales as sqrt(G_IIc) when D* is held fixed.
  * Pc scales as h^{3/2} under a SELF-SIMILAR thickness scaling, i.e.
    holding the layup fixed and scaling ply thickness (so every ply's
    relative z-position is preserved). Under CLT D_ij integrates
    Q_bar_k * (z_{k+1}^3 - z_k^3)/3, so D* only scales as h^3 when the
    ply positions scale uniformly with h. Concatenating a layup (e.g.
    ``[0,45,-45,90]*2``) does NOT preserve relative positions -- it
    rearranges plies in z and gives an effective exponent below 1.5
    (see issue #44).
  * Pc is a property of the laminate + impactor and is INDEPENDENT of
    panel size, so the same Laminate+ImpactorGeometry should yield the
    same Pc on a 100x100 mm and 200x150 mm panel.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.olsson import threshold_load


def _laminate(material_name: str, layup, t_ply: float = 0.152):
    return Laminate(MATERIAL_LIBRARY[material_name], list(layup), t_ply)


def _panel(Lx: float = 150.0, Ly: float = 100.0):
    return PanelGeometry(Lx_mm=Lx, Ly_mm=Ly, boundary="simply_supported")


def _impactor():
    return ImpactorGeometry(diameter_mm=16.0, shape="hemispherical")


def test_threshold_load_is_positive_and_finite():
    """A standard CFRP layup must produce a positive, finite Pc."""
    Pc = threshold_load(
        _laminate("IM7/8552", [0, 45, -45, 90, 90, -45, 45, 0]), _panel(), _impactor()
    )
    assert math.isfinite(Pc)
    assert Pc > 0


def _D_star(lam: Laminate) -> float:
    _, _, D = lam.abd_matrices()
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]
    D_geo = math.sqrt(D11 * D22)
    return math.sqrt(D11 * D22 * ((D12 + 2.0 * D66) / D_geo + 1.0) / 2.0)


def test_threshold_load_scales_as_sqrt_of_giic():
    """Pc ~ sqrt(G_IIc) when D* is fixed.

    Hold the laminate stack fixed (same E_ij, same ply count) and swap
    materials with different G_IIc only. We use IM7/8552 (G_IIc=0.79)
    and T700/2510 (G_IIc=0.60). The two materials have different E_ij
    too, so we isolate the G_IIc effect by comparing Pc / sqrt(D*)."""
    lam_a = _laminate("IM7/8552", [0, 45, -45, 90, 90, -45, 45, 0])
    lam_b = _laminate("T700/2510", [0, 45, -45, 90, 90, -45, 45, 0])
    pan, imp = _panel(), _impactor()
    Pc_a = threshold_load(lam_a, pan, imp)
    Pc_b = threshold_load(lam_b, pan, imp)
    # Pc / sqrt(D* * G_IIc) must be the same constant for both materials
    # (it equals pi * sqrt(32/3)).
    expected_const = math.pi * math.sqrt(32.0 / 3.0)
    actual_a = Pc_a / math.sqrt(_D_star(lam_a) * lam_a.material.G_IIc)
    actual_b = Pc_b / math.sqrt(_D_star(lam_b) * lam_b.material.G_IIc)
    assert actual_a == pytest.approx(expected_const, rel=1e-9)
    assert actual_b == pytest.approx(expected_const, rel=1e-9)


def test_threshold_load_uses_olsson_D_star_not_geometric_mean():
    """Olsson's D* = sqrt(D11 D22 (A + 1)/2) sets Pc, not sqrt(D11 D22).

    On a cross-ply A = (D12 + 2 D66)/sqrt(D11 D22) is about 0.17, so D* is
    about 0.77 sqrt(D11 D22) and Pc about 12% below the geometric-mean value.
    """
    lam = _laminate("IM7/8552", [0, 90, 0, 90, 90, 0, 90, 0])
    Pc = threshold_load(lam, _panel(), _impactor())
    G = lam.material.G_IIc
    Pc_geo = math.pi * math.sqrt(32.0 * lam.flexural_rigidity_Deff() * G / 3.0)
    assert Pc == pytest.approx(math.pi * math.sqrt(32.0 * _D_star(lam) * G / 3.0), rel=1e-12)
    assert 0.84 < Pc / Pc_geo < 0.92


def test_threshold_load_matches_isotropic_davies_form():
    """For an isotropic plate D* = D = E h^3 / (12 (1 - nu^2)) and the threshold
    is Davies' Pc^2 = 8 pi^2 E h^3 G_IIc / (9 (1 - nu^2))."""
    E, nu, G_IIc, h = 50_000.0, 0.3, 0.8, 4.0
    G = E / (2 * (1 + nu))
    iso = dataclasses.replace(
        MATERIAL_LIBRARY["IM7/8552"], E11=E, E22=E, nu12=nu, G12=G, G13=G, G23=G, G_IIc=G_IIc
    )
    lam = Laminate(iso, [0, 0, 0, 0], h / 4)
    Pc_davies = math.sqrt(8 * math.pi**2 * E * h**3 * G_IIc / (9 * (1 - nu**2)))
    assert threshold_load(lam, _panel(), _impactor()) == pytest.approx(Pc_davies, rel=1e-9)


def test_threshold_load_against_gonzalez_2011_measured_onset():
    """Gonzalez et al., Compos. Sci. Technol. 71 (2011) 805-817: AS4/8552
    [(45/0/-45/90)4]S, 32 plies, 5.8 mm, G_IIc 0.79 N/mm, drop-weight impacts
    to ASTM D7136. The measured delamination threshold load (the load drop),
    averaged over their three impact energies, is 9.94 kN. With the repo's
    AS4/8552 preset the closed form gives 8.58 kN (0.86x); the old
    pi * sqrt(8 G_IIc sqrt(D11 D22) / 9) gave 2.45 kN."""
    half = [45, 0, -45, 90] * 4
    lam = Laminate(MATERIAL_LIBRARY["AS4/8552"], half + half[::-1], 5.8 / 32)
    Pc = threshold_load(lam, _panel(), _impactor())
    assert Pc == pytest.approx(9940.0, rel=0.15)


def test_threshold_load_invariant_to_panel_size():
    """Pc depends on the LAMINATE, not the panel (per the Olsson model)."""
    lam = _laminate("IM7/8552", [0, 45, -45, 90, 90, -45, 45, 0])
    imp = _impactor()
    Pc_small = threshold_load(lam, _panel(100, 80), imp)
    Pc_large = threshold_load(lam, _panel(300, 200), imp)
    assert Pc_small == pytest.approx(Pc_large, rel=1e-12)


def test_threshold_load_scales_with_thickness_to_the_three_halves():
    """Pc scales as h^{3/2} under a self-similar thickness change.

    The closed form gives Pc = pi * sqrt(32 * D* * G_IIc / 3), with D*
    built from D11, D22, D12 and D66. The CLT D-matrix is
    ``D_ij = (1/3) * sum_k Q_bar_ij,k * (z_{k+1}^3 - z_k^3)``, so D*
    scales as h^3 -- and therefore Pc as h^{3/2} -- ONLY when each ply's
    relative through-thickness position is preserved as h changes. That
    means a self-similar scaling: keep the layup fixed and scale the ply
    thickness. Doubling the layup count (e.g. ``[0,45,-45,90] * 2``)
    does NOT preserve relative positions -- it rearranges where the 0/45
    /-45/90 plies sit in z, which changes their D_ij weighting and gives
    a numerical exponent <1.5 (originally tracked in issue #44).
    """
    base_layup = [0, 45, -45, 90]
    t_thin = 0.152
    t_thick = 2.0 * t_thin
    lam_thin = _laminate("IM7/8552", base_layup, t_ply=t_thin)
    lam_thick = _laminate("IM7/8552", base_layup, t_ply=t_thick)
    pan, imp = _panel(), _impactor()
    Pc_thin = threshold_load(lam_thin, pan, imp)
    Pc_thick = threshold_load(lam_thick, pan, imp)
    # h_thick / h_thin = 2 -> Pc_thick / Pc_thin = 2^{3/2} = 2.828...
    assert Pc_thick / Pc_thin == pytest.approx(2**1.5, rel=1e-9)
