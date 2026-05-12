"""Olsson threshold-load + onset-energy validation.

The Olsson (2001) quasi-static damage-threshold model predicts
   Pc = pi * sqrt(8 * G_IIc * D_eff / 9)
where D_eff is the geometric-mean flexural rigidity. These tests pin the
shape of the prediction (toughness scaling, thickness scaling, panel
size invariance) rather than absolute literature values — the absolute
calibration would need digitised Olsson plots, which is out of scope
for a self-contained regression test.

The shape relationships, in contrast, follow directly from the closed
form and must hold for every BVID-FE material:
  * Pc scales as sqrt(G_IIc) when D_eff is held fixed.
  * Pc scales as h^{3/2} when E_ij are held fixed AND the stack is thick
    enough that the discrete-CLT z-integration approaches the homogenised
    continuum (D_eff ~ E*h^3). For a short asymmetric stack like
    [0,45,-45,90] the 4->8 ply repeat is not self-similar (the new ply
    at z=+/-h/4 lands on a 0deg in the doubled stack, not the 90deg of
    the original); D_eff there only scales as ~h^{2.81}. The asymptote
    is recovered to within 0.1% by 32->64 plies of a quasi-iso
    [0,45,-45,90]_s sub-laminate. See issue #44.
  * Pc is a property of the laminate + impactor and is INDEPENDENT of
    panel size, so the same Laminate+ImpactorGeometry should yield the
    same Pc on a 100x100 mm and 200x150 mm panel.
"""

from __future__ import annotations

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


def test_threshold_load_scales_as_sqrt_of_giic():
    """Pc ~ sqrt(G_IIc) when D_eff is fixed.

    Hold the laminate stack fixed (same E_ij, same ply count) and swap
    materials with different G_IIc only. We use IM7/8552 (G_IIc=0.79)
    and T700/2510 (G_IIc=0.60) and check the ratio matches sqrt of the
    G_IIc ratio. The two materials have different E_ij too, so we
    isolate the G_IIc effect by comparing Pc / sqrt(D_eff)."""
    lam_a = _laminate("IM7/8552", [0, 45, -45, 90, 90, -45, 45, 0])
    lam_b = _laminate("T700/2510", [0, 45, -45, 90, 90, -45, 45, 0])
    pan, imp = _panel(), _impactor()
    Pc_a = threshold_load(lam_a, pan, imp)
    Pc_b = threshold_load(lam_b, pan, imp)
    # Pc / sqrt(D_eff * G_IIc) must be the same constant for both materials
    # (it equals pi * sqrt(8/9)).
    expected_const = math.pi * math.sqrt(8.0 / 9.0)
    Deff_a = lam_a.flexural_rigidity_Deff()
    Deff_b = lam_b.flexural_rigidity_Deff()
    actual_a = Pc_a / math.sqrt(Deff_a * lam_a.material.G_IIc)
    actual_b = Pc_b / math.sqrt(Deff_b * lam_b.material.G_IIc)
    assert actual_a == pytest.approx(expected_const, rel=1e-9)
    assert actual_b == pytest.approx(expected_const, rel=1e-9)


def test_threshold_load_invariant_to_panel_size():
    """Pc depends on the LAMINATE, not the panel (per the Olsson model)."""
    lam = _laminate("IM7/8552", [0, 45, -45, 90, 90, -45, 45, 0])
    imp = _impactor()
    Pc_small = threshold_load(lam, _panel(100, 80), imp)
    Pc_large = threshold_load(lam, _panel(300, 200), imp)
    assert Pc_small == pytest.approx(Pc_large, rel=1e-12)


def test_threshold_load_scales_as_h_three_halves_for_homogenised_layup():
    """Pc ~ h^{3/2} holds when the stack is repeated enough times that the
    discrete ply z-distribution approaches a homogenised section.

    Issue #44: doubling a 4-ply [0,45,-45,90] stack is NOT self-similar -
    the new ply at z=+/-h/4 in the doubled stack is a 0deg, not the 90deg
    of the original. D_eff there only scales as ~h^{2.81}, so the Pc
    exponent comes out at ~1.40 (observed 2.647 vs theoretical 2.828).
    For [0,45,-45,90]_s repeated 4x and 8x (32 vs 64 plies) the ply z
    distribution IS self-similar under doubling and D_eff scales as h^3
    to within 0.1%; the Pc ratio is 2.8262 (exponent 1.4989, within
    rel=2e-3 of 2^{3/2}).
    """
    sub = [0, 45, -45, 90, 90, -45, 45, 0]  # quasi-iso symmetric sub-laminate
    lam_thin = _laminate("IM7/8552", sub * 4)  # 32 plies, ~4.86 mm
    lam_thick = _laminate("IM7/8552", sub * 8)  # 64 plies, ~9.73 mm
    pan, imp = _panel(), _impactor()
    Pc_thin = threshold_load(lam_thin, pan, imp)
    Pc_thick = threshold_load(lam_thick, pan, imp)
    assert Pc_thick / Pc_thin == pytest.approx(2**1.5, rel=2e-3)


def test_threshold_load_matches_precomputed_reference_values():
    """Numerical-drift sentinel for the full Pc pipeline.

    Hard-codes Pc for two IM7/8552 quasi-iso laminates computed from the
    current verified implementation. Any change to the CLT D-matrix
    assembly, Q-bar transforms, G_IIc plumbing, or the Olsson formula
    prefactor will trip these. The ratio-based test above can mask a
    constant multiplicative bug (it cancels in Pc_thick / Pc_thin); this
    one cannot.
    """
    pan, imp = _panel(), _impactor()
    lam8 = _laminate("IM7/8552", [0, 45, -45, 90, 90, -45, 45, 0])
    lam16 = _laminate("IM7/8552", [0, 45, -45, 90, 90, -45, 45, 0] * 2)
    assert threshold_load(lam8, pan, imp) == pytest.approx(243.7361324522721, rel=1e-9)
    assert threshold_load(lam16, pan, imp) == pytest.approx(756.5779536988468, rel=1e-9)
