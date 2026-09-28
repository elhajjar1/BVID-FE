"""Whitney-Nuismer TAI analytical-limit validation.

The point-stress TAI formula in ``failure.soutis_openhole.whitney_nuismer_tai``
has two analytical limits that must hold for ANY material and any DPA:

  (a) DPA -> 0 : knockdown -> 1.0 (no hole, no strength loss).
  (b) DPA -> inf : xi -> 1, denom = 2 + 1 + 3 - (Kt_inf - 3)*(5 - 7) = 2 * Kt_inf;
      so knockdown -> 2 / (2 * Kt_inf) = 1 / Kt_inf. With the default
      Kt_inf = 3 the asymptote is exactly 1/3.

  (c) Monotonicity: knockdown decreases as DPA increases (the hole-effect
      stress concentration intensifies with hole size).

These three limits exercise pure algebra in the formula and should hold
to machine precision once xi is sufficiently far from its boundaries.
"""

from __future__ import annotations

import math
import warnings

import pytest

from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.soutis_openhole import lekhnitskii_kt_infinity, whitney_nuismer_tai


def test_wn_returns_pristine_at_zero_dpa():
    """The formula short-circuits to sigma_pristine at dpa <= 0."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert whitney_nuismer_tai(m, dpa_mm2=0.0, sigma_pristine_MPa=600.0) == 600.0


def test_wn_knockdown_approaches_unity_for_tiny_dpa():
    """As DPA -> 0+, knockdown should approach (but not exactly equal) 1.0."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    sigma_0 = 600.0
    # Pick DPA so xi = R/(R + d0) is small (<< 1)
    sigma = whitney_nuismer_tai(m, dpa_mm2=1e-4, sigma_pristine_MPa=sigma_0)
    assert sigma == pytest.approx(sigma_0, rel=1e-3)
    assert sigma <= sigma_0  # never exceeds pristine


def test_wn_knockdown_asymptote_at_infinite_dpa_equals_one_third():
    """For Kt_inf=3, the formula's asymptote is 1/3 as DPA -> inf (xi -> 1)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    sigma_0 = 600.0
    sigma = whitney_nuismer_tai(m, dpa_mm2=1e10, sigma_pristine_MPa=sigma_0, Kt_inf=3.0)
    expected_asymptote = sigma_0 / 3.0
    assert sigma == pytest.approx(expected_asymptote, rel=1e-4)


def test_wn_asymptote_scales_with_inverse_Kt_inf():
    """For arbitrary Kt_inf the asymptote is 1 / Kt_inf (algebraic, see
    module docstring). Verify with Kt_inf = 2 and Kt_inf = 5."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    sigma_0 = 600.0
    sigma2 = whitney_nuismer_tai(m, dpa_mm2=1e10, sigma_pristine_MPa=sigma_0, Kt_inf=2.0)
    sigma5 = whitney_nuismer_tai(m, dpa_mm2=1e10, sigma_pristine_MPa=sigma_0, Kt_inf=5.0)
    assert sigma2 == pytest.approx(sigma_0 / 2.0, rel=1e-4)
    # Convergence to the algebraic asymptote is slower at higher Kt_inf;
    # at dpa=1e10 the Kt_inf=5 residual is ~1.2e-4 relative.
    assert sigma5 == pytest.approx(sigma_0 / 5.0, rel=5e-4)


def test_wn_monotonically_decreases_with_dpa():
    """sigma_TAI must decrease as DPA grows."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    sigma_0 = 600.0
    series = [
        whitney_nuismer_tai(m, dpa_mm2=dpa, sigma_pristine_MPa=sigma_0)
        for dpa in (0.1, 1.0, 10.0, 100.0, 1000.0, 10_000.0)
    ]
    for prev, nxt in zip(series, series[1:]):
        assert prev > nxt, f"non-monotone: {series}"


def test_wn_knockdown_bounded_below_by_asymptote():
    """For Kt_inf=3 no finite DPA should yield knockdown < 1/3 (the
    asymptote is the lower bound)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    sigma_0 = 600.0
    for dpa in (1, 10, 100, 1000, 10_000, 1e6):
        sigma = whitney_nuismer_tai(m, dpa_mm2=dpa, sigma_pristine_MPa=sigma_0, Kt_inf=3.0)
        # Allow a tiny slack for floating-point drift near the asymptote.
        assert sigma >= sigma_0 / 3.0 - 1e-9, f"undershoot at DPA={dpa}: {sigma}"
        assert math.isfinite(sigma)


@pytest.mark.parametrize(
    "Kt_inf", [1.2, 2.0, 3.0, 5.0, 7.3, 9.0, 9.5, 12.0, 20.0, 20.4, 30.0, 60.0]
)
def test_wn_knockdown_bounded_in_unit_interval(Kt_inf):
    """0 < knockdown <= 1 for every Kt_inf and DPA, including past the
    Kt_inf ~ 9.2 onset of kd > 1 and the Kt_inf ~ 20.3 denominator pole."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    sigma_0 = 600.0
    for dpa in [10.0**e for e in range(-3, 9)] + [5.0, 20.0, 50.0, 200.0]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            sigma = whitney_nuismer_tai(m, dpa_mm2=dpa, sigma_pristine_MPa=sigma_0, Kt_inf=Kt_inf)
        assert math.isfinite(sigma)
        assert 0.0 < sigma <= sigma_0, f"Kt_inf={Kt_inf}, DPA={dpa}: {sigma}"


@pytest.mark.parametrize("Kt_inf", [12.0, 30.0])
def test_wn_warns_and_clamps_outside_valid_range(Kt_inf):
    """xi ~ 0.69 is where the stress-ratio polynomial dips lowest; at
    Kt_inf=12 the raw knockdown there is ~1.3 and at Kt_inf=30 the
    denominator is negative."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    d0 = m.wn_d0_mm
    R = 0.69 / (1 - 0.69) * d0  # xi = R / (R + d0) = 0.69
    dpa = math.pi * R**2
    with pytest.warns(UserWarning, match="outside its valid range"):
        sigma = whitney_nuismer_tai(m, dpa_mm2=dpa, sigma_pristine_MPa=600.0, Kt_inf=Kt_inf)
    assert sigma == 600.0


@pytest.mark.parametrize("name", sorted(MATERIAL_LIBRARY))
@pytest.mark.parametrize(
    "layup", [[0] * 8, [0, 90] * 4, [0, 45, -45, 90, 90, -45, 45, 0], [45, -45] * 4]
)
def test_wn_presets_never_hit_the_clamp(name, layup):
    """The clamp is a guard for custom cards; preset laminates stay in range."""
    m = MATERIAL_LIBRARY[name]
    Kt_inf = lekhnitskii_kt_infinity(Laminate(m, layup, 0.15))
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        for dpa in (1.0, 10.0, 100.0, 1000.0, 10_000.0):
            whitney_nuismer_tai(m, dpa_mm2=dpa, sigma_pristine_MPa=600.0, Kt_inf=Kt_inf)
