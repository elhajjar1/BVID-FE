"""Damage area by energy balance: DPA = alpha (E - E_onset) / G_IIc.

It used to be divided by the laminate thickness too, which left DPA in mm
and flat in thickness at constant energy per thickness. Measured areas grow
with thickness there (Hasebe et al. 2025; Lovejoy & Scotti 2019).
"""

from __future__ import annotations

import pytest

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage
from bvidfe.impact.olsson import onset_energy

_M = MATERIAL_LIBRARY["IM7/8552"]
_PANEL = PanelGeometry(150.0, 100.0)
_IMP = ImpactorGeometry(16.0)


def _dpa(n_rep, energy):
    lam = Laminate(_M, [45, 0, -45, 90] * n_rep + [90, -45, 0, 45] * n_rep, 0.19)
    ds = impact_to_damage(ImpactEvent(energy, _IMP, 5.5), lam, _PANEL)
    return lam, ds.projected_damage_area_mm2


def test_dpa_is_the_energy_balance_at_the_astm_reference():
    lam, dpa = _dpa(4, 54.0)
    e_onset = onset_energy(lam, _PANEL, _IMP)
    expected = _M.olsson_alpha * (54.0 - e_onset) * 1e3 / _M.G_IIc
    assert dpa == pytest.approx(expected, rel=0.01)


def test_dpa_grows_with_thickness_at_constant_energy_per_thickness():
    """Lovejoy & Scotti: 1.25x the thickness and energy at 8.9 J/mm gave
    1.21x the C-scan area."""
    lam_a, dpa_a = _dpa(4, 54.0)
    lam_b, dpa_b = _dpa(5, 54.0 * 5 / 4)
    assert lam_b.thickness_mm / lam_a.thickness_mm == pytest.approx(1.25)
    assert 1.15 < dpa_b / dpa_a < 1.3


def test_default_alpha_matches_lovejoy_scotti_c_scans():
    """Median measured / predicted DPA on the eight IM7/8552 coupons."""
    import importlib.util
    import statistics
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location(
        "validate_bvid_public", root / "validation" / "validate_bvid_public.py"
    )
    v = importlib.util.module_from_spec(spec)
    sys.modules["validate_bvid_public"] = v
    try:
        spec.loader.exec_module(v)
    finally:
        sys.modules.pop("validate_bvid_public", None)
    _, cases, _ = v.load_dataset(v.DATASET_DIR / "lovejoy_scotti_im7_8552_cai.json")
    ratios = [c.measured_dpa_mm2 / v.run_case(c, "empirical")["predicted_dpa_mm2"] for c in cases]
    assert statistics.median(ratios) == pytest.approx(1.0, abs=0.03)
