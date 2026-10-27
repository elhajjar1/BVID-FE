import dataclasses
import importlib.util
import statistics
import sys
import warnings
from pathlib import Path

import pytest

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.dent_model import dent_depth_mm, fiber_break_radius_mm
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage
from bvidfe.impact.olsson import onset_energy

ROOT = Path(__file__).resolve().parents[2]


def test_zero_below_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert dent_depth_mm(m, E_impact_J=1.0, E_onset_J=5.0, h_mm=2.0) == 0.0


def test_zero_at_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert dent_depth_mm(m, E_impact_J=5.0, E_onset_J=5.0, h_mm=2.0) == 0.0


def test_monotonic_in_energy():
    m = MATERIAL_LIBRARY["IM7/8552"]
    d1 = dent_depth_mm(m, E_impact_J=6.0, E_onset_J=5.0, h_mm=2.0)
    d2 = dent_depth_mm(m, E_impact_J=15.0, E_onset_J=5.0, h_mm=2.0)
    assert d2 > d1 > 0


@pytest.mark.parametrize("name", sorted(MATERIAL_LIBRARY))
def test_dent_constants_are_the_fitted_defaults(name):
    """One fitted pair for every preset (validation/README.md, dent fit)."""
    m = MATERIAL_LIBRARY[name]
    assert (m.dent_beta, m.dent_gamma) == (0.044, 0.57)


def test_dent_law_is_a_power_of_the_energy_above_onset():
    """d = dent_beta (E - E_onset) ** dent_gamma, independent of the
    thickness and the material's toughness below the h / 2 cap."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    expected = 0.044 * 10.0**0.57
    assert expected == pytest.approx(0.16348, abs=1e-5)
    for h in (2.0, 4.0, 8.0):
        assert dent_depth_mm(m, 15.0, 5.0, h) == pytest.approx(expected, rel=1e-12)
    tougher = dataclasses.replace(m, G_Ic=2 * m.G_Ic)
    assert dent_depth_mm(tougher, 15.0, 5.0, 4.0) == pytest.approx(expected, rel=1e-12)
    custom = dataclasses.replace(m, dent_beta=0.1, dent_gamma=1.0)
    assert dent_depth_mm(custom, 15.0, 5.0, 4.0) == pytest.approx(1.0)


def test_impactor_shape_factors():
    m = MATERIAL_LIBRARY["IM7/8552"]
    hemi = dent_depth_mm(m, 15.0, 5.0, 10.0)
    assert dent_depth_mm(m, 15.0, 5.0, 10.0, shape="hemispherical") == hemi
    assert dent_depth_mm(m, 15.0, 5.0, 10.0, shape="conical") == pytest.approx(8.5 * hemi)
    assert dent_depth_mm(m, 15.0, 5.0, 10.0, shape="flat") == pytest.approx(0.03 * hemi)


def test_impact_mapping_passes_the_impactor_shape():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 6, 0.152)
    pan = PanelGeometry(150, 100)
    dents = {}
    for shape in ("hemispherical", "conical"):
        ev = ImpactEvent(energy_J=30.0, impactor=ImpactorGeometry(16.0, shape=shape))
        dents[shape] = impact_to_damage(ev, lam, pan).dent_depth_mm
    assert 0 < dents["hemispherical"] < dents["conical"]


def test_fiber_break_radius_zero_when_eta_zero():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert m.fiber_break_eta == 0.0
    assert fiber_break_radius_mm(m, E_impact_J=100.0) == 0.0


def test_fiber_break_radius_zero_below_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    m2 = dataclasses.replace(m, fiber_break_eta=0.1, fiber_break_E_threshold=50.0)
    # Below threshold
    assert fiber_break_radius_mm(m2, E_impact_J=10.0) == 0.0
    # Above threshold: positive
    assert fiber_break_radius_mm(m2, E_impact_J=100.0) > 0


def test_dent_depth_capped_at_half_thickness():
    m = MATERIAL_LIBRARY["IM7/8552"]
    h = 1.0  # mm
    # Very high energy that would otherwise produce a huge dent
    d = dent_depth_mm(m, E_impact_J=10_000.0, E_onset_J=1.0, h_mm=h)
    assert d == pytest.approx(0.5 * h)


@pytest.mark.parametrize("n_plies", [8, 16, 24, 32, 40])
@pytest.mark.parametrize("J_per_mm", [2.0, 4.0, 6.7, 10.0])
def test_cap_does_not_bind_for_typical_bvid(n_plies, J_per_mm):
    """Hemispherical 16 mm impact on a 150 x 100 mm QI panel up to 10 J/mm
    (1.5x the ASTM D7136 6.7 J/mm): the dent stays well below h / 2, and in
    the 0.05-0.8 mm range of measured BVID dents. The old law was on the cap
    0.2 J above onset."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [45, 0, -45, 90] * (n_plies // 4), 0.19)
    h = lam.thickness_mm
    ev = ImpactEvent(energy_J=J_per_mm * h, impactor=ImpactorGeometry(15.875))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d = impact_to_damage(ev, lam, PanelGeometry(125, 75)).dent_depth_mm
    if d == 0.0:  # below onset
        return
    assert d < 0.35 * h
    assert 0.02 < d < 0.8


def _validator():
    spec = importlib.util.spec_from_file_location(
        "validate_bvid_public", ROOT / "validation" / "validate_bvid_public.py"
    )
    v = importlib.util.module_from_spec(spec)
    sys.modules["validate_bvid_public"] = v
    try:
        spec.loader.exec_module(v)
    finally:
        sys.modules.pop("validate_bvid_public", None)
    return v


def _predicted_dents(stem):
    """(predicted, measured, h, shape) for each case with a positive measured
    dent, with the onset energy of the impact panel at its centre, as
    ``impact_to_damage`` computes it (the full mapping is ~10x slower)."""
    v = _validator()
    _, cases, _ = v.load_dataset(v.DATASET_DIR / f"{stem}.json")
    out = []
    for c in cases:
        if not c.measured_dent_mm or c.measured_dent_mm <= 0:
            continue
        lam = Laminate(MATERIAL_LIBRARY[c.material], c.layup_deg, c.ply_thickness_mm)
        if c.impact_panel_Lx_mm is not None:
            pan = PanelGeometry(
                c.impact_panel_Lx_mm, c.impact_panel_Ly_mm, c.impact_boundary or c.boundary
            )
        else:
            pan = PanelGeometry(c.panel_Lx_mm, c.panel_Ly_mm, c.boundary)
        imp = ImpactorGeometry(c.impactor_diameter_mm, shape=c.impactor_shape)
        E0 = onset_energy(lam, pan, imp, location_xy_mm=(pan.Lx_mm / 2, pan.Ly_mm / 2))
        d = dent_depth_mm(lam.material, c.impact_energy_J, E0, lam.thickness_mm, imp.shape)
        out.append((d, c.measured_dent_mm, lam.thickness_mm, c.impactor_shape))
    return out


def test_dent_fit_lovejoy_scotti():
    """Median predicted / measured dent on the 8 IM7/8552 coupons (fitted:
    1.00; 8.2 with the old law). The band leaves room for onset-energy
    changes: the law is fitted on E - E_onset, and halving the onset moves
    this median by under 4%."""
    rows = _predicted_dents("lovejoy_scotti_im7_8552_cai")
    assert len(rows) == 8
    ratios = [p / m for p, m, _, _ in rows]
    assert statistics.median(ratios) == pytest.approx(1.0, abs=0.08)
    assert max(ratios) / min(ratios) < 1.3
    assert all(p < 0.5 * h for p, _, h, _ in rows)


def test_dent_fit_hasebe():
    """Median predicted / measured dent over the 308 Hasebe et al. specimens
    with a positive measured dent (fitted: 0.96 overall, 0.96 hemispherical,
    0.98 conical; 9.0 with the old law); 0.81-1.09 for onset energies 0.5-1.5x
    the current one."""
    rows = _predicted_dents("hasebe_t800s_3900_2b_cai")
    assert len(rows) == 308
    ratios = [p / m for p, m, _, _ in rows]
    assert 0.8 < statistics.median(ratios) < 1.15
    for shape in ("hemispherical", "conical"):
        r = [p / m for p, m, _, s in rows if s == shape]
        assert 0.8 < statistics.median(r) < 1.2, shape
    # The h / 2 cap binds on none of them.
    assert all(p < 0.5 * h for p, _, h, _ in rows)
