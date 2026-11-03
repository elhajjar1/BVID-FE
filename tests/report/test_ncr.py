"""Unit tests for ``bvidfe.report.ncr`` (NCR / MRB disposition, #172)."""

from __future__ import annotations

import pytest

import bvidfe.report.ncr as ncr_module
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.report.ncr import (
    BVID_DENT_THRESHOLD_MM,
    DISPOSITION_STATUS,
    KD_REPAIR,
    KD_SUBSTANTIATE,
    KD_USE_AS_IS,
    build_ncr,
    recommend_disposition,
)

UAI = "Use-As-Is (UAI)"
SUBSTANTIATE = "Use-As-Is (UAI) with engineering substantiation, or Repair"
REPAIR = "Repair"
SCRAP = "Repair (structural) or Scrap"


def _result(kd=0.9, dent=0.1, fbr=0.0, n_delam=0, dpa=500.0, tier="empirical") -> dict:
    """A minimal ``AnalysisResults.to_dict()``-shaped result."""
    ellipse = {
        "interface_index": 0,
        "centroid_mm": [75.0, 50.0],
        "major_mm": 10.0,
        "minor_mm": 5.0,
        "orientation_deg": 0.0,
    }
    return {
        "knockdown": kd,
        "residual_strength_MPa": 600.0 * kd,
        "pristine_strength_MPa": 600.0,
        "dpa_mm2": dpa,
        "tier_used": tier,
        "damage": {
            "dent_depth_mm": dent,
            "fiber_break_radius_mm": fbr,
            "delaminations": [ellipse] * n_delam,
        },
        "notes": [],
        "warnings": [],
    }


def _cfg(ply_thickness_mm=0.152) -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=ply_thickness_mm,
        panel=PanelGeometry(150.0, 100.0, boundary="clamped"),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(20.0),
    )


FORM = {
    "ncr_number": "NCR-20261103-042",
    "revision": "B",
    "originator": "A. Engineer",
    "date": "2026-11-03",
    "program": "Wing-Box-X",
    "part_number": "PN-12345-01",
    "part_name": "Lower skin panel",
    "serial_lot": "SN-0099",
    "work_order": "WO-777",
    "quantity": "1",
    "station": "Bay 4",
    "found_during": "Final inspection",
    "nc_type": "Impact damage (BVID/VID)",
    "dimensions": "25 mm dia. dent, 0.4 mm deep",
    "description": "Tool drop found on the outer mould line.",
}


def test_module_states_bands_are_placeholders() -> None:
    doc = ncr_module.__doc__.lower()
    assert "placeholder" in doc and "not taken from a cited standard" in doc
    assert (KD_USE_AS_IS, KD_SUBSTANTIATE, KD_REPAIR) == (0.97, 0.85, 0.70)
    assert BVID_DENT_THRESHOLD_MM == 0.30


@pytest.mark.parametrize("dent, flagged", [(0.30, True), (0.45, True), (0.2999, False)])
def test_dent_at_or_above_threshold_is_flagged_visible(dent, flagged) -> None:
    rec = recommend_disposition(_result(kd=0.99, dent=dent), dent_threshold_mm=0.30)
    assert any("VID, not BVID" in f for f in rec["flags"]) is flagged
    # A visible dent rules out a plain Use-As-Is even at a high knockdown.
    assert rec["recommended_path"] == (SUBSTANTIATE if flagged else UAI)


def test_dent_threshold_is_a_parameter() -> None:
    rec = recommend_disposition(_result(kd=0.99, dent=0.4), dent_threshold_mm=0.5)
    assert not any("VID, not BVID" in f for f in rec["flags"])
    assert rec["recommended_path"] == UAI
    assert any("0.500 mm design threshold" in c for c in rec["criteria"])


@pytest.mark.parametrize("kd", [1.0, 0.9, 0.5])
def test_fibre_break_forces_repair(kd) -> None:
    rec = recommend_disposition(_result(kd=kd, fbr=1.5))
    assert rec["recommended_path"] == REPAIR
    assert any("Fiber breakage present" in f for f in rec["flags"])
    assert "fiber breakage present" in rec["rationale"]


@pytest.mark.parametrize(
    "kd, expected",
    [
        (1.0, UAI),
        (KD_USE_AS_IS, UAI),
        (KD_USE_AS_IS - 1e-6, SUBSTANTIATE),
        (KD_SUBSTANTIATE, SUBSTANTIATE),
        (KD_SUBSTANTIATE - 1e-6, REPAIR),
        (KD_REPAIR, REPAIR),
        (KD_REPAIR - 1e-6, SCRAP),
        (0.3, SCRAP),
    ],
)
def test_knockdown_band_boundaries(kd, expected) -> None:
    rec = recommend_disposition(_result(kd=kd))
    assert rec["recommended_path"] == expected
    assert rec["status"] == DISPOSITION_STATUS
    assert rec["knockdown"] == kd


def test_delamination_rules_out_plain_use_as_is() -> None:
    rec = recommend_disposition(_result(kd=0.99, n_delam=2))
    assert rec["recommended_path"] == SUBSTANTIATE
    assert any("2 delamination(s) modelled" in f for f in rec["flags"])


def test_bands_follow_the_module_constants(monkeypatch) -> None:
    monkeypatch.setattr(ncr_module, "KD_REPAIR", 0.75)
    assert ncr_module.recommend_disposition(_result(kd=0.72))["recommended_path"] == SCRAP


def test_markdown_and_dict_contain_their_inputs() -> None:
    result = _result(kd=0.812, dent=0.35, n_delam=3, dpa=1234.0, tier="semi_analytical")
    result["notes"] = ["closed-form TAI note"]
    result["warnings"] = ["dpa_capped"]
    rec = recommend_disposition(result)
    md, d = build_ncr(FORM, _cfg(), result, rec, generated_utc="2026-11-03T12:00:00+00:00")

    for value in FORM.values():
        assert value in md, value
    for text in (
        "IM7/8552",
        "[0/45/-45/90/90/-45/45/0] (8 plies)",
        "0.152 mm",
        "150 × 100 mm, clamped",
        "semi_analytical",
        "600.0 MPa",
        "487.2 MPa",
        "0.812",
        "1234 mm²",
        "0.350 mm",
        "Delaminations modelled:** 3",
        "closed-form TAI note",
        "dpa_capped",
        rec["recommended_path"],
        rec["status"],
        "2026-11-03T12:00:00+00:00",
    ):
        assert text in md, text
    for flag in rec["flags"]:
        assert flag in md
    for criterion in rec["criteria"]:
        assert criterion in md

    assert d["ncr_number"] == FORM["ncr_number"]
    assert d["serial_or_lot"] == FORM["serial_lot"]
    assert d["material_configuration"]["panel_boundary"] == "clamped"
    assert d["engineering_analysis"]["n_delaminations"] == 3
    assert d["engineering_analysis"]["projected_damage_area_mm2"] == 1234.0
    assert d["proposed_disposition"] is rec
    assert d["generated_utc"] == "2026-11-03T12:00:00+00:00"


def test_build_ncr_without_analysis_marks_it_pending() -> None:
    md, d = build_ncr(FORM, _cfg(), None, None)
    assert "Engineering analysis pending" in md
    assert "No disposition recommendation" in md
    assert d["engineering_analysis"] is None and d["proposed_disposition"] is None
    assert d["generated_utc"]  # defaults to now


def test_build_ncr_is_reproducible_with_a_timestamp() -> None:
    args = (FORM, _cfg(), _result(), recommend_disposition(_result()))
    assert build_ncr(*args, generated_utc="t") == build_ncr(*args, generated_utc="t")


def test_build_ncr_per_ply_thickness() -> None:
    md, _ = build_ncr(FORM, _cfg([0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1]), None, None)
    assert "0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1 mm (per ply)" in md


def test_recommendation_from_a_real_result() -> None:
    cfg = _cfg()
    result = BvidAnalysis(cfg).run().to_dict()
    rec = recommend_disposition(result)
    assert rec["recommended_path"] in {UAI, SUBSTANTIATE, REPAIR, SCRAP}
    md, _ = build_ncr(FORM, cfg, result, rec)
    assert f"{result['knockdown']:.3f}" in md
