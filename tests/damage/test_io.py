import json

import pytest

from bvidfe.damage.io import (
    CScanSchemaError,
    damage_state_from_dict,
    damage_state_to_dict,
    load_cscan_json,
    save_cscan_json,
)
from bvidfe.damage.state import DamageState, DelaminationEllipse


def _make_state() -> DamageState:
    return DamageState(
        [
            DelaminationEllipse(3, (75, 50), 28, 18, 45),
            DelaminationEllipse(4, (78, 52), 32, 20, 50),
        ],
        dent_depth_mm=0.45,
        fiber_break_radius_mm=3.0,
    )


def test_round_trip_to_dict_and_back():
    ds = _make_state()
    ds2 = damage_state_from_dict(damage_state_to_dict(ds))
    assert ds2.dent_depth_mm == ds.dent_depth_mm
    assert ds2.fiber_break_radius_mm == ds.fiber_break_radius_mm
    assert len(ds2.delaminations) == 2
    assert ds2.delaminations[0].interface_index == 3
    assert ds2.delaminations[1].major_mm == 32


def test_round_trip_file(tmp_path):
    ds = _make_state()
    fp = tmp_path / "cscan.json"
    save_cscan_json(ds, fp)
    ds2 = load_cscan_json(fp)
    assert ds2.delaminations[1].orientation_deg == 50


def test_rejects_bad_schema_version(tmp_path):
    fp = tmp_path / "bad.json"
    fp.write_text(json.dumps({"schema_version": "99.0", "delaminations": [], "dent_depth_mm": 0}))
    with pytest.raises(CScanSchemaError):
        load_cscan_json(fp)


def test_rejects_missing_schema_version(tmp_path):
    fp = tmp_path / "bad.json"
    fp.write_text(json.dumps({"delaminations": [], "dent_depth_mm": 0}))
    with pytest.raises(CScanSchemaError):
        load_cscan_json(fp)


def test_rejects_negative_ellipse(tmp_path):
    fp = tmp_path / "bad.json"
    fp.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "dent_depth_mm": 0.0,
                "delaminations": [
                    {
                        "interface_index": 0,
                        "centroid_mm": [0, 0],
                        "major_mm": -1,
                        "minor_mm": 5,
                        "orientation_deg": 0,
                    }
                ],
            }
        )
    )
    with pytest.raises(CScanSchemaError):
        load_cscan_json(fp)


def test_rejects_negative_dent(tmp_path):
    fp = tmp_path / "bad.json"
    fp.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "dent_depth_mm": -0.1,
                "delaminations": [],
            }
        )
    )
    with pytest.raises(CScanSchemaError):
        load_cscan_json(fp)
