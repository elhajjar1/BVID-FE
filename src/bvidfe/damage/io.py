"""C-scan / NDE JSON I/O for BVID damage states.

Schema: `docs/cscan_schema.md`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Union

from bvidfe.damage.state import DamageState, DelaminationEllipse

SCHEMA_VERSION = "1.0"


class CScanSchemaError(ValueError):
    """Raised when a C-scan input does not conform to the BVID-FE schema."""


def damage_state_to_dict(ds: DamageState) -> Dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "dent_depth_mm": ds.dent_depth_mm,
        "fiber_break_radius_mm": ds.fiber_break_radius_mm,
        "delaminations": [
            {
                "interface_index": e.interface_index,
                "centroid_mm": list(e.centroid_mm),
                "major_mm": e.major_mm,
                "minor_mm": e.minor_mm,
                "orientation_deg": e.orientation_deg,
            }
            for e in ds.delaminations
        ],
    }


def _validate_dict(data: Dict[str, Any]) -> None:
    if not isinstance(data, dict):
        raise CScanSchemaError("top-level JSON must be an object")
    if "schema_version" not in data:
        raise CScanSchemaError("missing required field: schema_version")
    if data["schema_version"] != SCHEMA_VERSION:
        raise CScanSchemaError(
            f"unsupported schema_version {data['schema_version']!r}; expected {SCHEMA_VERSION!r}"
        )
    if "dent_depth_mm" not in data:
        raise CScanSchemaError("missing required field: dent_depth_mm")
    if data["dent_depth_mm"] < 0:
        raise CScanSchemaError(f"dent_depth_mm must be >= 0 (got {data['dent_depth_mm']})")
    if "delaminations" not in data or not isinstance(data["delaminations"], list):
        raise CScanSchemaError("delaminations must be a list")
    for i, d in enumerate(data["delaminations"]):
        for k in ("interface_index", "centroid_mm", "major_mm", "minor_mm", "orientation_deg"):
            if k not in d:
                raise CScanSchemaError(f"delaminations[{i}] missing field {k!r}")
        if d["major_mm"] <= 0 or d["minor_mm"] <= 0:
            raise CScanSchemaError(
                f"delaminations[{i}] has non-positive axis (major={d['major_mm']}, minor={d['minor_mm']})"
            )
        if d["interface_index"] < 0:
            raise CScanSchemaError(f"delaminations[{i}].interface_index must be >= 0")
        c = d["centroid_mm"]
        if not (isinstance(c, (list, tuple)) and len(c) == 2):
            raise CScanSchemaError(f"delaminations[{i}].centroid_mm must be [x, y]")


def damage_state_from_dict(data: Dict[str, Any]) -> DamageState:
    _validate_dict(data)
    try:
        dels = [
            DelaminationEllipse(
                interface_index=int(d["interface_index"]),
                centroid_mm=(float(d["centroid_mm"][0]), float(d["centroid_mm"][1])),
                major_mm=float(d["major_mm"]),
                minor_mm=float(d["minor_mm"]),
                orientation_deg=float(d["orientation_deg"]),
            )
            for d in data["delaminations"]
        ]
        return DamageState(
            delaminations=dels,
            dent_depth_mm=float(data["dent_depth_mm"]),
            fiber_break_radius_mm=float(data.get("fiber_break_radius_mm", 0.0)),
        )
    except (ValueError, KeyError, TypeError) as exc:
        raise CScanSchemaError(f"invalid damage record: {exc}") from exc


def save_cscan_json(ds: DamageState, path: Union[str, Path]) -> None:
    Path(path).write_text(json.dumps(damage_state_to_dict(ds), indent=2))


def load_cscan_json(path: Union[str, Path]) -> DamageState:
    try:
        data = json.loads(Path(path).read_text())
    except json.JSONDecodeError as exc:
        raise CScanSchemaError(f"invalid JSON: {exc}") from exc
    return damage_state_from_dict(data)
