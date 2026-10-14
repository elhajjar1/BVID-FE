# C-scan / NDE JSON schema (v1.0)

BVID-FE reads inspection data describing damage in a panel using the JSON format below.
`bvidfe.damage.io.load_cscan_json(path)` returns a `DamageState`, and
`save_cscan_json(damage_state, path)` writes one in the same format (the
CLI reads it with `--cscan path.json`). A complete example ships as
`examples/sample_cscan.json`.

## Example

```json
{
  "schema_version": "1.0",
  "dent_depth_mm": 0.45,
  "fiber_break_radius_mm": 3.0,
  "delaminations": [
    {
      "interface_index": 3,
      "centroid_mm": [75.0, 50.0],
      "major_mm": 28.0,
      "minor_mm": 18.0,
      "orientation_deg": 45.0
    }
  ]
}
```

## Fields

- `schema_version` (string, required): must be `"1.0"`.
- `dent_depth_mm` (number, required): residual dent depth on the impact face, millimeters. Must be >= 0.
- `fiber_break_radius_mm` (number, optional, default 0.0): radius of the fiber-break core in millimeters. Must be >= 0.
- `delaminations` (array, required): one object per observed delamination. Must be a list (can be empty).

Each delamination object (all five fields required):
- `interface_index` (integer, >= 0): ply interface index (`0` = between plies 0 and 1, so the last interface of an `n`-ply layup is `n - 2`). Plies are numbered in `layup_deg` order.
- `centroid_mm` (array `[x, y]`): ellipse center in the panel frame, millimeters. The panel spans `0 <= x <= Lx_mm` (x is the load direction) and `0 <= y <= Ly_mm`.
- `major_mm` (number, > 0): semi-major axis length in millimeters.
- `minor_mm` (number, > 0): semi-minor axis length in millimeters.
- `orientation_deg` (number): rotation of the `major_mm` axis from the panel x axis, counterclockwise, in degrees.

Every number must be finite (`NaN` and `Infinity` are rejected), and
`interface_index` must be a JSON integer (`3`, not `3.0` or `true`).
Unknown fields, at the top level or in a delamination, are ignored with a
`UserWarning`, so a typo such as `"dent_dept_mm"` is reported rather than
silently dropped.

## How the tiers use it

- **Projected damage area (DPA)** is the area of the union of all ellipse
  footprints (`DamageState.projected_damage_area_mm2`), so overlapping
  delaminations at different interfaces count once. The empirical tier
  uses only this area.
- The semi-analytical and fe3d CAI growth check uses, at each interface,
  the ellipse with the largest area: its sublaminate buckles over a
  rectangle of `2 · major_mm` along the load direction (x) by
  `2 · minor_mm`. `orientation_deg` does not enter this check.
- fe3d places the damage on its mesh by the ellipses' positions, sizes and
  orientations, and softens the elements within `fiber_break_radius_mm` of
  any delamination centroid.
- `dent_depth_mm` is stored and reported but does not enter any
  residual-strength calculation.

The file does not carry the layup, so `interface_index` is not checked
against it when loading. An index past the last interface gives a
degenerate sublaminate (no growth bound) in the closed-form tiers and can
make fe3d raise an `IndexError`.

## Validation errors

Files that fail validation raise `bvidfe.damage.io.CScanSchemaError` (a
subclass of `ValueError`) with a message indicating the offending field.
Malformed JSON raises the same error.
