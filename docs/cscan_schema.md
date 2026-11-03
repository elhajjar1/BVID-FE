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
- `major_mm` (number, > 0): semi-major axis length in millimeters (half the delamination's length along that axis, not the full length). The ellipse area is `π · major_mm · minor_mm` mm².
- `minor_mm` (number, > 0): semi-minor axis length in millimeters. Nothing requires `minor_mm <= major_mm`; the names only label the axis that `orientation_deg` rotates.
- `orientation_deg` (number): angle in degrees from the panel +x axis to the `major_mm` axis, positive from +x towards +y. This is the same sense as the ply angles in `layup_deg` (a fibre at +45° and an ellipse at +45° are parallel). Any real value is accepted; 0° and 180° describe the same ellipse.

### Which face is which

The file has no field for the impacted face; the convention comes from
the layup order. Ply 0 is the first entry of `layup_deg`. It lies at
`z = −h/2` in the classical-lamination stack and at `z = 0` in the fe3d
mesh. The impact
model (`impact_to_damage`) treats ply 0 as the impact face: its
delaminations grow from interface 0 (smallest) to interface `n − 2`
(largest, back face). So list `layup_deg` starting from the impacted
face, and number the C-scan interfaces from that face: interface `i`
lies between the `i`-th and `(i+1)`-th plies counted from the impact
side, both counted from 0. The semi-analytical growth check (which takes
the thinner of plies `0 … i` and `i+1 … n−1`) and fe3d both read
`interface_index` against `layup_deg`, so the file and the layup must
count from the same face.

### Units

Lengths in millimeters throughout (`dent_depth_mm`, `fiber_break_radius_mm`,
`centroid_mm`, `major_mm`, `minor_mm`); angles in degrees. Inches or
in² from an inspection report must be converted (1 in = 25.4 mm).

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
  the ellipse with the largest area: its sublaminate buckles over the
  ellipse's bounding box in the panel frame, `2 · ex` along the load
  direction (x) by `2 · ey`, where `(ex, ey)` are the rotated ellipse's
  half-extents (`DelaminationEllipse.half_extents_mm()`; `major_mm` by
  `minor_mm` at `orientation_deg = 0`, swapped at 90°). The TAI hole
  width reads the same `ey`.
- fe3d places the damage on its mesh by the ellipses' positions, sizes and
  orientations, and softens the elements within `fiber_break_radius_mm` of
  any delamination centroid.
- `dent_depth_mm` is stored and reported but does not enter any
  residual-strength calculation.

The file does not carry the layup, so `interface_index` is checked when
the damage is given to `AnalysisConfig`: an index past the last interface
(`n_plies − 2`) raises `ValueError`.

## Validation errors

Files that fail validation raise `bvidfe.damage.io.CScanSchemaError` (a
subclass of `ValueError`) with a message indicating the offending field.
Malformed JSON raises the same error.
