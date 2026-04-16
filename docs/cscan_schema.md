# C-scan / NDE JSON schema (v1.0)

BVID-FE reads inspection data describing damage in a panel using the JSON format below.

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
- `fiber_break_radius_mm` (number, optional, default 0.0): radius of the fiber-break core in millimeters.
- `delaminations` (array, required): one object per observed delamination. Must be a list (can be empty).

Each delamination object:
- `interface_index` (integer, >= 0): ply interface index (`0` = between plies 0 and 1).
- `centroid_mm` (array `[x, y]`): ellipse center in the panel frame, millimeters.
- `major_mm` (number, > 0): semi-major axis length in millimeters.
- `minor_mm` (number, > 0): semi-minor axis length in millimeters.
- `orientation_deg` (number): ellipse major-axis rotation in degrees, measured in the panel XY plane.

## Validation errors

Files that fail validation raise `bvidfe.damage.io.CScanSchemaError` with a message indicating the offending field.
