# BVID-FE Validation

This directory anchors BVID-FE's validation credibility against published
composite-after-impact test data.

## Current state (v0.1.0)

The `datasets/synthetic_selfcheck.json` file contains a tautological
self-check dataset: each case is produced by running the empirical tier
at a known input. By construction the MAE is ~0%. This exists only to
exercise the harness end-to-end and gate regressions in the core analysis
pipeline.

## Roadmap (v0.2.0)

The following published datasets must be digitized by hand (from PDFs)
and added here, each with its own JSON file in `datasets/`. Target MAE%:

| Dataset | Loading | Target MAE | Reference |
|---|---|---|---|
| Soutis AS4/3501-6 | CAI | <12% | Soutis & Curtis (1996) |
| Caprino AS4/epoxy | TAI + CAI | <15% | Caprino (1984) |
| Sanchez-Saez IM7/8552 | CAI | <12% | Sanchez-Saez et al. (2005) |
| NASA round-robin | CAI | <15% | NASA/TM-2007 |

Each case record requires: material name, layup, ply thickness, panel
dimensions, impact energy, measured CAI/TAI strength (MPa), and optionally
the measured dent depth and DPA (for calibration of `olsson_alpha`,
`soutis_k_s`, `dent_beta`).

## Running

```bash
python validation/validate_bvid_public.py              # run all datasets, print table
python validation/validate_bvid_public.py --gate       # exit non-zero if MAE exceeds 1.25 * target
python validation/validate_bvid_public.py --dataset synthetic_selfcheck
```

## Adding a new dataset

1. Create `datasets/<name>.json` following the schema in
   `validate_bvid_public.py` (`DatasetCase` dataclass).
2. The validator auto-discovers any `*.json` file in `datasets/`.
3. Tune material calibration constants in
   `src/bvidfe/core/material.py` until MAE% meets the target.
4. Commit the dataset + any material tuning in a single PR.
