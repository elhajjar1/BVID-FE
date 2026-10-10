# BVID-FE Validation

This directory anchors BVID-FE's validation credibility against published
composite-after-impact test data.

## Datasets

| Dataset | Material | Loading | Cases | Target MAE | Empirical | Semi-analytical |
|---|---|---|---|---|---|---|
| `ncamp_as4_8552_cai` | AS4/8552 | CAI | 7 | 15% | 15.9% | 11.6% |
| `lovejoy_scotti_im7_8552_cai` (advisory) | IM7/8552 | CAI | 8 | 15% | 39.9% | 23.5% |
| `synthetic_selfcheck` | IM7/8552 | CAI + TAI | 12 | 1% | ~0% (by construction) | n/a |

**`ncamp_as4_8552_cai`** is measured data: NCAMP CAM-RP-2010-002 Rev A
(2011), the Hexcel 8552 AS4 unitape qualification report, sections 2.3.31
and 4.31. It has seven [45/0/-45/90]3s coupons (24 plies, 4.5 mm), tested
to ASTM D7136/D7137 at RTD and impacted at 1500 in-lbf/in (about 30 J). The
measured CAI is 171–182 MPa and every coupon failed through the impact damage
(failure code LDM). Each case uses its own measured thickness and impact
energy. The report does not give the impactor mass or tup, so the ASTM
D7136 standard values (5.5 kg, 0.625 in hemispherical) are used. It reports
no dent depth or damage area. The fe3d tier gives 169 MPa (5% low) on the
first coupon, but at about 3.5 minutes per coupon it is not run in CI.

**`lovejoy_scotti_im7_8552_cai`** is measured data: Lovejoy & Scotti, NASA
LaRC 2019 (NTRS 20200002432), Tables 1 and 4. It has the eight standard-ply
(190 gsm) coupons, 4 x 6 in, impacted and tested at NIAR to ASTM D7136/D7137:
- four quasi-isotropic [45/0/-45/90]4s coupons (32 plies, 6.1 mm) at 40 ft-lb
  (54 J), 133–147 MPa;
- four hard [45/0/90/0/-45/0/45/0/-45/0]2s coupons (40 plies, 7.6 mm) at
  50 ft-lb (68 J), 236–249 MPa.

Each case carries the measured thickness, dent depth and C-scan damage area.
The impactor is not described, so the ASTM D7136 standard is assumed. The
hybrid thin-ply coupons are omitted, because a case has a single ply
thickness.

This dataset is **advisory** (`"gate": false`): the model misses it.
- The quasi-isotropic coupons are predicted 36–56% too strong.
- The damage area is overpredicted: about 8,900 mm² against 5,200–7,500 mm²
  measured.
- Rescaling the predicted delaminations to the measured areas makes the
  strengths higher still.
- The undamaged strength is right (601 MPa against 600 MPa measured by
  NCAMP).

So the damage-to-knockdown map is too lenient: `soutis_k_s` = 2.5 is an
uncalibrated default, and these coupons imply about 4–5.

**`synthetic_selfcheck`** is tautological: each case is produced by running
the empirical tier at a known input, so its MAE is ~0% by construction. It
exercises the harness end-to-end and gates regressions in the empirical
pipeline. Regenerate it when the empirical tier changes on purpose.

CI gates the empirical tier on every dataset and the semi-analytical tier on
`ncamp_as4_8552_cai`, at 1.25x each dataset's target. A dataset with
`"gate": false` is still run and its error printed, with its `gate_note`, but
it cannot fail the gate. The fe3d tier is run on the
synthetic set as an advisory check.

## Next datasets

From the literature survey of open CAI/TAI data, in priority order. Each
needs the inputs listed before it can be scored.

| Dataset | Material | Loading | Needs |
|---|---|---|---|
| Hasebe et al. 2025, *Data in Brief* (Mendeley, CC BY 4.0) | T800S/3900-2B | CAI | T800S/3900-2B preset; separate impact (60x60 mm window) and test (80x50 mm) geometry; cone angle; an assumed impactor mass |
| Sánchez-Sáez et al. 2008 (UC3M open manuscript) | AS4/3501-6 | CAI | Thin, stability-limited 78 mm coupons |
| NASA TP-3102, drop-weight point | T800/3900-2 | CAI | Single point |
| Körbelin 2022, TU Hamburg thesis | M21/T800S | TAI | M21/T800S preset; separate support frame and tension coupon |
| González et al. 2011 and Falcó et al. 2014 (Girona) | AS4/8552 | CAI | Ready now that the AS4/8552 preset exists |

The v0.2.0 roadmap this table replaces named four sources that cannot be
used as listed:
- Soutis & Curtis 1996 is paywalled and re-analyses earlier studies.
- Caprino 1984 is paywalled and covers tension only.
- Sánchez-Sáez et al. 2005 tested AS4/3501-6, not IM7/8552.
- No "NASA/TM-2007" CAI round robin was found on NTRS.

Each case record requires: material name, layup, ply thickness, panel
dimensions, impactor diameter and mass (no defaults: a dataset that omits
them must say what it assumes), impact energy, and measured CAI/TAI
strength (MPa). Optional: `boundary` (`simply_supported`, `clamped`,
`free`), `impactor_shape` (`hemispherical`, `flat`, `conical`), measured
dent depth and DPA.

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
