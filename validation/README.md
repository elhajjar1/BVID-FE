# BVID-FE Validation

This directory anchors BVID-FE's validation credibility against published
composite-after-impact test data.

## Datasets

| Dataset | Material | Loading | Cases | Target MAE | Empirical | Semi-analytical |
|---|---|---|---|---|---|---|
| `ncamp_as4_8552_cai` | AS4/8552 | CAI | 7 | 15% | 1.5% | 1.5% |
| `lovejoy_scotti_im7_8552_cai` | IM7/8552 | CAI | 8 | 15% | 9.0% | 9.0% |
| `hasebe_t800s_3900_2b_cai` (not gated) | T800S/3900-2B | CAI | 343 | 15% | 21.3% | 23.3% |
| `synthetic_selfcheck` | IM7/8552 | CAI + TAI | 12 | 1% | ~0% (by construction) | n/a |

**`ncamp_as4_8552_cai`** is measured data: NCAMP CAM-RP-2010-002 Rev A
(2011), the Hexcel 8552 AS4 unitape qualification report, sections 2.3.31
and 4.31. It has seven [45/0/-45/90]3s coupons (24 plies, 4.5 mm), tested
to ASTM D7136/D7137 at RTD and impacted at 1500 in-lbf/in (about 30 J). The
measured CAI is 171–182 MPa and every coupon failed through the impact damage
(failure code LDM). Each case uses its own measured thickness and impact
energy. The report does not give the impactor mass or tup, so the ASTM
D7136 standard values (5.5 kg, 0.625 in hemispherical) are used, and the
plate is impacted on the standard's 125 × 75 mm support window
(`impact_panel_*`), which also confines the predicted delaminations. It
reports no dent depth or damage area. The fe3d tier gives 171 MPa on the
first coupon (3.5% low; its damage-zone factor is fitted to this coupon,
see below). At about 3 minutes per coupon it is not run in CI.

**`lovejoy_scotti_im7_8552_cai`** is measured data: Lovejoy & Scotti, NASA
LaRC 2019 (NTRS 20200002432), Tables 1 and 4. It has the eight standard-ply
(190 gsm) coupons, 4 x 6 in, impacted and tested at NIAR to ASTM D7136/D7137:
- four quasi-isotropic [45/0/-45/90]4s coupons (32 plies, 6.1 mm) at 40 ft-lb
  (54 J), 133–147 MPa;
- four hard [45/0/90/0/-45/0/45/0/-45/0]2s coupons (40 plies, 7.6 mm) at
  50 ft-lb (68 J), 236–249 MPa.

Each case carries the measured thickness, dent depth and C-scan damage area.
The impactor is not described, so the ASTM D7136 standard is assumed.
Like the NCAMP cases, each is impacted on the ASTM D7136 125 × 75 mm
support window (`impact_panel_*`), which also confines the predicted
delaminations. The hybrid thin-ply coupons are omitted, because a case has
a single ply thickness.

### fe3d damage-zone factor

fe3d softens the in-plane stiffness of every ply under the damage
footprint by `DAMAGE_ZONE_INPLANE_FACTOR`. It was an unfitted 0.30. With
the energy-balance damage area the zone shrank, and 0.30 put the
quasi-isotropic coupons 30–36% high. A first sweep on three coupons (error
of the fe3d strength against the measured CAI), with the delaminations
not yet confined:

| Factor | NCAMP #0 | Lovejoy & Scotti QI #0 | Lovejoy & Scotti hard #4 |
|---|---|---|---|
| 0.30 | +30.2% | +35.5% | −14.1% |
| 0.25 | +18.3% | +16.7% | −31.1% |
| 0.20 | +4.8% | −2.0% | −47.1% |
| 0.15 | −8.9% | −20.4% | −62.1% |

No factor fitted both families: the hard laminate wanted more than 0.30
while the quasi-isotropic coupons wanted about 0.20. The cause was the
footprint, not the factor. The union of the template's differently
oriented ellipses spreads well beyond the circle of the same area: on the
hard laminate it was 108 × 108 mm (same-area circle 92 mm) on a 102 mm
wide coupon, so fe3d's soft zone spanned the whole section and the 0°
fibres failed at the coupon's side edge (2.4 mm from it). On the
quasi-isotropic coupons it was 99 × 52 mm (NCAMP, circle 61 mm) and
134 × 70 mm (Lovejoy & Scotti, circle 82 mm). The delaminations are now
kept inside the impacted panel, here the ASTM D7136 125 × 75 mm support
window, and the sweep with that confinement gives:

| Factor | NCAMP #0 | Lovejoy & Scotti QI #0 | Lovejoy & Scotti hard #4 |
|---|---|---|---|
| 0.20 | +4.7% | +6.6% | +12.4% |
| 0.17 (default) | −3.5% | −4.4% | +1.9% |
| 0.15 | −9.0% | −11.7% | not run |

At 0.17 the three coupons give 171.4 MPa (177.7 measured), 140.4 MPa
(146.8) and 240.2 MPa (235.7): both laminate families fit within 5%, with
a mean error of 3.3% (18.0% at 0.20 before the confinement). The three
coupons are in-sample. The closed-form tiers barely move with the
confinement: NCAMP 1.54% → 1.50%, Lovejoy & Scotti 8.97% → 8.95%.

### Calibration of `soutis_k_s`

The empirical CAI knockdown is `1 / (1 + k_s (DPA/A_ref)^m)`, with
`A_ref` = 15,000 mm², the 150 × 100 mm ASTM D7137 coupon. `m` = 0.5 has
been kept; `soutis_k_s` has been calibrated twice.

**2.5 to 3.5.** `soutis_k_s` was an uncalibrated 2.5, which overpredicted
the Lovejoy & Scotti strengths by 29–56% although the undamaged strength
was right (601 MPa against 600 MPa measured by NCAMP). The damage area
was then `α (E − E_onset) / (G_IIc h)` with α = 0.8, and the knockdown
scaled it by the analysed panel's area. 3.5 fitted both datasets:

| `k_s` | NCAMP empirical | Lovejoy empirical | NCAMP semi-analytical | Lovejoy semi-analytical |
|---|---|---|---|---|
| 2.5 | 15.9% | 39.9% | 11.6% | 23.5% |
| 3.5 | 7.3% | 10.9% | 7.3% | 10.8% |

Fitted on NCAMP alone (3.1) it predicted Lovejoy & Scotti at 16–21% MAE,
and fitted on Lovejoy & Scotti alone (4.0) it predicted NCAMP at 15.7%.
The Lovejoy & Scotti damage area was overpredicted (about 8,900 mm²
against 5,190–7,530 mm²).

**3.5 to 4.7.** Two changes to the model made the datasets agree:
- The damage area is an energy balance, `α (E − E_onset) / G_IIc` in mm²,
  with α = 0.078 fitted to the median of Lovejoy & Scotti's eight C-scan
  areas. Dividing by the thickness as well left the area flat in
  thickness at constant J/mm, where the measured area grows with it:
  Lovejoy & Scotti's 7.6 mm coupons had 1.21× the C-scan area of their
  6.1 mm coupons at the same 8.9 J/mm. The predicted areas are now
  5,304 mm² (quasi-isotropic) and 6,616 mm² (hard), against 5,190–5,740
  and 6,260–7,530 mm² measured. NCAMP reports no damage area, so its
  predicted area (about 2,900 mm², was about 6,650 mm²) is an independent
  prediction.
- The knockdown scales the damage area by the fixed coupon area instead of
  the analysed panel's. The panel area sent the knockdown to 1 as the
  panel grew, although impact damage is a local notch, and over-penalised
  Hasebe et al.'s 80 × 50 mm coupons.

With both, the empirical MAE against `k_s` (the semi-analytical tier gives
the same on these datasets, because the Soutis bound governs every
coupon):

| `k_s` | NCAMP | Lovejoy & Scotti |
|---|---|---|
| 4.0 | 11.8% | 15.8% |
| 4.5 | 3.5% | 9.2% |
| 4.7 (default) | 1.5% | 9.0% |
| 5.0 | 3.6% | 8.6% |
| 5.5 | 9.8% | 10.1% |

4.7 is the pooled optimum of both tiers. Neither change was enough alone.
With the thickness-divided area and the fixed reference area the two
datasets wanted different `k_s` (NCAMP 7.3% and Lovejoy & Scotti 9.6% at
3.5; 15.7% and 4.6% at 4.0). With the energy-balance area and the panel
area, `k_s` ≈ 4.7 fitted both (3.5% and 9.7%) but put the Hasebe set at
31–34%.

The in-sample numbers flatter the model. The fair measure is
cross-dataset:
- `k_s` fitted on NCAMP alone (4.7) predicts Lovejoy & Scotti at 9.0% MAE.
- `k_s` fitted on Lovejoy & Scotti alone (5.2) predicts NCAMP at 6.2%.

Expect errors of about 6–9% on a new material or layup tested to ASTM
D7136/D7137. `olsson_alpha` was fitted on Lovejoy & Scotti's areas, so
only the NCAMP damage area is out of sample.

An equivalent-open-hole model (Whitney-Nuismer point stress with a
finite-width correction) was also tried, with the thickness-divided damage
area. Its characteristic distance calibrates cleanly from the NCAMP
open-hole compression data (d0 = 0.7–1.8 mm). But the predicted damage
spanned 78–126 mm of a 100 mm wide coupon, so the equivalent hole filled
the specimen. It has not been retried with the energy-balance area (about
2,900 mm² on the NCAMP coupons, against about 6,650 mm² then).

**`hasebe_t800s_3900_2b_cai`** is measured data the model is known to miss,
so it is not gated: Hasebe, Higuchi, Yokozeki & Takeda, *Data in Brief* 60
(2025) 111509, with the per-specimen spreadsheets on Mendeley Data
(10.17632/8scdmfdcfb.3, CC BY 4.0). Toray T800S/3900-2B in 8 layups:
cross-ply c8/c16/c24 [0/90]ns, quasi-isotropic q8/q16/q24 [45/0/-45/90]ns,
hard r0 and soft r45 (20 plies). Plates 80 × 80 mm were impacted clamped on
a 60 × 60 mm window (`impact_panel_*`, `impact_boundary`), then trimmed to
80 × 50 mm and tested to JIS K 7089. The 343 specimens struck by a
hemispherical impactor (12.7, 25.4 and 31.8 mm), the 15.9 mm flat-ended
impactor or the 60° cone (2–21 J) are carried, each with its measured
thickness, coupon size, C-scan area and dent. Notes:
- The impactor mass is not reported; 5.5 kg is assumed (DPA scales with
  `(5.5/m)^0.1`).
- In the impact sheet the `[J/mm]` and `[J]` labels are swapped;
  `impact_energy_J` is the column equal to energy per mm × thickness.
- The undamaged coupons fail at 0.4–0.6× the unnotched laminate strength,
  and lower the thinner they are (261–277 MPa for 8 plies, 434–528 MPa for
  24), so the JIS fixture, not the material, limits them. Each case
  therefore carries `measured_pristine_MPa` (the layup's mean undamaged
  strength, the basis of the dataset's own reduction ratio) and is scored
  on predicted knockdown × that strength.
- The 120° cone is left out (the model's cone has a 30° half-angle). The
  cones' shank diameter is not given; 12.7 mm is entered, which the model
  does not use for a conical impactor.
- The model's impactor factors are fitted to these data, the only open set
  that varies the impactor on one material. Over same-layup, same-energy
  groups, against the 12.7 mm hemisphere, the 25.4 and 31.8 mm hemispheres
  left 0.88× and 0.87× the C-scan area (diameter exponent 0.2, was 0.3),
  the 15.9 mm flat-ended impactor 0.20× and the 60° cone 1.20×. Relative
  to a 16 mm hemisphere that is a shape factor of 0.2 (flat) and 1.3
  (conical), where the unsourced previous values were 1.4 and 0.7.

The `T800S/3900-2B` preset takes its elastic constants and toughness from
Hasebe's Table 1 and its strengths from the Toray 3900 data sheet (Hasebe's
Xc of 719 MPa is under half the data sheet's ASTM D6641 1779 MPa).
`olsson_alpha` = 0.047 is fitted to the C-scan areas of the
quasi-isotropic laminates struck by the 12.7 mm hemisphere (median of 38
specimens; the default 0.078 overpredicts them about 1.7×, as the
material is interlayer-toughened). With the energy-balance area it
predicts the area of all 343 specimens at a median 1.03× (0.81–1.31× per
layup), and the mean area per impactor within 8%, except the 31.8 mm
hemisphere (19% low): 368 against 364 mm² for the 60° cone, 48 against
52 mm² for the flat-ended impactor, and 277 against 281, 266 against 277
and 204 against 251 mm² for the 12.7, 25.4 and 31.8 mm hemispheres.

Why the model misses: the damage area is right, but the knockdown map is
not. A constant knockdown (the set's mean, 0.73) scores 30.3%, against
21.3% (empirical) and 23.3% (semi-analytical). The empirical MAE per
layup is 13% (c8), 33% (c16), 19% (c24), 19% (q8), 23% (q16), 36% (q24),
14% (r0) and 5% (r45). The predicted knockdown averages 0.64 against 0.73
measured. Two things remain:
- Layup. At similar delamination area, cross-ply c16 keeps 0.64 of its
  strength and quasi-isotropic q16 0.84; the DPA-based knockdown gives
  both about 0.64.
- Thin, fixture-limited coupons. The undamaged reference is itself
  limited by the JIS fixture (see above), and a `k_s` fitted per layup
  rises steeply with thickness within each family (q8 0.8, q16 1.2,
  q24 2.5). That looks like the fixture, not the material.

A layup factor from Lekhnitskii's `Kt_inf` was tried and helps little
(about 1.6 points with the measured areas) and points the wrong way for
the soft r45 laminate, so no layup term was added. Replicate scatter alone
accounts for about 6–9%.

**`synthetic_selfcheck`** is tautological: each case is produced by running
the empirical tier at a known input, so its MAE is ~0% by construction. It
exercises the harness end-to-end and gates regressions in the empirical
pipeline. Regenerate it when the empirical tier changes on purpose; there
is no generator script, so rerun each case through the empirical tier and
write the predicted residual strength into `measured_strength_MPa`.

CI gates the empirical tier on every dataset and the semi-analytical tier on
the measured ones, at 1.25x each dataset's target. A dataset may set
`"gate": false` with a `gate_note` for a known miss: it is still run and its
error printed, but it cannot fail the gate. `hasebe_t800s_3900_2b_cai` does. The fe3d tier
is run on the first three synthetic cases as a smoke run, not a gate: the
synthetic set is empirical-tier output, so its error there (23.8%) measures
the difference between the tiers, not fe3d's accuracy.

## Unnotched (pristine) strength

Every tier divides by the same pristine strength, the unnotched laminate
strength from ply-discount last-ply failure on classical lamination
theory. `tests/validation/test_unnotched_strength_ncamp.py` checks it
against the NCAMP RTD unnotched tension and compression strengths of
AS4/8552 (CAM-RP-2010-002) and IM7/8552 (CAM-RP-2009-015) in the
25/50/25, 10/80/10 and 50/40/10 layups, twelve values in all
(predicted / measured):

| Pristine strength | Range | Mean | Mean \|log error\| |
|---|---|---|---|
| Thickness-weighted ply-strength average (v0.2.0) | 1.17–2.85× | 1.75× | |
| Laminate strength, unidirectional-coupon `Xt` | 0.83–1.19× | 1.00× | 10.9% |
| Laminate strength, NCAMP cross-ply back-out `Xt` (current) | 0.83–1.14× | 0.97× | 8.4% |

The 0.2.0 average treated the plies as sharing stress; they share strain,
so the off-axis plies were credited with strength they never reach before
the 0° fibres break. With the unidirectional-coupon `Xt` (AS4/8552
1996 MPa, IM7/8552 2560 MPa) the 0° plies outlived the ±45 cracking and
the fibre-governed tension strengths came out 5–19% high; NCAMP's
cross-ply back-out value (1872 and 2238 MPa) is the 0° ply strength
realised inside a laminate and puts the 50/40/10 laminates at
0.99–1.02×. Where the ±45 plies cracking in shear sets the strength, one
`S12` cannot fit both families: 25/50/25 tension is 4% (IM7/8552) and 14%
(AS4/8552) high, 10/80/10 tension 10% low and 10/80/10 compression
15–17% low.

## Next datasets

From the literature survey of open CAI/TAI data, in priority order. Each
needs the inputs listed before it can be scored.

| Dataset | Material | Loading | Needs |
|---|---|---|---|
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

Each case record requires: `material` (a preset name), `layup_deg`,
`ply_thickness_mm`, `panel_Lx_mm`, `panel_Ly_mm`, `impactor_diameter_mm`
and `impactor_mass_kg` (no defaults: a dataset whose source omits them
must say in its `reference` what it assumes), `impact_energy_J`, and
`measured_strength_MPa`. Optional: `loading` (`compression`, the default,
or `tension`), `boundary` (`simply_supported`, the default, `clamped`,
`free`), `impactor_shape` (`hemispherical`, the default, `flat`,
`conical`), `measured_dent_mm` and `measured_dpa_mm2` (recorded, not
scored), `measured_pristine_MPa` (scores the case on predicted knockdown ×
this undamaged strength, for fixtures that cap the undamaged coupon),
`impact_panel_Lx_mm`, `impact_panel_Ly_mm` and `impact_boundary` (the
window the plate was impacted on, when it is not the test panel), and
`note`.

## Running

```bash
python validation/validate_bvid_public.py              # empirical tier, all datasets, print table
python validation/validate_bvid_public.py --gate       # exit non-zero if a gated MAE exceeds 1.25 * target
python validation/validate_bvid_public.py --dataset synthetic_selfcheck
python validation/validate_bvid_public.py --tier semi_analytical --dataset ncamp_as4_8552_cai lovejoy_scotti_im7_8552_cai hasebe_t800s_3900_2b_cai
python validation/validate_bvid_public.py --tier fe3d --dataset synthetic_selfcheck --max-cases 3
```

`--tier` is `empirical` (default), `semi_analytical` or `fe3d`. The gate
multiplier is 1.25 for the closed-form tiers and 2.0 for fe3d.
`--dataset` takes one or more file stems; a name that matches no file is
an error (exit 1). `--max-cases N` runs only the first N cases of each
dataset.

## Adding a new dataset

1. Create `datasets/<name>.json` with top-level `name`, `reference` (the
   citation, tables used, and every assumed value), `target_mae_pct`
   (defaults to 100 if omitted) and `cases` (fields above; they map onto
   the `DatasetCase` dataclass in `validate_bvid_public.py`). The material
   must be a preset in `MATERIAL_LIBRARY`; add one first if needed (cited
   constants only).
2. The validator auto-discovers any `*.json` file in `datasets/`, so the
   CI empirical gate picks it up. Run both closed-form tiers on it and
   record the MAE in the table above.
3. If the model is known to miss it, set `"gate": false` and a
   `gate_note` saying why, rather than loosening `target_mae_pct`.
4. Do not tune calibration constants (`soutis_k_s`, ...) on the new
   dataset alone. If it motivates a recalibration, fit on part of the data
   and report the error on data held out from the fit (the cross-dataset
   check above), and rerun every other dataset.
5. For a measured dataset, add its name to the semi-analytical step in
   `.github/workflows/tests.yml` and to
   `test_measured_cai_gates_pass_for_the_closed_form_tiers` in
   `tests/test_validation_smoke.py`.
6. Commit the dataset, any recalibration, and the updated numbers in a
   single PR.
