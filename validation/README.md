# BVID-FE Validation

This directory anchors BVID-FE's validation credibility against published
composite-after-impact test data.

## Datasets

| Dataset | Material | Loading | Cases | Target MAE | Empirical | Semi-analytical | Calibration role |
|---|---|---|---|---|---|---|---|
| `ncamp_as4_8552_cai` | AS4/8552 | CAI | 7 | 15% | 1.5% | 1.5% | In-sample for `soutis_k_s` (pooled fit) and the fe3d zone factor (coupon #0); damage area out of sample |
| `lovejoy_scotti_im7_8552_cai` | IM7/8552 | CAI | 8 | 15% | 8.8% | 8.8% | In-sample for `soutis_k_s`, `olsson_alpha` (areas) and the fe3d zone factor (QI #0, hard #4) |
| `hasebe_t800s_3900_2b_cai` (not gated) | T800S/3900-2B | CAI | 343 | 15% | 21.5% | 23.3% | Areas in-sample (T800S `olsson_alpha`, impactor shape and diameter factors); strengths not fitted |
| `korbelin_m21_t800s_tai` (not gated) | M21/T800S (T800S/3900-2B card) | TAI | 5 | 15% | 3.5% | 3.5% | In-sample for `tai_notch_fraction` (the only TAI fit) |
| `synthetic_selfcheck` | IM7/8552 | CAI + TAI | 12 | 1% | 0.0% (by construction) | 1.2% (not gated) | Tautological (empirical-tier output) |

These are the MAEs `validate_bvid_public.py` prints (rounded). Three
things the gate does not show:
- **The gated strength MAEs are in-sample.** `soutis_k_s` = 4.8 is the
  pooled optimum over NCAMP and Lovejoy & Scotti, so both gated MAEs are
  in-sample; the fair measure is the cross-dataset check under
  [Calibration of `soutis_k_s`](#calibration-of-soutis_k_s) (5.0–8.8%).
  Hasebe et al.'s strengths are not fitted, but the model is known to
  miss them and they are not gated.
- **The semi-analytical growth channel is not exercised by the gate.**
  On every NCAMP and Lovejoy & Scotti coupon the Soutis bound is lower
  than the delamination growth stress (176–177 against 222–224 MPa on
  NCAMP; 156 against 218–219 MPa and 224 against 251–252 MPa on Lovejoy &
  Scotti), so the two tiers give identical numbers there. On Hasebe et
  al.'s set growth governs 238 of 343 coupons (153 of 157 cross-ply
  c8/c16/c24, all 27 hard r0, 58 of 132 quasi-isotropic, no soft r45),
  and on those 238 the semi-analytical MAE is 24.9% against 22.3%
  empirical (both 19.8% on the other 105). The growth channel is therefore unvalidated, and
  where it acts it makes the prediction worse.
- **TAI is not validated.** The only measured tension data are one
  impact condition (`korbelin_m21_t800s_tai`, five coupons), and
  `tai_notch_fraction` is fitted to it; nothing out of sample tests the
  TAI model.

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
reports no dent depth or damage area. The fe3d tier gives 141.0 MPa on
the first coupon (20.6% low; it is one of the three coupons the
damage-zone factor is fitted to, see below). At about 40 s per coupon it
is not run in CI.

**`lovejoy_scotti_im7_8552_cai`** is measured data: Lovejoy & Scotti, NASA
LaRC 2019 (NTRS 20200002432), Tables 1 and 4. It has the eight standard-ply
(190 gsm) coupons, 4 x 6 in, impacted and tested at NIAR to ASTM D7136/D7137:
- four quasi-isotropic [45/0/-45/90]4s coupons (32 plies, 6.1 mm) at 40 ft-lb
  (54 J), 133–147 MPa;
- four hard [45/0/90/0/-45/0/45/0/-45/0]2s coupons (40 plies, 7.6 mm) at
  50 ft-lb (68 J), 236–249 MPa.

Each case carries the measured thickness, dent depth and C-scan damage area.
The areas and dents are recorded, not scored. The predicted areas
(about 5,270 and 6,540 mm²) are in-sample, since `olsson_alpha` is fitted
to their median. The predicted dents (0.40 and 0.45 mm) are 0.95–1.10×
the measured ones, also in-sample (see *Dent depth* below).
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
| 0.17 (then the default) | −3.5% | −4.4% | +1.9% |
| 0.15 | −9.0% | −11.7% | not run |

At 0.17 the three coupons gave 171.4 MPa (177.7 measured), 140.4 MPa
(146.8) and 240.2 MPa (235.7): both laminate families fit within 5%, with
a mean error of 3.3% (18.0% at 0.20 before the confinement). The three
coupons are in-sample. The closed-form tiers barely move with the
confinement: NCAMP 1.54% → 1.50%, Lovejoy & Scotti 8.97% → 8.95%.

Since then the zone is weighted by element area fraction and the failure
index is averaged over a 2.5 mm disk (see `docs/physics_models.md`).
0.17 still fitted: the coupons give 171.9 MPa
(−3.2%), 145.1 MPa (−1.2%) and 246.8 MPa (+4.7%), a mean error of 3.0%,
still in-sample. fe3d has not been run on the other 12 measured NCAMP
and Lovejoy & Scotti coupons, so it has no out-of-sample check.

Three physics fixes then moved the coupons. The first two, measured at
0.17 before the Olsson threshold was corrected (*4.7 to 4.8* below): the
damage template now
takes ply angles modulo 180° (a −45/90 interface is a 45° mismatch, not
135°; each ellipse lies along the back-face-side ply; no delamination
between plies of the same direction): 156.3 MPa (−12.0%), 174.6 MPa
(+18.9%) and 265.6 MPa (+12.7%). The 3D ply stiffness then uses
`ν23 = E22 / (2 G23) − 1` instead of `ν21` for the 2-3 coupling:
152.1 MPa (−14.4%), 162.6 MPa (+10.8%) and 261.2 MPa (+10.8%), a mean
error of 12.0%. The NCAMP coupon now wants a stiffer zone and both
Lovejoy & Scotti coupons a softer one. With the corrected Olsson
threshold, the membrane stiffening in the onset energy and the
mid-surface `u_z` restraint as well (the final code; coupon damage areas
2,783, 5,269 and 6,542 mm²) the factor was swept again (5 mm mesh):

| Factor | NCAMP #0 | Lovejoy & Scotti QI #0 | Lovejoy & Scotti hard #4 | Mean error |
|---|---|---|---|---|
| 0.13 | −26.7% | −6.8% | −7.9% | 13.8% |
| 0.15 (default) | −20.6% | +1.2% | +0.3% | 7.4% |
| 0.17 | −14.6% | +9.0% | +7.9% | 10.5% |
| 0.20 | −5.6% | +20.6% | +11.1% | 12.4% |

0.15 is kept, with the lowest mean error, 7.4%: 141.0 MPa (177.7
measured), 148.6 MPa (146.8) and 236.4 MPa (235.7). No single factor fits
both families within 5%: NCAMP wants about 0.22 and Lovejoy & Scotti
about 0.15, and any factor that moves one family towards its measured
strength moves the other away. The earlier fit of all three within 5%
(0.17: −3.2%, −1.2%, +4.7%) rested on the old template, whose ply
mismatch was not wrapped modulo 180° and whose −45/90 interfaces were
elongated ellipses at 22.5°; with the near-round footprints of the
corrected template the two families part. Which part of the model
separates them (material, thickness, energy or damage area: 2,783 mm² on
the NCAMP coupon against 5,269 and 6,542 mm²) is not established, and
with fe3d run on a single NCAMP coupon the factor cannot be split by
family without over-fitting. The fit is in-sample on three coupons:
Lovejoy & Scotti within about 1%, NCAMP about 21% low. Treat fe3d
strengths as qualitative. fe3d takes about 40 s on the 24-ply NCAMP
coupon and about 80 s on the 40-ply hard laminate (5 mm mesh).

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

**3.5 to 4.7.** Two changes to the model made the datasets agree
(numbers in this step are with the onset energy of the time, about 10×
too low; see *4.7 to 4.8* below):
- The damage area is an energy balance, `α (E − E_onset) / G_IIc` in mm²,
  with α = 0.078 fitted to the median of Lovejoy & Scotti's eight C-scan
  areas. Dividing by the thickness as well left the area flat in
  thickness at constant J/mm, where the measured area grows with it:
  Lovejoy & Scotti's 7.6 mm coupons had 1.21× the C-scan area of their
  6.1 mm coupons at the same 8.9 J/mm. The predicted areas became
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
coupon), rerun with the delaminations confined to the impact window:

| `k_s` | NCAMP | Lovejoy & Scotti |
|---|---|---|
| 4.0 | 11.4% | 15.5% |
| 4.5 | 3.2% | 9.2% |
| 4.7 | 1.5% | 9.0% |
| 5.0 | 4.0% | 8.6% |
| 5.2 | 6.5% | 8.4% |
| 5.5 | 10.2% | 10.4% |

4.7 was the pooled optimum of both tiers. Neither change was enough alone.
With the thickness-divided area and the fixed reference area the two
datasets wanted different `k_s` (NCAMP 7.3% and Lovejoy & Scotti 9.6% at
3.5; 15.7% and 4.6% at 4.0). With the energy-balance area and the panel
area, `k_s` ≈ 4.7 fitted both (3.5% and 9.7%) but put the Hasebe set at
31–34%.

**4.7 to 4.8.** The Olsson threshold load was 3.5× too low (see
`docs/physics_models.md`), so the onset energy was about 10× too low:
0.4 J on the NCAMP coupons and 0.8 J on Lovejoy & Scotti's, now 4.7 J
and 6.1–8.3 J. With the corrected onset `olsson_alpha` was refitted to the
median of Lovejoy & Scotti's C-scan areas (0.078 → 0.0865), the
T800S/3900-2B value and the cone factor to Hasebe et al. (below), and
`k_s` again to the pooled empirical and semi-analytical MAE over NCAMP (7
coupons) and Lovejoy & Scotti (8), weighted by coupon count (the two tiers
agree on both datasets):

| `k_s` | NCAMP | Lovejoy & Scotti | Pooled |
|---|---|---|---|
| 4.5 | 4.9% | 9.1% | 7.1% |
| 4.6 | 3.4% | 9.0% | 6.4% |
| 4.7 | 2.3% | 8.8% | 5.8% |
| 4.8 (default) | 1.6% | 8.7% | 5.4% |
| 4.9 | 1.8% | 8.6% | 5.4% |
| 5.0 | 2.3% | 8.5% | 5.6% |
| 5.2 | 4.9% | 8.3% | 6.7% |
| 5.5 | 8.5% | 10.1% | 9.4% |

Keeping 4.7 with the refitted `olsson_alpha` would cost 0.8 points on
NCAMP. The NCAMP damage area falls from about 2,900 to 2,790 mm².

**4.8 kept (membrane-stiffened onset).** The onset energy then gained the
plate's membrane stiffness (`P = k_b w + k_m w³`) and the exact contact
energy (2/5 of `P_c δ` for Hertz instead of 1/2), see
`docs/physics_models.md`. On these thick coupons (deflection at `P_c`
0.1–0.3 of the thickness) the onset falls by 8–9%, mostly from the
contact energy: 4.65 → 4.27 J (NCAMP) and a median 7.2 → 6.5 J (Lovejoy &
Scotti). With the same procedure `olsson_alpha` refits to 0.0855 (from
0.0865; median of Lovejoy & Scotti's areas), the T800S/3900-2B value to
0.0605 and the cone factor to 1.9 (below), and the `k_s` scan again
puts the pooled optimum at 4.8:

| `k_s` | NCAMP | Lovejoy & Scotti | Pooled |
|---|---|---|---|
| 4.5 | 4.8% | 9.1% | 7.1% |
| 4.6 | 3.3% | 9.0% | 6.3% |
| 4.7 | 2.2% | 8.9% | 5.8% |
| 4.8 (default) | 1.5% | 8.8% | 5.4% |
| 4.9 | 1.8% | 8.7% | 5.5% |
| 5.0 | 2.4% | 8.6% | 5.7% |
| 5.2 | 5.0% | 8.3% | 6.8% |
| 5.5 | 8.6% | 10.1% | 9.4% |

Before the refit (linear plate, the constants above) the MAEs were
1.55% (NCAMP), 8.73% (Lovejoy & Scotti) and 21.19% / 23.07% (Hasebe et
al., empirical / semi-analytical); after it they are 1.54%, 8.78% and
21.53% / 23.35%. The NCAMP damage area is about 2,790 mm² either way.

The in-sample numbers flatter the model. The fair measure is
cross-dataset:
- `k_s` fitted on NCAMP alone (4.8) predicts Lovejoy & Scotti at 8.8% MAE.
- `k_s` fitted on Lovejoy & Scotti alone (5.2) predicts NCAMP at 5.0%.

Expect errors of about 5–9% on a new material or layup tested to ASTM
D7136/D7137, and more on a matrix other than 8552: both datasets share
it. `olsson_alpha` was fitted on Lovejoy & Scotti's areas, so only the
NCAMP damage area is out of sample, and NCAMP reports no area to check
it against.

An equivalent-open-hole model (Whitney-Nuismer point stress with a
finite-width correction) was also tried, with the thickness-divided damage
area. Its characteristic distance calibrates cleanly from the NCAMP
open-hole compression data (d0 = 0.7–1.8 mm). But the predicted damage
spanned 78–126 mm of a 100 mm wide coupon, so the equivalent hole filled
the specimen. It has not been retried with the energy-balance area (about
2,800 mm² on the NCAMP coupons, against about 6,650 mm² then).

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
  (conical), where the unsourced previous values were 1.4 and 0.7. With
  the corrected onset energy the tip also changes the onset (through the
  contact stiffness; the linearised cone contact is the softest, so the
  cone's onset is the highest), and the factors were refitted so that,
  in each same-layup, same-energy group, the model's predicted / measured
  area for a tip matches that of the 12.7 mm hemisphere (median over
  groups). That ratio moved from 0.99 to 0.84 for the cone (37 groups)
  and from 0.83 to 0.95 for the flat-ended impactor (11 groups), giving
  1.54 (then 1.5) and 0.21 (kept at 0.2). For the 25.4 and 31.8 mm
  hemispheres it moved by under 5% (1.01 → 1.04, 30 groups; 1.25 → 1.27,
  6 cross-ply groups), so the exponent stays 0.2. With the
  membrane-stiffened onset and the exact contact energy (the elastic
  cone's contact now stores a third of `P_c δ`, not a quarter, which
  raises the cone's onset further) the same procedure moved the cone's
  ratio from 0.97 to 0.81 and the flat-ended impactor's from 0.95 to
  0.92, giving 1.86 (now 1.9) and 0.22 (kept at 0.2); the 25.4 and
  31.8 mm hemispheres moved to 1.03 and 1.26 (the fitted exponent would be
  0.30, but as before a change under 5% keeps 0.2).

The `T800S/3900-2B` preset takes its elastic constants and toughness from
Hasebe's Table 1 and its strengths from the Toray 3900 data sheet (Hasebe's
Xc of 719 MPa is under half the data sheet's ASTM D6641 1779 MPa).
`olsson_alpha` = 0.0605 is fitted to the C-scan areas of the
quasi-isotropic laminates struck by the 12.7 mm hemisphere (median of 38
specimens; the default 0.0855 overpredicts them about 1.4×, as the
material is interlayer-toughened; it was 0.0648 against a default of
0.0865 with the linear plate, and 0.047 against 0.078 while the onset
energy was about 10× too low). It predicts the area of the 329 specimens
with a measured area (11 flat-impactor and 3 25.4 mm hemisphere
specimens show none) at a median 1.03× (per-layup medians 0.77–1.26×,
the highest on the 8-ply cross-ply c8; per-impactor 0.94–1.07×), and the
mean area per impactor over all 343 within 2% for the flat-ended
impactor and the 12.7 and 25.4 mm hemispheres, 13% high for the cone and
18% low for the 31.8 mm hemisphere: 410 against 364 mm² for the 60° cone,
52 against 52 mm² for the flat-ended impactor, and 282 against 281, 280
against 277 and 205 against 251 mm² for the 12.7, 25.4 and 31.8 mm
hemispheres. The onset energy is a median 2.3 J here (0.9–7.3 J): three
24-ply cross-ply specimens struck by the cone at 1.1 J/mm (5.0 J) fall
below it (6.9–7.1 J, most of it the elastic cone's contact energy) and
are predicted undamaged, against 62–72 mm² measured. The 14 specimens
with no measured damage were struck at 2.2–5.9× the predicted onset,
so the data do not test the onset itself.
These area checks are in-sample: `olsson_alpha` and the impactor shape
and diameter factors are fitted to the same specimens. The predicted
dent is a median 0.94× the measured dent over the 308 with a positive one
(in-sample, see *Dent depth* below).

Why the model misses: the damage area is right, but the knockdown map is
not. A constant knockdown (the set's mean, 0.73) scores 30.3%, against
21.5% (empirical) and 23.3% (semi-analytical). The empirical MAE per
layup is 14% (c8), 31% (c16), 18% (c24), 20% (q8), 24% (q16), 37% (q24),
15% (r0) and 4% (r45). The predicted knockdown averages 0.63 against 0.73
measured. Two things remain:
- Layup. At similar delamination area, cross-ply c16 keeps 0.64 of its
  strength and quasi-isotropic q16 0.84; the DPA-based knockdown gives
  both about 0.64.
- Thin, fixture-limited coupons. The undamaged reference is itself
  limited by the JIS fixture (see above), and a `k_s` fitted per layup
  rises with thickness within each family (q8 1.4, q16 1.5, q24 3.2).
  That looks like the fixture, not the material.

A layup factor from Lekhnitskii's `Kt_inf` was tried and helps little
(about 1.6 points with the measured areas) and points the wrong way for
the soft r45 laminate, so no layup term was added. Replicate scatter alone
accounts for about 6–9%.

**`korbelin_m21_t800s_tai`** is measured TAI data: Körbelin, Dreiner &
Fiedler, *Composites Part C: Open Access* 3 (2020) 100074
(doi:10.1016/j.jcomc.2020.100074, CC BY 4.0; TUHH Open Research
hdl 11420/8990). Hexcel M21/T800S [45/0/−45/90]2s, 4.08 mm (16 plies of
0.255 mm). The tension coupons are 250 × 60 mm, impacted on a 40%
scaled-down ASTM D7136 support (a 50 × 30 mm window, `impact_panel_*`)
by a 5.34 kg, 12 mm hemispherical impactor at 25 J (Table 2; the text
says 12.5 mm), and tested at 20 °C in a jig that clamps 5 mm of each
edge against free-edge delamination. The five coupons impacted at 20 °C
are scored; each strength, and the mean of the five unimpacted coupons
(1045.3 MPa, `measured_pristine_MPa`), is read from the coordinates of
the vector markers of the article's Fig. 17, with the y axis set by its
700–1100 MPa tick labels. The digitized means reproduce the reductions
the text states (11.1% at 20 °C, 18.58% at 80 °C, and 25.01% for the
three lowest 80 °C coupons). The mean C-scan area (1,199 mm², Fig. 8)
and dent (0.151 mm, Fig. 7) are recorded on every case, not scored. The
five coupons impacted at 80 °C, three with fibre failure on the impacted
face (736.8–958.5 MPa), are left out, because the model has no impact
temperature.

There is no M21/T800S preset, so the T800S/3900-2B card stands in (same
fibre, also interlayer-toughened). The cases are scored on knockdown ×
`measured_pristine_MPa`, so the card enters only through the damage
mapping and `K_T∞` (3.0, quasi-isotropic), and the predicted footprint
(878 mm², against the measured 1,199) fills the 30 mm window either way:
with `G_IIc` = 1.0 N/mm and the default `olsson_alpha` the area is
1,119 mm² and the width and knockdown do not change.

The coupons keep 0.889 of their strength on average (0.85–0.97).
`tai_notch_fraction` is fitted to that mean: with `wn_d0_mm` = 1 mm the
30 mm footprint needs a 1.4 mm hole, `f` = 0.046; 0.05 is adopted,
knockdown 0.875, MAE 3.5% for both closed-form tiers (the coupons'
scatter about their mean is 3.7%). The fit is in-sample on one energy,
one layup and one thickness, and with the footprint fixed by the window
it constrains only `f · w / d₀`. The previous TAI model, a hole of the
whole damage area in an infinite plate, gave 0.380 here (MAE 57.2%).

The fe3d tier gives 0.335 here (MAE 62.2%). fe3d TAI is not validated
and not recommended: its damage zone's strain-equivalent strength
reduction applies in every in-plane direction, so in tension the soft
zone fails as if its fibres were cut, while delaminations cut none, and
its damage-zone factor is fitted on CAI coupons only. On 150 × 100 mm
IM7/8552 panels its tension knockdown is about one third of the
closed-form TAI at 15–30 J.

Other TAI sources looked at and not used: NASA TP-2981 (Nettles 1990)
tests eight-ply coupons punctured over a 10.3 mm hole with a 4.2 mm tup,
reports its results only as plots and uses no preset material; Lovejoy &
Scotti 2019 has no TAI tests; Kötter's 2021 TUHH thesis has fatigue after
impact at one energy only; Itakura et al. 2022 (JSME) is gas-gun
penetration; Poe's NASA reports cover thick filament-wound rocket motor cases;
Körbelin's 2022 TUHH thesis ("Damage tolerance of high-performance
composites", Band 41) was not located online.

**`synthetic_selfcheck`** is tautological: each case is produced by running
the empirical tier at a known input, so its MAE is ~0% by construction.
The semi-analytical tier scores 1.2% on it (TAI is the same model, and
Soutis governs every compression case but the 5 J one, where the
delamination growth stress is 14.5% lower), but CI does not run that
tier on it. It
exercises the harness end-to-end and gates regressions in the empirical
pipeline. Regenerate it when the empirical tier changes on purpose; there
is no generator script, so rerun each case through the empirical tier and
write the predicted residual strength into `measured_strength_MPa`.

### Dent depth

The validator prints, for each dataset that records `measured_dent_mm`,
the median predicted / measured dent (and its 10–90% band) over the cases
with a positive measured dent. It is reported, not gated.

The dent law `d = dent_beta (E − E_onset)^dent_gamma · S_shape`
(`src/bvidfe/impact/dent_model.py`, `docs/physics_models.md`) is fitted
as follows:

1. Take every Lovejoy & Scotti and Hasebe et al. case with a measured dent
   (351; 316 positive, 34 of the 40 flat-ended Hasebe dents and one
   hemispherical one are recorded as 0). Compute `E_onset` with the code
   (`bvidfe.impact.olsson.onset_energy` on the case's impact panel,
   boundary and impactor, at the panel centre), so a later onset change
   is refitted by rerunning the fit.
2. `dent_beta`, `dent_gamma`: least squares on
   `log d = log dent_beta + dent_gamma log(E − E_onset)` over the
   hemispherical cases with a positive dent above onset (8 Lovejoy &
   Scotti, 210 Hasebe), each dataset weighted equally (weight
   1 / n_dataset). Result 0.0441 mm/J^0.569 with the linear-plate onset
   energy, used as 0.044 and 0.57 (unweighted: 0.0416, 0.609). With the
   membrane-stiffened onset the refit is 0.0390 mm/J^0.601 (unweighted
   0.0338, 0.686); the constants as used still score within 6% of 1 per
   dataset (table below), so they are kept.
3. Shape factors (Hasebe only), relative to the hemispherical law: the 60°
   cone's median measured / predicted (8.57 with the linear-plate onset,
   used as 8.5; 9.83 against the refitted law above), and the flat-ended
   impactor's mean measured over mean predicted (0.027, used as 0.03;
   0.028 refitted; most of its dents are 0, so a median or log fit is
   undefined).
4. Score: median predicted / measured and log-sd per dataset, material,
   impactor and layup; fit on one dataset and predict the other; refit
   with `E_onset` scaled by 0.5 and 1.5.

Results with the constants as used:

| Cases | n | Median pred / meas | 10–90% | log-sd |
|---|---|---|---|---|
| Lovejoy & Scotti, IM7/8552, hemispherical | 8 | 1.01 | 0.96–1.10 | 0.06 |
| Hasebe et al., T800S/3900-2B, all positive dents | 308 | 0.94 | 0.52–1.44 | 0.57 |
| Hasebe, 12.7 / 25.4 / 31.8 mm hemisphere | 98 / 95 / 17 | 0.78 / 1.12 / 1.44 | | 0.45 / 0.50 / 0.62 |
| Hasebe, 60° cone | 92 | 0.91 | 0.47–1.12 | 0.44 |
| Hasebe, flat-ended (positive dents only) | 6 | 0.18 | | |

Per Hasebe layup (hemispherical) the median is 0.74 (q8) to 1.13 (c8).
The mean flat-ended dent is 0.0039 mm predicted against 0.0034 mm
measured. Cross-dataset: fitted on Hasebe alone (0.0313, 0.732) the law
predicts Lovejoy & Scotti at 1.38×; fitted on Lovejoy & Scotti with
Hasebe's exponent (0.0225), Hasebe's hemispherical dents at 0.71×
(Lovejoy & Scotti's two energies alone give an exponent of 0.81 and
Hasebe at 0.61×). The earlier form scaled by `G_Ic` and thickness,
`d / h = β ((E − E_onset) / (G_Ic h²))^γ`, refitted the same way, gives
1.54× on Lovejoy & Scotti from Hasebe, so neither is in the law and one
pair serves both materials. With `E_onset` halved the refit is (0.034,
0.63) and the constants as used score 1.04 / 1.09 (Lovejoy & Scotti /
Hasebe); with it 1.5× higher, (0.048, 0.55) and 0.97 / 0.77. The `h / 2`
cap binds on 7 of the 351 cases, all 8-ply specimens struck by the cone at
5.1–5.2 J (uncapped, up to 1.04× the cap), whose measured dents
(1.3–1.9 mm) are themselves deeper than `h / 2`; 18 thin-laminate cone
dents exceed `h / 2`.

Before the fit the law was `h · 0.05 · ((E − E_onset) · 10³ / (G_Ic h²))^0.5`,
on the `h / 2` cap for 345 of the 351 cases (the other six were three
below onset and three 4.6 mm cross-ply specimens at 0.95 of the cap):
median 8.2× (7.5–8.8×) on Lovejoy & Scotti and 9.0× on Hasebe. Its units
were consistent (J → N·mm against `G_Ic` in N/mm); the constants had never
been fitted, and uncapped it was 33× (Hasebe) and 51× (Lovejoy & Scotti)
the measured dents.

CI gates the empirical tier on every dataset and the semi-analytical tier on
the measured ones, at 1.25x each dataset's target. A dataset may set
`"gate": false` with a `gate_note` for a known miss: it is still run and its
error printed, but it cannot fail the gate. `hasebe_t800s_3900_2b_cai` does. The fe3d tier
is run on the first three synthetic cases as a smoke run, not a gate: the
synthetic set is empirical-tier output, so its error there (29.2%) measures
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
python validation/validate_bvid_public.py --tier semi_analytical --dataset ncamp_as4_8552_cai lovejoy_scotti_im7_8552_cai hasebe_t800s_3900_2b_cai korbelin_m21_t800s_tai
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
   `.github/workflows/tests.yml`; if it is gated (no `"gate": false`),
   also add it to `test_measured_cai_gates_pass_for_the_closed_form_tiers`
   in `tests/test_validation_smoke.py`.
6. Commit the dataset, any recalibration, and the updated numbers in a
   single PR.
