# BVID-FE Validation

This directory anchors BVID-FE's validation credibility against published
composite-after-impact test data.

## Datasets

| Dataset | Material | Loading | Cases | Target MAE | Empirical | Semi-analytical | Calibration role |
|---|---|---|---|---|---|---|---|
| `ncamp_as4_8552_cai` | AS4/8552 | CAI | 7 | 15% | 1.5% | 1.5% | In-sample for `soutis_k_s` (pooled fit) and the fe3d zone factor (coupon #0); damage area out of sample |
| `lovejoy_scotti_im7_8552_cai` | IM7/8552 | CAI | 8 | 15% | 8.8% | 8.8% | In-sample for `soutis_k_s`, `olsson_alpha` (areas) and the fe3d zone factor (QI #0, hard #4) |
| `hasebe_t800s_3900_2b_cai` (not gated) | T800S/3900-2B | CAI | 343 | 15% | 21.5% | 21.5% | Areas in-sample (T800S `olsson_alpha`, impactor shape and diameter factors); strengths not fitted |
| `korbelin_m21_t800s_tai` (not gated) | M21/T800S (T800S/3900-2B card) | TAI | 5 | 15% | 3.5% | 3.5% | In-sample for `tai_notch_fraction` (the only TAI fit) |
| `gonzalez_as4_8552_cai` (reported, not gated) | AS4/8552 | CAI | 19 | 15% | 14.3% | 14.3% | Out of sample: nothing is fitted to it |
| `synthetic_selfcheck` | IM7/8552 | CAI + TAI | 12 | 1% | 0.0% (by construction) | 0.0% (not gated) | Tautological (empirical-tier output) |

These are the MAEs `validate_bvid_public.py` prints (rounded). Three
things the gate does not show:
- **The gated strength MAEs are in-sample.** `soutis_k_s` = 4.8 is the
  pooled optimum over NCAMP and Lovejoy & Scotti, so both gated MAEs are
  in-sample; the fair measure is the cross-dataset check under
  [Calibration of `soutis_k_s`](#calibration-of-soutis_k_s) (5.0–8.8%).
  Hasebe et al.'s strengths are not fitted, but the model is known to
  miss them and they are not gated.
- **The semi-analytical growth channel is not exercised by any
  dataset.** On every NCAMP and Lovejoy & Scotti coupon the Soutis bound
  is lower than the delamination growth stress (176–177 against
  285–287 MPa on NCAMP; 156 against 307–308 MPa and 224 against 408 MPa
  on Lovejoy & Scotti), and with the energy release rate of the 2D
  blister (the post-buckled elliptical film) it is lower on every Hasebe
  et al. coupon too, so the two tiers give identical numbers on all
  three. The growth channel is therefore unvalidated. With the
  straight-sided strip's energy release rate the growth stresses were
  217–219, 225–226 and 271 MPa, and growth governed 176 of the 343
  Hasebe coupons (94 of 157 cross-ply c8/c16/c24, 15 of 27 hard r0, 67
  of 132 quasi-isotropic, no soft r45), where it made the prediction
  worse: 28.6% against 23.9% empirical (both 19.1% on the other 167),
  24.0% over the set. Per layup the semi-analytical MAE went from 19.8%
  to 14.1% (c8), 33.0% to 31.1% (c16), 20.7% to 18.2% (c24), 21.6% to
  19.8% (q8), 26.5% to 24.2% (q16), 39.0% to 37.5% (q24) and 17.8% to
  15.4% (r0), the empirical numbers; r45 stayed 4.0%. (Before the
  clamped-ellipse sublaminate of #132 the growth stresses were 205–206,
  182 and 227–228 MPa, and growth governed 256 Hasebe coupons: 25.6%
  against 23.0%, 17.2% on the other 87.)
- **TAI is not validated.** The only measured tension data are one
  impact condition (`korbelin_m21_t800s_tai`, five coupons), and
  `tai_notch_fraction` is fitted to it; nothing out of sample tests the
  closed-form TAI model. fe3d scores 3.6% there (`--tier=fe3d`, not run
  in CI), also partly in-sample: its fibre-damage hole is the closed-form
  one.

The Lovejoy & Scotti strength MAE also rests on a cancellation between
the knockdown and the pristine strength. Against the measured pristine
strength its knockdown MAE is 13.6%, not 8.8%; see [Knockdown against the
measured pristine strength](#knockdown-against-the-measured-pristine-strength).

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
damage-zone factor is fitted to, see below). At about 20 s per coupon it
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
| 0.20 | −5.6% | +20.6% | +15.0% | 13.7% |

0.15 is kept, with the lowest mean error, 7.4%: 141.0 MPa (177.7
measured), 148.6 MPa (146.8) and 236.4 MPa (235.7). First-ply failure
sets all three up to 0.17; at 0.20 the hard laminate is set by its
delamination growth stress (270.9 MPa), which the factor does not change.
With the enclosing-rectangle sublaminate before issue #132 that growth
stress was 227.5 MPa and set the hard laminate from 0.15 up (−3.5% at
0.15 / 0.17 / 0.20, mean errors 8.4% / 9.0% / 9.9%; the factor was not
refitted, and 0.15 stays the lowest either way); before issue #166 fe3d
rescaled it by σ0 / FPF_undamaged and first-ply failure set it. No single factor fits
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
strengths as qualitative. fe3d takes about 20 s on the 24-ply NCAMP
coupon and about 45 s on the 40-ply hard laminate (5 mm mesh).

**Why fe3d's knockdown does not fall with damage area.** fe3d's CAI
knockdown is set by the strain concentration at the soft zone's flank,
where a 0° ply fails in fibre compression. That concentration depends on
the zone's width and aspect ratio, not its area. In a wide panel the same
zone gives the same knockdown from 700 to 5,269 mm². On the NCAMP coupon
the template grows self-similarly until its width reaches the 75 mm
impact window at about 4,200 mm². Beyond that, all further area goes
along the load direction, which lowers the flank concentration. fe3d on
the NCAMP coupon therefore gives a knockdown with a minimum at the cap:

| Damage area (mm²) | 1,400 | 2,783 | 3,500 | 4,200 | 5,269 | 6,023 |
|---|---|---|---|---|---|---|
| fe3d knockdown | 0.336 | 0.260 | 0.240 | **0.217** | 0.253 | 0.239 |

- **The coincidence:** the NCAMP coupon (2,783 mm²) and the Lovejoy &
  Scotti coupons (5,269–6,991 mm²) sit on either side of the minimum, so
  their near-equal fe3d knockdowns are a coincidence, not a material or
  laminate effect.
- **Without the window cap** (the zone bounded only by the coupon edges)
  the knockdown keeps falling: 0.260 / 0.188 / 0.133 at 2,783 / 5,269 /
  6,991 mm².
- **What is not the cause:**
  - The 2.5 mm averaging disk: radii of 0–20 mm change the 2,783 to
    5,269 mm² ratio only to 0.97–1.17, against 1.34 measured.
  - A characteristic-length (Whitney-Nuismer) hole: it is flat above
    about 1,400 mm² for d0 = 1–4 mm.
- **What the data cannot separate:** with the width pinned by the window,
  no elastic zone model makes the knockdown fall with area. The measured
  fall needs either damage wider than the window or damage severity that
  rises with energy (6.7 J/mm on the NCAMP coupon, 8.9 J/mm on the Lovejoy
  & Scotti QI). The two datasets move both together and cannot tell them
  apart.
- **What would decide it:** measured damage widths from large-area ASTM
  C-scans, and an energy series on a single laminate and coupon.
- **What is not adopted:**
  - Bounding the footprint by the coupon edges (with the factor refitted
    to 0.22) brings back the area dependence. It scores NCAMP +0.2%,
    Lovejoy & Scotti QI +7.0% and hard −14.2%, a 7.1% mean against 7.4%
    today. It is not adopted: it undoes the confinement that fixed the
    hard laminate, and no measured widths support it.
  - Softening only the plies bounding a delamination gives no size
    effect in 3D.

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
21.53% / 23.35% (23.45% semi-analytical once the sublaminate buckling
rectangle followed the ellipse orientation, #165, 23.97% with the
clamped-ellipse sublaminate, #132, and 21.53% with the 2D blister energy
release rate). The NCAMP damage area
is about 2,790 mm² either way.

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
21.5% (empirical and semi-analytical; 24.0% semi-analytical with the
straight-sided strip's energy release rate). The empirical MAE per
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

The fe3d tier gives 0.870 here (909.0 MPa against 884–1014 MPa, MAE
3.6%, 5 mm mesh). In tension its damage zone keeps its fibre-direction
stiffness and strength (delaminations and matrix cracks cut no fibres)
and loses only its matrix-dominated moduli, and the broken fibres are
the closed-form equivalent hole (1.5 mm across here, below the mesh),
applied as a Whitney-Nuismer notch on the stress of the elements under
it. The result is therefore partly in-sample: the hole size uses
`tai_notch_fraction`, fitted to these coupons; only the zone around it
(fitted on CAI) is out of sample, and it costs little here (0.870
against the closed form's 0.875). Without the hole fe3d gave 0.965 (MAE
9.0%): the zone's matrix cracking alone, nearly flat in energy. Before
the zone kept its fibres, it was softened in every direction in tension
too and failed as if its fibres were cut: 0.335 (MAE 62.2%). See
`docs/physics_models.md`, *Limitations*.

Other TAI sources looked at and not used: NASA TP-2981 (Nettles 1990)
tests eight-ply coupons punctured over a 10.3 mm hole with a 4.2 mm tup,
reports its results only as plots and uses no preset material; Lovejoy &
Scotti 2019 has no TAI tests; Kötter's 2021 TUHH thesis has fatigue after
impact at one energy only; Itakura et al. 2022 (JSME) is gas-gun
penetration; Poe's NASA reports cover thick filament-wound rocket motor cases;
Körbelin's 2022 TUHH thesis ("Damage tolerance of high-performance
composites", Band 41) was not located online.

**`gonzalez_as4_8552_cai`** is measured data, reported but not gated
(`"gate": false`) until the maintainer decides: E.V. González's PhD
thesis (Universitat de Girona, 2010; open on TDX,
[hdl.handle.net/10803/22834](http://hdl.handle.net/10803/22834)), which
tabulates per coupon the tests behind González et al., *Compos Sci
Technol* 71 (2011) 805–817 (ply clustering) and *Compos Struct* 94 (2012)
3364–3378 (simulation). Only the thesis was read in full; the two papers
were read as abstracts. Hexply AS4/8552 tape, 150 × 100 mm coupons
impacted on the ASTM D7136 125 × 75 mm window (four rubber-tipped
clamps, `impact_panel_*`) by a 16 mm hemispherical tup of 5 kg, then
tested in an ASTM D7137 fixture. 19 coupons on 5 laminates:
- L02 quasi-isotropic [(45/0/-45/90)4]s, 32 plies, 5.8 mm, at 38.6, 28.6
  and 19.3 J (4 coupons);
- L03 [(45₂/0₂/-45₂/90₂)2]s and L04 [45₄/0₄/-45₄/90₄]s, the same
  in-plane layup with plies clustered in pairs and fours (4 each);
- L06 cross-ply [(0/90)9]s, 36 plies, 6.8 mm (4);
- L05 quasi-isotropic [(45/0/-45/90)6]s, 48 plies, 8.8 mm, at 38.6 J
  (5 kg) and 58 J (6 and 15 kg) (3).

What is measured and what is not:
- Every value is transcribed from a table or the text; nothing is
  digitized. CAI loads: Tables 4.4 and 4.8, L05 from Table 4.13 and
  section 4.4.4. Laminates, cured ply and plate thickness: Table 3.12.
  C-scan projected areas: Table 4.3, for L02–L04 only (12 coupons).
  Dents: the maximum indentation printed over each contour plot
  (Figs. 4.17 and 4.32), carried on the 9 coupons where the plotted
  specimen is the CAI specimen.
- The strength is the failure load over the nominal 100 mm width times
  the laminate's cured thickness, as the thesis does; per-coupon
  dimensions are not reported. The thesis's Tables 4.9 and 4.13 print
  the result a factor 10 low (16.6 for 96 kN / (100 × 5.8 mm) =
  165.5 MPa).
- No undamaged compression strength is reported, and the L05 and L06
  C-scans are images without a scale. The L05 coupons at 38.6 J and the
  58 J pair are one value each, with the specimen not named.
- The 16-ply L01 (3.0 mm) is left out: the thesis reports global
  buckling in the CAI fixture, and its 197 MPa is near the plate's
  buckling stress.

The thesis's AS4/8552 card (Table 7.1, partly measured on the same
prepreg roll) differs from the `AS4/8552` preset (NCAMP), which is not
changed: E22 7.63 against 9.24 GPa, G12 4.36 against 4.83 GPa, ν12 0.35
against 0.30, Xt 2300 against 1872 MPa (the preset uses NCAMP's
cross-ply back-out), Xc 1531 against 1398 MPa, Yt 26 against 64 MPa (the
thesis applies in-situ factors to it), Yc 200 against 268 MPa and S12 78
against 92 MPa. E11 (128.0 against 127.3 GPa), ν23 (0.45), G_Ic (0.28),
G_IIc (0.79 N/mm) and the density agree.

Scores (MAE, signed bias, `soutis_k_s` = 4.8; nothing is fitted to this
dataset). "Measured area" puts the C-scan area in place of the predicted
one; "vs NCAMP pristine" is the knockdown error against NCAMP's measured
25/50/25 unnotched compression (560.7 MPa; the model gives 543.2 MPa):

| Laminate | n | h (mm) | Empirical MAE (bias) | Semi-analytical MAE | Predicted / C-scan area | Bias with measured area | Knockdown bias vs NCAMP pristine | One scale per laminate |
|---|---|---|---|---|---|---|---|---|
| L02 quasi-isotropic | 4 | 5.8 | 3.9% (−2.5%) | 3.9% | 1.25× | +5.3% | +0.6% | 1.00 |
| L03 clustered ×2 | 4 | 5.8 | 5.2% (+1.3%) | 5.2% | 0.80× | −5.9% | +4.6% | 0.95 |
| L04 clustered ×4 | 4 | 5.8 | 7.1% (+2.8%) | 7.1% | 0.52× | −19.3% | +6.1% | 0.98 |
| L06 cross-ply | 4 | 6.8 | 35.9% (+35.9%) | 35.9% | no C-scan | – | no pristine | 0.74 |
| L05 quasi-isotropic, 48 plies | 3 | 8.8 | 21.2% (−21.2%) | 21.2% | no C-scan | – | −18.7% | 1.25 |
| All | 19 | | 14.3% (max 41.2%) | 14.3% (max 41.2%) | | | | 3.8% MAE |

The replicate floor (coupons at the same laminate and energy) is 2.1%
over 10 coupons. Which channel governs: in the semi-analytical tier the
Soutis bound governs all 19 coupons, so the two tiers agree. (With the
straight-sided strip's energy release rate, before the 2D blister,
delamination growth governed 3: L03 at 28.6 and 19.3 J and L06 at
19.3 J, making L03 at 19.3 J worse, −17% against −8%, and the other two
3–4 points better; semi-analytical MAE 14.4%.) fe3d at the default 5 mm
mesh, on three coupons (not run in CI; first-ply failure governs all
three): L02-S01 119.7 against 160.3 MPa (−25%), L06-S01 110.1 against
166.2 MPa (−34%), L05 at 38.6 J 126.1 against 206 MPa (−39%), about 1 to
1.5 minutes each.

The per-layup pattern:
- **Clustering.** The 5.8 mm laminates score within 4–7%, but for the
  wrong reason. The predicted area is the same for all three
  (3,670 mm² at 38.6 J) while the C-scans grow from 3,022 to 3,955 and
  6,158 mm² as the plies are clustered, and the CAI strength does not
  change (the papers' conclusion). With the measured area, the
  knockdown of L04 is 19% too low. L04-S05 (19.3 J, +17.5%) is the one
  coupon the thesis singles out as low.
- **Cross-ply.** The model keeps too much strength: 35.9% high, one
  scale of 0.74 fits all four coupons (1.7%). The measured CAI of L06
  is the same as L02's (166–167 against 160–171 MPa at 38.6 J), while
  the model's pristine strength is 749 against 543 MPa and its
  knockdowns nearly equal. This agrees with Hasebe et al., whose
  cross-ply laminates need about 25% less strength than their
  quasi-isotropic ones relative to the model, and here the coupons are
  ASTM coupons, so it is not the JIS fixture. It is the opposite of
  Lovejoy & Scotti's hard 50/40/10 laminate, which needs about 26% more
  than its quasi-isotropic coupons. fe3d does not show the cross-ply
  excess: its damaged first-ply failure stress is 110 MPa on L06 and
  117 MPa on L02, so it puts the cross-ply 8% below the
  quasi-isotropic coupon, where the measurement is 4% above and the
  closed-form tiers 39% above.
- **Thickness.** The 48-ply laminate keeps more strength than the model
  gives (21% low). Its C-scans are not tabulated; the thesis says that
  at equal impact they are "quite similar" to L02's, and the predicted
  area at 38.6 J (3,404 mm²) is close to L02's measured 3,022 mm², so
  the miss is in the knockdown map, not the area.

Cross-validation of the layup terms with this dataset as a fourth set
(the CAI study's levers, refitted leave-one-dataset-out: fit on three
datasets with equal weight, score the fourth; validator scoring,
predicted area). The last column fits on two of the three ASTM datasets
(NCAMP, Lovejoy & Scotti, González) and scores the third, leaving
Hasebe's fixture-limited coupons out:

| Knockdown map | New constants | Held out: NCAMP / L&S / Hasebe / González | Mean | Exponent by fold (without NCAMP / L&S / Hasebe / González) | ASTM-only mean |
|---|---|---|---|---|---|
| Today (`k_s` 4.8, m 0.5) | 0 | 1.5 / 8.8 / 21.5 / 14.3 | 11.5 | – | 8.2 |
| Refit `k_s` | 0 | 2.2 / 8.7 / 21.6 / 14.2 | 11.7 | – | 8.4 |
| `k_s (Kt∞/3)^p` | 1 | 2.2 / 13.9 / 22.8 / 17.9 | 14.2 | +0.40 / +0.76 / +0.75 / −0.37 | 11.8 |
| `k_s (f0/0.25)^p`, 0° fraction | 1 | 2.2 / 21.0 / 22.5 / 15.7 | 15.4 | −0.15 / +0.56 / −0.15 / −0.12 | 14.0 |
| `k_s (Ex/Ex_QI)^p` | 1 | 4.2 / 24.2 / 22.3 / 15.3 | 16.5 | −0.31 / +1.16 / −0.24 / −0.18 | 15.1 |
| `k_s (1 + p f90)`, 90° interfaces | 1 | 1.5 / 12.3 / 22.1 / 16.2 | 13.0 | +0.44 / +0.51 / +0.48 / −0.10 | 13.7 |
| `k_s (h/4.5 mm)^q`, thickness | 1 | 3.5 / 10.1 / 23.0 / 14.8 | 12.8 | −0.09 / +0.22 / −0.20 / −0.14 | 14.0 |

No layup term survives: each one makes the held-out mean worse, by 1.3
to 5.0 points with four datasets and by 3.6 to 6.9 points on the ASTM
sets alone, and each makes the held-out González error worse too (14.8–17.9%
against 14.3%). Fitted on all four, `Kt∞` and the 90°-interface fraction
take González from 14.3% to 10.3% and 7.9% in-sample, but Lovejoy &
Scotti from 8.8% to 11.6–11.7%. The exponents still change sign between
folds: the `Kt∞` and 90°-interface exponents are positive whenever
González's cross-ply is in the fit and negative without it, and the 0°
fraction and Ex exponents are negative whenever Lovejoy & Scotti's hard
laminate is in the fit and positive without it. In every candidate the
hard and the cross-ply laminates sit on the same side of the
quasi-isotropic one (both have more 0° plies, a higher `Kt∞` and a higher
Ex), and the data want them moved in opposite directions, the hard
laminate up and the cross-ply down, so no monotone term in these
quantities can fit both. The result holds with
the measured areas (held-out mean 11.5% today against 11.9–16.7%), and
without L04-S05. Knockdown against NCAMP's measured pristine cannot test
the layup terms here, as the cross-ply has no measured pristine strength.
(Computed with the CAI accuracy study's lever scripts, which are not in
the repository; with the original three datasets they give its 10.6%
held-out mean.)

**`synthetic_selfcheck`** is tautological: each case is produced by running
the empirical tier at a known input, so its MAE is ~0% by construction.
The semi-analytical tier scores 0.0% on it (TAI is the same model, and
Soutis governs every compression case; before the clamped-ellipse
sublaminate of #132 the delamination growth stress was 14.5% lower at
5 J, 1.2% MAE), but CI does not run that tier on it. It
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

### Knockdown against the measured pristine strength

The gated MAE scores a strength. A case without `measured_pristine_MPa`
(NCAMP, Lovejoy & Scotti) is scored on the predicted knockdown times the
model's own pristine (unnotched) strength, so an error in that pristine
strength can hide an error in the knockdown, or the other way round. For
each case with a measured undamaged strength `P`, the validator therefore
also reports, without gating them:
- **knockdown error** = predicted knockdown / (measured residual / `P`) − 1;
- **pristine error** = model pristine strength / `P` − 1.

Each is printed as MAE and mean signed bias, per dataset and per
`layup_label` (a `knockdown (reported, not gated)` line after the MAE
line, one `layup` line per group), next to the gated strength MAE and its
bias. Where no case has a `P` it prints `n/a`. The per-case table gains
the columns `layup`, `predicted_pristine_MPa`, `measured_pristine_MPa`
(the `P` used), `measured_knockdown`, `knockdown_error_pct` and
`pristine_error_pct` (signed), and `knockdown_summary(df)` returns the
same numbers as a dict. `P` is `measured_pristine_MPa` when the case has
one (an undamaged coupon of the same test), else `measured_unnotched_MPa`.
That field is reported only, never scored: the measured unnotched strength
of the same material and layup family from another source. NCAMP and
Lovejoy & Scotti test no undamaged coupons, so their cases carry NCAMP's
RTD unnotched compression, normalized to the cured ply thickness (the
values of `tests/validation/test_unnotched_strength_ncamp.py`): AS4/8552
25/50/25 560.7 MPa (CAM-RP-2010-002), IM7/8552 25/50/25 600.2 MPa and
50/40/10 833.2 MPa (CAM-RP-2009-015). The Lovejoy & Scotti hard laminate
has the 50/40/10 ply fractions in another stacking, which the ply-discount
pristine strength does not depend on. The reference comes from other
coupons and another laboratory than the CAI tests, so the knockdown error
carries a few percent of that uncertainty.

Empirical tier (semi-analytical identical on every dataset):

| Dataset / layup | Cases | Strength MAE / bias (gated) | Knockdown MAE / bias vs measured pristine | Pristine error |
|---|---|---|---|---|
| NCAMP, 25/50/25 | 7 | 1.5% / +0.4% | 3.6% / +3.6% | −3.1% |
| Lovejoy & Scotti, all | 8 | 8.8% / +1.4% | 13.6% / −3.7% | MAE 6.1% / bias +6.1% |
| Lovejoy & Scotti, QI 25/50/25 | 4 | 10.2% / +10.2% | 10.0% / +10.0% | +0.2% |
| Lovejoy & Scotti, hard 50/40/10 | 4 | 7.4% / −7.4% | 17.3% / −17.3% | +11.9% |
| Hasebe et al. (scored on `P`) | 343 | 21.5% / −6.6% | the same | +101% (−7% r45 to +237% c8) |
| Körbelin et al., TAI (scored on `P`) | 5 | 3.5% / −1.4% | the same | −3.9% (proxy card) |

On Lovejoy & Scotti the 8.8% is a cancellation. The model puts the hard
laminate's pristine strength at 932.6 MPa, 11.9% above NCAMP's 833.2 MPa,
and its knockdown (0.240) 17.3% below the measured one (0.283–0.299
against NCAMP's pristine strength); the product is 7.4% low. The
quasi-isotropic coupons have the right pristine strength (+0.2%) and a
knockdown 10.0% high. Against the measured pristine strength the dataset's
knockdown MAE is 13.6%, and the knockdown map has a +10% / −17% split
between the two layups that the strength MAE shows only as +10% / −7%.
On NCAMP the two errors offset a little the other way: pristine −3.1%,
knockdown +3.6%, strength bias +0.4%. On a case scored on
`measured_pristine_MPa` the knockdown error is the signed strength error,
and the pristine error measures the test rather than the model: Hasebe
et al.'s undamaged coupons are fixture-limited (above), and only the soft
r45 laminate, which the fixture does not limit, has a pristine strength
near the model's.

The IM7/8552 card's `Xc` (1590 MPa, a generic value) is not NCAMP's, and
the pristine compression strength of every layup scales with it. The
NCAMP data do not give one `Xc` for both laminates: 1420 MPa puts the
50/40/10 laminate on 833.2 MPa but the 25/50/25 laminate 10.5% low, and
1587 MPa fits 25/50/25 only. The pristine error of the hard laminate is
therefore a layup error of the ply-discount strength as much as a
material constant.

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

### Compression: a cap on the 0° ply strain

In compression no ply cracks before the 0° fibres fail, so the model is a
cap on the 0° ply strain. Every laminate fails at the card's `Xc / E11`
(1.098% for AS4/8552 and 0.964% for IM7/8552; within 0.6% in the NCAMP
layups and the [90/0/90]ns cross-ply). Its strength is therefore the
CLT modulus times one strain.

NCAMP's measured unnotched compression moduli split each error into a
stiffness part and a failure-strain part. Measured failure strain is
strength over the measured modulus; the brackets give it relative to the
material's 25/50/25:

| Material | Cross-ply (UNC0) | 25/50/25 | 50/40/10 | 10/80/10 | Card `Xc / E11` | NCAMP `F1cu / E1c` |
|---|---|---|---|---|---|---|
| AS4/8552 | 1.288% (1.02) | 1.265% (1) | 1.361% (1.08) | 1.452% (1.15) | 1.098% | 1.254% |
| IM7/8552 | 1.265% (1.14) | 1.108% (1) | 1.015% (0.92) | 1.356% (1.22) | 0.964% | 1.242% |

- **The stiffness part is not the problem.** With NCAMP's compression
  constants (`E1c` 111.5 and 138.2 GPa, where the cards use the tension
  `E11`), CLT reproduces the measured moduli within 2% on the
  multidirectional laminates and within 4% on the cross-plies. The cards'
  moduli are 9–19% stiff and their failure strains low, and the two
  roughly cancel.
- **The failure strain depends on layup and material together.**
  - IM7/8552's 50/40/10 fails at 0.92× its 25/50/25 strain, the lowest of
    any IM7/8552 laminate. AS4/8552's fails at 1.08×, with the same
    stacking, resin, lab and method. That is about a 17% gap, against
    standard errors of the means of 0.7–2.3%, and NCAMP lists no outliers
    for either.
  - The 10/80/10 laminates fail 15–22% above their 25/50/25 strain. That
    is the opposite of what kinking under the transverse tension their
    Poisson mismatch puts in the 0° plies predicts.
- **So no layup-only correction can fix it.** None of the following
  improves on the current model's 10.2% mean |log error| over the six
  laminates, with literature constants and no per-layup fitting:
  - ply-discount factor, first-fibre against last-ply failure, in-situ
    matrix strengths, or a Budiansky-Fleck shear term: unchanged, because
    the 0° plies carry no shear;
  - LaRC03 kinking: 14.1%;
  - LaRC03 with in-situ shear strength: 15.8%;
  - NCAMP's full cards: 10.8%;
  - one-constant laws in the ±45 fraction, the 0° fraction or the Poisson
    mismatch: 8.4% at best held out, against 8.6% unchanged, and each
    leaves IM7/8552 50/40/10 13–28% high. The Poisson-mismatch law also
    puts the cross-plies 16–20% low.

**What a pristine correction would do to CAI.** The knockdown error
against the measured pristine strength does not depend on the model's
pristine strength. Lovejoy & Scotti stays at 13.6% (QI +10.0%, hard
−17.3%), so correcting the pristine strength alone only moves the gated
strength MAE:
- with the measured pristine strength in place of the model's, NCAMP goes
  from 1.5% to 3.6% and Lovejoy & Scotti from 8.8% to 13.6%;
- with NCAMP-consistent IM7/8552 constants (`E1c` 138.2 GPa, `F1cu` 1716
  MPa), the IM7/8552 pristine strengths become 10% (QI) and 21% (hard)
  high, and the synthetic self-check fails its gate.

The fix belongs in the knockdown map. It would need IM7/8552's hard
laminate about a third higher relative to its quasi-isotropic laminate
(measured knockdown against NCAMP's pristine 0.290 hard, 0.237 QI; the
model gives 0.240 and 0.260). The leave-one-dataset-out cross-validation of the layup terms (in the
González dataset's section) rejects such a layup term, because the
cross-plies want the opposite move. `Xc` stays at the card's 1590 MPa.

**Data that would decide it:**
- Undamaged D7137 coupons of Lovejoy & Scotti's two layups, from the same
  panels.
- A second-source IM7/8552 50/40/10 unnotched compression strength with
  its modulus.
- The failure modes of NCAMP's 50/40/10 laminates.

Sources: NCAMP NCP-RP-2009-028 Rev B (IM7/8552) and NCP-RP-2010-008 Rev D
(AS4/8552) statistical analysis reports, and the NCP-RP-2013-015 IM7/8552
shared database, whose 11 units reproduce the cross-ply strength.

## Next datasets

From the literature survey of open CAI/TAI data, in priority order. Each
needs the inputs listed before it can be scored.

| Dataset | Material | Loading | Needs |
|---|---|---|---|
| Sánchez-Sáez et al. 2008 (UC3M open manuscript) | AS4/3501-6 | CAI | Thin, stability-limited 78 mm coupons |
| NASA TP-3102, drop-weight point | T800/3900-2 | CAI | Single point |
| Falcó et al. 2014, *Compos Struct* 116, 94–103 (Girona, tow-drop gaps in variable-stiffness panels) | AS4/8552? (only the abstract was seen) | CAI | The full text (the Girona repository blocks automated downloads); only its straight-fibre baseline laminates can be modelled |

González et al. 2011 is now `gonzalez_as4_8552_cai`, from the thesis
that tabulates its coupons. A Girona dataset with a hard (50/40/10-type)
laminate and a measured undamaged strength would be the most useful next
step: it is the combination that separates the cross-ply and hard
laminate corrections above.

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
`measured_unnotched_MPa` (the measured unnotched strength of the same
material and layup family from another source, for a dataset with no
undamaged coupons; reported against the knockdown, not scored),
`layup_label` (the group the knockdown errors are broken down by),
`impact_panel_Lx_mm`, `impact_panel_Ly_mm` and `impact_boundary` (the
window the plate was impacted on, when it is not the test panel), and
`note`. Give a new dataset a measured pristine strength (either field)
where the source has one, so its knockdown is checked as well as its
strength.

## Running

```bash
python validation/validate_bvid_public.py              # empirical tier, all datasets, print table
python validation/validate_bvid_public.py --gate       # exit non-zero if a gated MAE exceeds 1.25 * target
python validation/validate_bvid_public.py --dataset synthetic_selfcheck
python validation/validate_bvid_public.py --tier semi_analytical --measured   # every dataset with "kind": "measured"
python validation/validate_bvid_public.py --tier fe3d --dataset synthetic_selfcheck --max-cases 3
```

`--tier` is `empirical` (default), `semi_analytical` or `fe3d`. The gate
multiplier is 1.25 for the closed-form tiers and 2.0 for fe3d.
`--dataset` takes one or more file stems; a name that matches no file is
an error (exit 1). `--max-cases N` runs only the first N cases of each
dataset.

## Adding a new dataset

1. Create `datasets/<name>.json` with top-level `name`, `kind`
   (`"measured"` for test data from a cited source, `"synthetic"` for data
   generated by the model; required, the validator refuses a file without
   it), `reference` (the citation, tables used, and every assumed value),
   `target_mae_pct` (defaults to 100 if omitted) and `cases` (fields above;
   they map onto the `DatasetCase` dataclass in `validate_bvid_public.py`).
   The material must be a preset in `MATERIAL_LIBRARY`; add one first if
   needed (cited constants only).
2. The validator auto-discovers any `*.json` file in `datasets/`, so the
   CI empirical gate picks it up, and a `"kind": "measured"` file is also
   picked up by `--measured` (the CI semi-analytical gate). Run both
   closed-form tiers on it and record the MAE in the table above.
3. If the model is known to miss it, set `"gate": false` and a
   `gate_note` saying why, rather than loosening `target_mae_pct`.
4. Do not tune calibration constants (`soutis_k_s`, ...) on the new
   dataset alone. If it motivates a recalibration, fit on part of the data
   and report the error on data held out from the fit (the cross-dataset
   check above), and rerun every other dataset.
5. There is no list to edit: `tests/test_validation_smoke.py` reads the
   same files and runs the closed-form gates on every gated measured
   dataset.
6. Commit the dataset, any recalibration, and the updated numbers in a
   single PR.
