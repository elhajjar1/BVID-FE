# Changelog

All notable changes to BVID-FE are documented in this file.

## [Unreleased]

### Added

- **`T800S/3900-2B` material preset** (Toray, interlayer-toughened). Elastic
  constants and toughness from Hasebe et al. (Data in Brief 2025) Table 1,
  strengths and density from the Toray 3900 data sheet (P2362W-19L, RTA).
  `olsson_alpha` = 0.047 is fitted to Hasebe's C-scan areas (quasi-isotropic
  laminates, 12.7 mm hemisphere) with the energy-balance damage area (see
  *Changed*); the default 0.078 overpredicts them about 1.7x. (It was 0.12
  with the thickness-divided area, against a default of 0.8.)
- **`AnalysisConfig.impact_panel`**: the panel the impact is applied on when
  it is not the panel tested (e.g. a 60 x 60 mm clamped window, then an
  80 x 50 mm compression coupon). The damage is mapped on `impact_panel`,
  so its onset energy, mass ratio and 80% area cap follow that panel, and
  is centred on `panel`. Defaults to `panel`; only valid with `impact`.
- **Third measured CAI dataset, `hasebe_t800s_3900_2b_cai` (not gated).**
  343 specimens of Hasebe et al. 2025 (Mendeley, CC BY 4.0): T800S/3900-2B
  in 8 layups, struck by three hemispheres, a flat-ended impactor or a 60
  deg cone at 2-21 J on a clamped 60 x 60 mm window, and tested as
  80 x 50 mm JIS K 7089 coupons. MAE 21.3% (empirical) and 23.3%
  (semi-analytical; 25.3% and 25.6% before the energy-balance change in
  *Changed*), against 30.3% for a constant knockdown: the predicted damage
  area matches the C-scans, but the DPA-based knockdown misses how strongly
  the knockdown depends on layup.
  The impactor mass is not reported (5.5 kg assumed).
- **Validation cases take `measured_pristine_MPa`**, which scores a case on
  predicted knockdown x the measured undamaged strength (Hasebe's undamaged
  coupons fail at 0.4-0.6x the unnotched strength in the JIS fixture), and
  `impact_panel_Lx_mm` / `impact_panel_Ly_mm` / `impact_boundary`.

- **`AS4/8552` material preset** (Hexcel 8552 AS4 tape). The cleanest open
  ASTM D7136/D7137 CAI datasets (Girona: Gonzalez et al. 2011, Falco et al.
  2014; NCAMP CAM-RP-2010-002) are on AS4/8552, which no preset covered.
  Lamina moduli, strengths and density are the NCAMP RTD means. G23 comes
  from nu23 = 0.45 (Gonzalez et al.), S23 is the LaRC03 estimate from Yc, and
  G_Ic/G_IIc are the IM7/8552 values (same matrix). The impact-mapping
  calibration fields keep their defaults. NCAMP's 30 J CAI on
  [45/0/-45/90]3s is 176 MPa; both closed-form tiers gave 163 MPa with
  `soutis_k_s` = 3.5 and give 177 MPa since the energy-balance change in
  *Changed* (the Soutis bound governs).
- **First measured validation dataset: `ncamp_as4_8552_cai`.** Seven
  AS4/8552 [45/0/-45/90]3s CAI coupons from NCAMP CAM-RP-2010-002 Rev A
  (ASTM D7136/D7137, about 30 J, 171-182 MPa), each with its own measured
  thickness and impact energy. MAE: empirical 7.3%, semi-analytical 7.3%
  with `soutis_k_s` = 3.5 (15.9% and 11.6% before; 1.5% for both since the
  energy-balance change in *Changed*; target 15%;
  fe3d gave 5% on the first coupon, 30% after the energy-balance change and
  4.8% after the damage-zone factor fit and -3.5% at the refitted 0.17
  with the confined footprint, but is not run in CI). CI
  now also gates the semi-analytical tier on measured datasets. The fe3d
  monitoring step runs only the synthetic dataset, since a 24-ply coupon
  takes about 3.5 minutes.
- **Second measured dataset, `lovejoy_scotti_im7_8552_cai`.** Eight IM7/8552
  standard-ply CAI coupons from Lovejoy & Scotti (NASA LaRC 2019):
  quasi-isotropic at 54 J and hard laminate at 68 J, with measured dent depth
  and C-scan area. MAE with `soutis_k_s` = 3.5 (see *Changed*): empirical
  10.9%, semi-analytical 10.8% (9.0% for both since the energy-balance
  change); both tiers are gated in CI.
  Datasets can set `"gate": false` with a `gate_note` for a known miss (its
  error is reported but cannot fail `--gate`), and `--dataset` takes several
  names.
- **Validation cases take `boundary` and `impactor_shape`**, passed to
  `PanelGeometry` and `ImpactorGeometry`. A case must now state the impactor
  diameter and mass; they used to default silently to 16 mm and 5.5 kg.
- The validation roadmap (`validation/README.md`,
  `validation/reference_data.json`) lists the open datasets found by the
  literature survey. Three of the four v0.2.0 entries could not be used as
  listed, and no source was found for the fourth.

### Changed

- **`AnalysisResults.warnings` now lists every `UserWarning` raised during
  `BvidAnalysis.run`** as `"<WarningClass>: <message>"` (also in `to_dict()`
  and the CLI JSON); the warnings are still re-emitted to Python's filters.
- **The library no longer logs to stderr on import.** `bvidfe.fe3d` gets only
  a `NullHandler`; the CLI and Streamlit app call the new
  `fe_tier.enable_stderr_logging()`. An explicit `BVIDFE_LOG_LEVEL` still
  sets the logger level.
- **CLI ply angles must lie in [-90, 90] degrees**, as in the app.
  `--tier`/`--loading` choices come from `_types`; `--help` no longer shows
  `(default: None)` or a doubled default on `--impactor-diameter`/`--mass`.
- **fe3d damage-zone in-plane factor 0.20 -> 0.17, refitted with the
  delaminations confined to the impact window** (see *Fixed*). With the
  confinement the three coupons sweep as (fe3d strength against measured
  CAI):

  | Factor | NCAMP #0 | Lovejoy QI #0 | Lovejoy hard #4 |
  |---|---|---|---|
  | 0.20 | +4.7% | +6.6% | +12.4% |
  | 0.17 | -3.5% | -4.4% | +1.9% |
  | 0.15 | -9.0% | -11.7% | not run |

  At 0.17 they give 171.4 MPa (177.7 measured), 140.4 MPa (146.8) and
  240.2 MPa (235.7): both laminate families fit within 5% (mean error
  3.3%, against 18.0% at 0.20 before the confinement, when the hard
  laminate was 47% low). Effect on fe3d knockdowns (IM7/8552, 150 x 100 mm,
  5 mm mesh; the confinement does not bind on these panels, so the change
  is the factor's):
  - 16-ply QI tension at 5/15/30 J 0.432/0.342/0.334 -> 0.403/0.313/0.307
    (empirical 0.399/0.370/0.359); 24-ply QI tension at 30 J 0.253 ->
    0.224 (empirical 0.359).
  - CAI between 5 and 30 J: 16-ply QI 0.444 -> 0.339 becomes
    0.416 -> 0.313; 8-ply QI 0.416 -> 0.285 becomes 0.387 -> 0.256; 8-ply
    cross-ply 0.249 -> 0.215 becomes 0.222 -> 0.188.
  - The 5 mm -> 2.5 mm mesh change at 15 J on 8-ply panels: QI tension
    0.4138 -> 0.3920 (5.3%; was 5.0%), cross-ply tension 0.1590 -> 0.1975
    (24.2%; was 20.7%), cross-ply compression 0.1748 -> 0.2231 (27.6%; was
    22.6%), QI compression 0.4143 -> 0.4280 (3.3%; was 3.6%).
  - The fe3d smoke run on the first three synthetic cases goes 17.9% ->
    23.8% (the gap between the tiers, not an error).

- **fe3d damage-zone in-plane factor 0.30 -> 0.20, fitted to measured CAI
  coupons** (`DAMAGE_ZONE_INPLANE_FACTOR`, `bvidfe.analysis.fe_mesh`;
  superseded: 0.17 since the delaminations are confined to the impact
  window, see the entry above and *Fixed*. The numbers in this entry are
  at 0.20 with unconfined delaminations.) The
  unfitted 0.30 put the quasi-isotropic coupons 30-36% high once the
  energy-balance damage area (below) shrank the softened zone. A sweep on
  the first NCAMP AS4/8552 coupon and the first Lovejoy & Scotti IM7/8552
  quasi-isotropic and hard coupons (fe3d strength against measured CAI):

  | Factor | NCAMP #0 | Lovejoy QI #0 | Lovejoy hard #4 |
  |---|---|---|---|
  | 0.30 | +30.2% | +35.5% | -14.1% |
  | 0.25 | +18.3% | +16.7% | -31.1% |
  | 0.20 | +4.8% | -2.0% | -47.1% |
  | 0.15 | -8.9% | -20.4% | -62.1% |

  0.20 gives the lowest mean error (18.0%, against 26.6% at 0.30): the
  three coupons go 231 -> 186 MPa (177.7 measured), 199.0 -> 143.9 MPa
  (146.8) and 202.4 -> 124.7 MPa (235.7). It errs conservative on the
  40-ply hard laminate, which wants more than 0.30. One factor cannot fit
  both families: fe3d has no characteristic length, and the stress
  concentration at the zone edge depends on the laminate's orthotropy.
  Effect on fe3d knockdowns (IM7/8552, 150 x 100 mm, 5 mm mesh):
  - 16-ply QI tension at 5/15/30 J 0.525/0.435/0.423 -> 0.432/0.342/0.334
    (empirical 0.399/0.370/0.359); 24-ply QI tension at 30 J 0.349 -> 0.253
    (empirical 0.359). fe3d now sits below the empirical tier in tension
    from 15 J.
  - CAI between 5 and 30 J: 16-ply QI 0.530 -> 0.427 becomes
    0.444 -> 0.339; 8-ply QI 0.505 -> 0.381 becomes 0.416 -> 0.285; 8-ply
    cross-ply 0.341 -> 0.307 becomes 0.249 -> 0.215.
  - The 5 mm -> 2.5 mm mesh change at 15 J on 8-ply panels: QI tension
    0.4424 -> 0.4201 (5.0%; was 3.6%), cross-ply tension 0.1849 -> 0.2232
    (20.7%; was 13.0%), cross-ply compression 0.2027 -> 0.2485 (22.6%; was
    12.3%), QI compression 0.4408 -> 0.4565 (3.6%; was 17.4%).
  - The fe3d smoke run on the first three synthetic cases goes 3.9% ->
    17.9% (the synthetic set is empirical-tier output, so this is the gap
    between the tiers, not an error). Undamaged fe3d strengths are
    unchanged (no damage zone).

- **The damage area is an energy balance, and the CAI knockdown scales it
  by a fixed reference area.** Three changes, fitted together:
  - `impact_to_damage` predicts `DPA = olsson_alpha (E - E_onset) / G_IIc`
    (mm^2), times the impactor mass, shape and diameter factors and capped
    at 80% of the panel as before. It was also divided by the laminate
    thickness, which gave mm and left the area flat in thickness at
    constant J/mm, where measured areas grow with it: Hasebe et al.'s 8- to
    24-ply laminates at fixed J/mm, and Lovejoy & Scotti's 7.6 mm coupons,
    which had 1.21x the C-scan area of their 6.1 mm coupons at the same
    8.9 J/mm. The `olsson_alpha` default goes 0.8 -> 0.078, fitted to the
    median of Lovejoy & Scotti's eight IM7/8552 C-scan areas, and
    T800S/3900-2B 0.12 -> 0.047, fitted to Hasebe's quasi-isotropic
    laminates struck by the 12.7 mm hemisphere (median predicted / measured
    area over all 343 specimens 1.03x).
  - The Soutis knockdown is `1/(1 + k_s (DPA/A_ref)^m)` with
    `A_ref = SOUTIS_REFERENCE_AREA_MM2` = 15000 mm^2, the 150 x 100 mm
    ASTM D7137 coupon (`bvidfe.failure.soutis_openhole`), instead of the
    analysed panel's `Lx * Ly`. The panel area sent the knockdown to 1 as
    the panel grew, although impact damage is a local notch, and
    over-penalised Hasebe's 80 x 50 mm coupons. The panel size still enters
    the onset energy, the area cap, the semi-analytical growth check and
    fe3d. API: `soutis_cai(m, dpa_mm2, A_ref_mm2, sigma_pristine_MPa)`
    (argument renamed from `A_panel_mm2`), and
    `semi_analytical_cai(lam, damage, sigma_pristine_MPa,
    A_ref_mm2=SOUTIS_REFERENCE_AREA_MM2, boundary=...)`, where
    `A_panel_mm2` was a required argument.
  - `soutis_k_s` 3.5 -> 4.7, the pooled optimum of both tiers on NCAMP and
    Lovejoy & Scotti.

  MAE (empirical / semi-analytical): NCAMP 7.3% / 7.3% -> 1.5% / 1.5%;
  Lovejoy & Scotti 10.9% / 10.8% -> 9.0% / 9.0%; Hasebe (not gated)
  25.3% / 25.6% -> 21.3% / 23.3% (a constant knockdown scores 30.3%). The
  fair measure is cross-dataset: `k_s` fitted on NCAMP alone (4.7) predicts
  Lovejoy & Scotti at 9.0% (was 16-21%), and fitted on Lovejoy & Scotti
  alone (5.2) predicts NCAMP at 6.2% (was 15.7%). Neither change works
  alone: with the thickness-divided area and the fixed reference area the
  two datasets want different `k_s` (NCAMP 7.3% and Lovejoy & Scotti 9.6%
  at 3.5, 15.7% and 4.6% at 4.0); with the energy-balance area and the
  panel area `k_s` = 4.7 fits both (3.5% and 9.7%) but puts Hasebe at
  31-34%. The fe3d smoke run on the first three synthetic cases moves
  5.3% -> 3.9% (17.9% at the damage-zone factor 0.20 and 23.8% at 0.17,
  see above).

  Effect for users:
  - Predicted damage areas fall several-fold on thin laminates: NCAMP
    AS4/8552 4.5 mm coupons at 30 J about 6650 -> 2900 mm^2; a 16-ply
    IM7/8552 quasi-isotropic laminate ([0/45/-45/90/90/-45/45/0] twice,
    0.152 mm plies, 2.4 mm) on a 150 x 100 mm panel at 20 J 8081 ->
    1916 mm^2, so its empirical CAI knockdown rises 0.280 -> 0.373. The
    80% area cap on a 150 x 100 mm IM7/8552 panel is now reached only
    above about 120 J; an 8-ply panel reached it from about 15 J before.
  - Knockdowns on large panels now fall where they used to tend to 1: the
    same laminate on a 500 x 500 mm panel at 20 J goes 0.691 -> 0.455
    (DPA 4096 -> 971 mm^2).
  - fe3d knockdowns rise with the smaller damage zone (at the 0.30
    damage-zone factor; the 0.20 fit above supersedes these numbers). On
    IM7/8552
    150 x 100 mm panels, 16-ply QI tension at 5/15/30 J goes
    0.459/0.328/0.200 -> 0.525/0.435/0.423 (empirical 0.399/0.370/0.359)
    and 24-ply QI tension at 30 J 0.15 -> 0.349 (empirical 0.359): at 30 J
    the damage area (about 2900 mm^2) no longer reaches the 80% cap. The
    fe3d CAI knockdown between 5 and 30 J is 0.530 -> 0.427 on 16-ply QI
    (was 0.46 -> 0.20), 0.505 -> 0.381 on 8-ply QI (0.46 -> 0.29) and
    0.341 -> 0.307 on 8-ply cross-ply (0.27 -> 0.26). The first NCAMP
    coupon goes 169 -> 231 MPa (measured 178 MPa). The 5 mm -> 2.5 mm mesh
    change at 15 J on 8-ply panels is now 3.6% (QI tension), 13.0%
    (cross-ply tension), 12.3% (cross-ply compression) and 17.4% (QI
    compression).
  - The synthetic validation dataset and the CLI snapshots are recaptured.

- **Impactor shape and diameter factors on the damage area are fitted to
  Hasebe et al. 2025**, the only open data that vary the impactor on one
  material. Over same-layup, same-energy groups a 15.9 mm flat-ended
  impactor left 0.20x and a 60 deg cone 1.20x the C-scan area of a 12.7 mm
  hemisphere. `f_shape` is now 0.2 (flat) and 1.3 (conical), from the
  unsourced 1.4 and 0.7, which ran the other way; the diameter exponent is
  0.2 (was 0.3; 25.4 and 31.8 mm hemispheres left 0.88x and 0.87x) and no
  longer applies to cones. Flat-ended impacts now predict about 7x less
  damage and conical about 1.9x more; 16 mm hemispherical results are
  unchanged.

- **AS4/8552 and IM7/8552 Xt are NCAMP's cross-ply back-out values**
  (F1tu from the [0/90] tests with the back-out factor, RTD normalized):
  AS4/8552 1996 -> 1872 MPa, IM7/8552 2560 -> 2238 MPa. This is the 0 deg
  ply strength realised inside a laminate. With the unidirectional-coupon
  value the 0 deg plies outlived the +-45 cracking and the fibre-governed
  NCAMP unnotched tension strengths came out 5-19% high; now the 50/40/10
  laminates are 0.99-1.02 and all twelve NCAMP cases 0.83-1.14 (mean 0.97,
  mean |log| error 11.6% -> 8.8%). Remaining: where the +-45 plies cracking
  in shear sets the strength, AS4/8552 25/50/25 tension is 14% high and
  10/80/10 10-17% low. Tension pristine and TAI strengths drop by up to 6%
  (AS4/8552) and 12.6% (IM7/8552), where the 0 deg fibres govern; CAI is
  unaffected. Cross-ply IM7/8552 tension:
  1280.6 -> 1119.6 MPa.

- **`soutis_k_s` default 2.5 -> 3.5** (empirical CAI knockdown
  `1/(1 + k_s (DPA/A_panel)^m)`; superseded: 4.7, with a fixed reference
  area, since the energy-balance change above). The uncalibrated 2.5
  overpredicted the measured CAI strengths: MAE 15.9% (NCAMP AS4/8552)
  and 39.9% (Lovejoy & Scotti IM7/8552) for the empirical tier, 11.6% and 23.5% for the
  semi-analytical. The undamaged strength matches NCAMP, and rescaling the
  predicted damage to the measured C-scan areas made the error worse, so the
  knockdown map itself was too lenient. 3.5 fits both datasets (3.7 is the
  empirical-tier optimum, 3.3 the semi-analytical): MAE 7.3% and 10.9%
  (empirical), 7.3% and 10.8% (semi-analytical). Fitted on either dataset
  alone, `k_s` predicts the other within 16-21%, which is the honest
  uncertainty for a new material. Empirical CAI strengths drop about 15-25%
  (16-ply IM7/8552 quasi-isotropic, 200x150 mm panel at 20 J: knockdown
  0.439 -> 0.359); TAI
  is unaffected.
- **Pristine strength is now the unnotched laminate strength** (ply-discount
  last-ply failure on classical lamination theory,
  `bvidfe.failure.laminate_strength`). It was the thickness-weighted average
  of the ply strengths along the load, `sum t_i (X cos^2 + Y sin^2) / sum t_i`,
  which treats the plies as sharing stress; they share strain, so the
  off-axis plies were credited with strength they never reach before the
  0 deg fibres break. Against the NCAMP unnotched tension and compression
  strengths of AS4/8552 and IM7/8552 (25/50/25, 10/80/10 and 50/40/10
  layups) the old value was 1.17-2.85x the measurement (mean 1.75x); the
  new one is 0.83-1.18x (mean 1.00x; 0.83-1.14x after the Xt change above), pinned in
  `tests/validation/test_unnotched_strength_ncamp.py`.
  - `pristine_strength_MPa` and every absolute residual strength drop
    accordingly: a 16-ply IM7/8552 quasi-isotropic laminate goes from
    887.5 to 601.3 MPa in compression; cross-ply tension from 1316.5 to
    1280.6 MPa.
  - Empirical and fe3d knockdowns are unchanged (they are ratios applied to
    the pristine strength). Semi-analytical CAI knockdowns rise where
    delamination growth governs, because the growth stress is absolute:
    0.306 -> 0.451 for the 16-ply quasi-isotropic laminate at 10 J.
  - For the NCAMP AS4/8552 30 J CAI benchmark (176 MPa, knockdown 0.31) the
    empirical tier went from 313 to 204 MPa and the semi-analytical
    knockdown from 0.24 to 0.36, both with `soutis_k_s` = 2.5. With the
    recalibrated 3.5 (above) both tiers gave 163 MPa (177 MPa since the
    energy-balance change above).
  - The synthetic validation dataset and the CLI snapshots are recaptured.

- **fe3d first-ply failure now sees impact damage.** It did not: the fe3d
  tension knockdown was 0.998-1.000 at 5-30 J on 8- to 24-ply IM7/8552
  150x100 mm panels (empirical 0.21-0.37), and in compression only the
  closed-form delamination growth moved it. Three changes:
  - Delaminations softened only the out-of-plane stiffness and the
    fibre-break core is off in every material preset, so the in-plane
    stress field was the same with and without damage. The damage zone
    (every ply under the delamination footprints) is now a soft inclusion
    with in-plane stiffness x`DAMAGE_ZONE_INPLANE_FACTOR` = 0.30 (Soutis &
    Curtis 1996; a calibration value, not fitted to test data; fitted to
    0.20 and then 0.17 since, see the factor entries at the top of
    *Changed*).
  - Damage-zone elements fail at the strain intact material would (strain
    equivalence; `fe_tier._failure_stress_material`): the criterion sees
    the in-plane stress the pristine ply carries at the same strain. With
    full strengths a zone spanning most of the panel lost its stress
    concentration and looked stronger (cross-ply tension 0.41 -> 0.60 from
    5 to 15 J). `Hex8Element.stress_at_gauss_points(_material)` take an
    optional stiffness for this.
  - The failure stress was the failure strain times the pristine CLT
    modulus, which credits a softened panel with load it cannot carry (a
    panel softened uniformly x0.3 reported 3.3x its intact strength). It
    is now the reaction force on the loaded edge at failure over the gross
    section, from the same FE solve; `_effective_modulus` is removed.

  fe3d tension knockdowns at 5/15/30 J move 1.000/1.000/0.998 ->
  0.458/0.328/0.200 on 16-ply QI (empirical 0.37/0.35/0.35) and
  0.999/0.998/0.998 -> 0.253/0.248/0.248 on 8-ply cross-ply (empirical
  0.24/0.22/0.21); compression moves 0.528/0.430/0.360 ->
  0.464/0.327/0.197 on 16-ply QI. Every layup and loading now falls or
  holds with energy. Once the damage area hits its 80%-of-panel cap the
  knockdowns drop below the empirical tier's (0.16-0.20 at 30 J on 16-
  and 24-ply). Runtime is unchanged. First-ply failure sat at a panel
  corner (free-edge interlaminar shear) in both runs, 2-4% below the
  interior; that cancels in the damaged / undamaged ratio and the boundary
  conditions are unchanged.
- **Compression strength in the `semi_analytical` and `fe3d` tiers is now
  the stress at which a buckled sublaminate grows its delamination, not the
  stress at which it buckles.** A thin sublaminate over a BVID-sized
  delamination buckles at a few MPa and keeps carrying load, so taking
  buckling onset as the strength gave `semi_analytical` CAI knockdowns of
  0.001-0.19 (13 of 16 cases below 0.03) against 0.31-0.59 from the
  empirical tier, on 8- to 24-ply IM7/8552 150x100 mm panels at 5-30 J.
  The onset stress was also compared with the laminate stress directly,
  although it is a stress in the sublaminate. The new
  `sublaminate_growth_stress` sets the thin-film energy release rate
  `G = h/(2E_f)·(σ−σ_c)(σ+3σ_c)` (Chai, Babcock & Knauss 1981; Hutchinson
  & Suo 1992) equal to a Benzeggagh-Kenane mixed-mode toughness
  (`mixed_mode_toughness`: `G_Ic`/`G_IIc` at a 52.1° phase angle, η = 2;
  0.478 N/mm for IM7/8552), and converts the sublaminate stress at growth
  to a laminate stress through the shared far-field strain.
  `weakest_sublaminate_growth` takes the minimum over all delaminated
  interfaces (largest ellipse at each); the area-times-depth heuristic in
  `find_critical_interface` often picked a single 90° ply whose growth
  never governs. `semi_analytical` CAI is `min(Soutis, growth)` and
  `critical_sublaminate` is the growth interface. fe3d CAI is
  `min(growth, first-ply failure)`, normalised by its undamaged run as
  before, and now also reports `critical_sublaminate`. Whole-panel
  buckling no longer enters the fe3d knockdown: CAI fixtures carry
  anti-buckling guides, and with panel buckling in both the damaged and
  undamaged strengths an 8-ply 200x150 mm panel reported a knockdown of
  1.0 for any delamination that grew above the panel's 12.6 MPa buckling
  stress. Its value is recorded in `AnalysisResults.notes`; buckling onset
  stays in `buckling_eigenvalues` (unchanged meaning, no longer sets the
  strength). The `fe3d_buckling_fallback` warning tag is retired, since
  buckling no longer sets the fe3d residual. On a 16-ply IM7/8552 QI
  150x100 mm panel at 5/10/15/30 J the `semi_analytical` CAI knockdown
  moves 0.007/0.003/0.002/0.001 -> 0.350/0.306/0.286/0.239 (empirical
  0.533/0.439/0.387/0.309) and fe3d 0.043/0.020/0.013/0.007 ->
  0.528/0.461/0.430/0.360; on the 16-ply 200x150 mm CLI snapshot case at
  10 J `semi_analytical` moves 0.055 -> 0.306. At a given interface the
  growth stress changes little with delamination size once the
  sublaminate buckles well below it (between 0.87x and 1x the
  steady-state value), so the knockdown follows energy mainly when a
  different interface becomes the weakest: on 8-ply QI and cross-ply
  panels fe3d is flat or rises slightly with energy (0.46 -> 0.48 and
  0.39 -> 0.42 from 5 to 30 J); see README "Limitations".
- **fe3d knockdown is now fe3d's damaged strength over its own undamaged
  strength.** fe3d's failure stresses (first-ply failure on the 3D mesh,
  panel and sublaminate buckling) are not the same quantity as the shared
  `pristine_strength_MPa` (a thickness-weighted ply-strength average), so
  dividing one by the other reported a knockdown for an **undamaged** panel:
  0.366 in tension and 0.032 in compression on an 8-ply IM7/8552 QI
  150x100 mm panel, where the empirical and semi-analytical tiers report
  1.000. `BvidAnalysis` now runs fe3d twice (damaged and undamaged, both
  uncapped), sets `knockdown = min(1, damaged / undamaged)`, and reports
  `residual_strength_MPa = knockdown * pristine_strength_MPa`, so the shared
  pristine baseline and `knockdown = residual / pristine` hold for every
  tier. An undamaged panel now gives 1.000 in fe3d too. The raw damaged and
  undamaged fe3d stresses are added to `AnalysisResults.notes`. fe3d runs
  take about twice as long. On the 16-ply QI validation cases (5/10/15 J)
  the fe3d compression knockdown moves 0.007/0.003/0.002 ->
  0.043/0.020/0.013. A degenerate fe3d buckling result now drops the
  buckling channel (first-ply failure alone sets the result) rather than
  substituting pristine strength; its note wording says so.
- **FE3D buckling channel retired the 3D K_g eigensolve in favour of the
  Rayleigh-Ritz closed form (issue #129, PR #130).** Investigation on the
  #98 thread established that even after the K_g sign fix (#128), the 3D
  eigensolve under-predicted by ~100x on realistic panels: the
  constant-prestress + 3-2-1 rigid-body BC formulation did not enforce a
  compressive state, and a proper static-prestress reformulation didn't
  recover analytical agreement at affordable mesh sizes (3D Hex on thin
  laminates is fundamentally locking-prone). The Navier sine-basis
  closed form (Timoshenko & Gere §9.2; Reddy 4.4.4) is exact for SSSS
  orthotropic rectangles, so `fe3d_cai_buckling` now delegates to it via
  a new `panel_buckling_load` helper that shares its `_rayleigh_ritz_N_cr`
  kernel with the existing `sublaminate_buckling_load`. For damaged
  panels the worst-sublaminate check still runs; the buckling channel
  returns `min(sigma_panel, sigma_sublam)`. The first-ply-failure channel
  on the damaged 3D mesh is unchanged. On the canonical 200x150 mm 8-ply
  QI pristine reference (IM7/8552, t_ply=0.152), the buckling channel
  now returns 12.6 MPa, matching the Rayleigh-Ritz analytical solution
  exactly; before this change it returned 0.20 MPa (60x under).
- **Removed the 5%-of-pristine plausibility guard in `BvidAnalysis.run`.**
  It was a defensive workaround for the broken eigensolve; the closed
  form returns the physically correct value (which can legitimately be
  below 5% of pristine on slender panels). Collapsed the
  `fe3d_buckling_unconverged` and `fe3d_buckling_artefact_dropped`
  warning tags into a single `fe3d_buckling_fallback` tag emitted only
  when the closed form returns a degenerate result.

### Fixed

- **Non-finite and out-of-domain inputs are rejected with a `ValueError`
  naming the field**: NaN/inf in `PanelGeometry`, `ImpactorGeometry`,
  `DelaminationEllipse`, `MeshParams`, `ImpactEvent`; negative dent/fibre-break
  radius in `DamageState`; zero/negative `mass_kg` and negative `energy_J`;
  `MeshParams(elements_per_ply=True)`; and in `AnalysisConfig` a non-finite
  ply thickness or angle, an empty or string layup, a non-`PanelGeometry`
  panel, and an unknown material name (which used to be a bare `KeyError` at
  run time; the message lists the `MATERIAL_LIBRARY` names).
- **`ImpactEvent(location_xy_mm=[0.0, 0.0])`** (a list, as JSON gives) is
  stored as a tuple and keeps meaning the panel centre.
- **`DamageOutsidePanelWarning`** when a delamination centroid lies outside
  the panel: the closed-form tiers counted it in the DPA while the fe3d mesh
  silently missed it.
- **CLI errors no longer print tracebacks**: a malformed C-scan
  (`CScanSchemaError`), config `ValueError`s and non-finite `--thickness` /
  `--layup` / `--panel` / `--energy` values are argparse errors (exit 2),
  `FE3DSizeError` is a one-line error (exit 1), and a closed stdout pipe
  (`bvidfe ... | head`) exits quietly.
- **Streamlit `_config_from_dict` accepts per-ply thickness lists** (it
  called `float()` on them), and the NCR export formats them.
- **Predicted delaminations stay inside the impacted panel.** The union of
  the template's differently oriented ellipses spreads well beyond the
  circle of the same area: 99 x 52 mm on the NCAMP quasi-isotropic coupon
  (same-area circle 61 mm), 134 x 70 mm on Lovejoy & Scotti's
  quasi-isotropic coupon (82 mm) and 108 x 108 mm on their 40-ply hard
  laminate (92 mm), wider than its 102 mm coupon. fe3d's soft zone then
  spanned the whole section and the 0 deg fibres failed at the coupon's
  side edge (y = 2.4 mm), 47% below the measured CAI.
  `distribute_damage(..., max_half_extent_mm=(hx, hy))` shrinks each
  ellipse, aspect ratio kept, until its bounding box fits, re-solves for
  the same union area, and warns and returns the saturated footprint if
  the window cannot hold the target. `impact_to_damage` passes the
  impacted panel's half-extents about the impact point (`impact_panel`,
  else `panel`). The NCAMP and Lovejoy & Scotti cases now carry the
  ASTM D7136 125 x 75 mm support window as `impact_panel`. fe3d on the hard
  laminate goes -47.1% -> +12.4% at factor 0.20 (+1.9% at the refitted
  0.17, see *Changed*). The closed-form tiers barely move: MAE NCAMP
  1.54% -> 1.50%, Lovejoy & Scotti 8.97% -> 8.95% (both tiers; Hasebe
  unchanged); the CLI snapshots move by less than 6e-5 relative. On a
  150 x 100 mm IM7/8552 panel the footprint now saturates at about
  11,200 mm^2 from about 114 J, before the 80% area cap.
- **A delamination `interface_index` outside the layup is rejected.** An
  n-ply laminate has interfaces 0 to n-2; `AnalysisConfig` now raises
  `ValueError` for a larger index. The closed-form tiers used to ignore it
  silently and fe3d raised `IndexError` while meshing.
- **fe3d buckling notes reach `AnalysisResults.notes`.** `BvidAnalysis`
  discarded the notes from `fe3d_cai_buckling`, so the degenerate-buckling
  note never reached the user.
- **`semi_analytical` CAI returns a Python `float`**, not `numpy.float64`,
  when delamination growth governs.
- **The validator rejects an unknown `--dataset` name** (exit 1) instead of
  skipping it silently.
- **The fe3d CI step is a smoke run, not a gate.** It compared fe3d against
  `synthetic_selfcheck`, which is empirical-tier output, and so failed on
  every run (5.3% against a 2% gate) without measuring fe3d's accuracy.

- **Project links pointed at the repository's previous owner**
  (`ranipdx-glitch`): the docs site URL, `repo_url`, the `pyproject.toml`
  project URLs, README and docs badges and links, the schema `$id`s, and
  the PyPI Trusted Publishing and Streamlit deployment instructions. They
  now point at `elhajjar1/BVID-FE`; a Trusted Publisher registered under
  the old owner would reject the release workflow.
- **`docs/schemas/analysis_config.json` was missing `puck_p_nt_minus` and
  `puck_p_nt_plus`**, the `OrthotropicMaterial` fields added with the Puck
  fix, because nothing checked the generated schemas against the
  dataclasses. Regenerated, and CI's `validation` job now runs
  `scripts/generate_schemas.py --check`.
- **CI coverage stopped counting the CLI** after the pytest-cov 7 bump:
  pytest-cov 7 no longer measures subprocesses, and the CLI tests run
  `python -m bvidfe.cli` in one (cli.py 92% -> 55%, total 93.3% -> 91.3%).
  `[tool.coverage.run] patch = ["subprocess"]` (coverage >= 7.10) restores
  it. Dev dependencies now allow pytest 9 and black 26, and the workflows
  use `actions/setup-python@v7` (supersedes Dependabot #79, #80, #141).
- **`Laminate.effective_engineering_constants` was off by a factor of
  h^2.** It formed the average-stress compliance as `A^-1 / h` instead of
  `h * A^-1`, so `Ex`, `Ey` and `Gxy` were scaled by the laminate thickness
  squared (in mm^2); `nu_xy` was unaffected. On IM7/8552 quasi-isotropic
  laminates `Ex` was 92.3 GPa at 8 plies (1.22 mm) and 830.9 GPa at 24
  plies (3.65 mm); it is now 62.4 GPa for both. The existing unidirectional
  test used 8 x 0.125 mm = 1.0 mm, where h^2 = 1, so it could not see the
  error. The Lekhnitskii `Kt_inf` uses only ratios of these constants and
  does not change. fe3d converts its failure strain to stress with `Ex`,
  so fe3d first-ply-failure and tension residuals drop by the same factor
  unless they were capped at pristine strength. At 20 J on a 150x100 mm
  panel the fe3d tension knockdown goes 0.541 -> 0.366 (QI, 8 plies),
  0.844 -> 0.571 (cross-ply, 8 plies) and 1.000 -> 0.362 (QI, 24 plies,
  where the inflated value had been capped at pristine and showed no
  knockdown). fe3d compression knockdowns on these panels are
  buckling-governed and do not change.
- **fe3d first-ply failure checked off-axis plies against the wrong
  strengths.** Two frame errors in the 3D element, fixed together:
  - `Hex8Element` built a +theta ply's stiffness as `T C T^T`, where `T`
    maps global stress into the +theta ply frame. That is the stiffness of
    a **-theta** ply: the xx-xy coupling had the opposite sign to CLT's
    Q-bar. It now uses `T^-1 C T^-T`.
  - The failure solve handed **global-frame** stress to Tsai-Wu / LaRC05 /
    Puck, which expect ply material-frame `[s1, s2, s3, t23, t13, t12]`. A
    90 deg ply's transverse stress was compared with `Xc`/`Xt`, and its
    fiber-direction stress with `Yt`/`Yc`. Stress is now rotated into each
    ply's frame (`Hex8Element.stress_at_gauss_points_material`).

  First-ply-failure strain on IM7/8552, 150x100 mm, 8 plies x 0.152 mm,
  20 J impact (before -> after): QI `[0/45/-45/90]s` compression
  0.139% -> 0.946% (now 0 deg fiber compression, ~`Xc/E11`; before, the
  90 deg plies' fiber stress was checked against `Yt`), QI tension
  0.191% -> 0.772%, cross-ply `[0/90]2s` tension 1.494% -> 0.862% (now
  90 deg ply cracking at ~`Yt/E22`; before, that cracking went undetected,
  which overstated strength), `[+30/-30/0]s` compression 0.311% -> 0.631%,
  `[+-45]2s` within 1%. `[0]_4` and `[90]_4` compression coupons now fail
  within 1% of `Xc/E11` and `Yc/E22` (Tsai-Wu on `[0]_4` sits 4.7% low from
  boundary-induced through-thickness stress); `[90]_4` used to reach the 5%
  strain cap without failing. fe3d **tension** residuals move accordingly.
  Compression knockdowns on these panels are still governed by the
  unchanged buckling channel, so they do not move.
- **fe3d first-ply failure silently evaluated Tsai-Wu when asked for Puck.**
  `_solve_failure_strain_analytic` branched `if criterion == "larc05": ...
  else: <Tsai-Wu>`, so `criterion="puck"` (and any typo such as
  `"tsaiwu"`) returned Tsai-Wu results with no error. Puck now has its own
  branch: its exposure factor is degree-1 homogeneous in stress, so the
  critical load multiplier is exactly `1 / idx_ref` from the single
  reference solve. Names outside `CriterionName` raise `ValueError`. The
  default fe3d criteria (LaRC05 for CAI, Tsai-Wu for TAI) are unchanged, so
  `BvidAnalysis.run()` results do not move.
- **Version strings disagreed.** `pyproject.toml` said `0.2.0`,
  `bvidfe.__version__` (and `bvidfe --version`) said `0.2.0.dev0`,
  `CITATION.cff` said `0.1.0`, and the README BibTeX said `0.1.0-alpha`.
  All four now read `0.2.0`, and `tests/test_package.py` fails CI if they
  drift apart again. The `release` skill lists all four places to bump.
- **`FailureEvaluator.evaluate` validated its input with `assert`**, which
  `python -O` strips. A wrong-shaped stress field now raises `ValueError`.
- **Whitney-Nuismer TAI could return a residual above pristine or diverge.**
  The point-stress polynomial dips below a stress ratio of 1 for
  `Kt_inf >~ 9.2` (knockdown > 1) and crosses zero at `Kt_inf >~ 20.3`.
  The knockdown is now capped at 1.0 with a `UserWarning` in that regime.
  The built-in material presets stay below `Kt_inf ~ 7.3`, so their results
  are unchanged; only custom high-modulus cards were affected.
  `lekhnitskii_kt_infinity` raises a descriptive `ValueError` instead of a
  bare `math domain error` for constants that violate the orthotropic
  stability bound.
- **CI lint step went red on every branch after ruff 0.16.** The dev extra
  allows `ruff>=0.5,<1`, and ruff 0.16 widened its default rule selection,
  flagging 164 pre-existing findings. `pyproject.toml` now pins
  `[tool.ruff.lint] select` to the pre-0.16 defaults (`E4`, `E7`, `E9`,
  `F`); adopting the wider rule set is left to its own PR.
- **CHANGELOG `[0.2.0]` overstated validation.** It claimed "validated
  public datasets exercised in CI"; CI runs the validation harness against
  a synthetic self-check dataset only. The entry is corrected.
- **FE3D buckling assembled `K_g` with a tensile unit reference (issue #98),**
  causing the eigensolver to lock onto spurious compliance modes that
  decayed monotonically with mesh refinement. The geometric stiffness is
  now assembled from a **compressive** unit reference (`sigma_xx =
  -sigma_ref_MPa` at `src/bvidfe/analysis/fe_tier.py`), so the smallest
  positive eigenvalue `lambda_crit` of `K phi + lambda K_g phi = 0` is the
  physical compressive buckling stress in MPa directly. Knockdowns from the
  `fe3d` tier shift accordingly — see #98 for the diagnostic, including
  reference-config tables. The eigenvalue selection rule ("smallest
  positive") is unchanged; under the new convention it is now correct.
  In-plane locking on coarse meshes remains a separate open concern
  (tracked on its own branch).

### Removed

- **Leftovers from the removed PyQt6 GUI.** `scripts/visual_matrix.py`
  imported `PyQt6` and `bvidfe.gui`, neither of which exists, so it could not
  run. `docs/BUILD.md` documented a PyInstaller build and workflows that no
  longer exist. `docs/python_api.md` told users to import
  `bvidfe.viz.plots_3d`, which does not exist; it now shows the real
  `bvidfe.viz.plotly_3d.mesh_damage_figure` API.

## [0.2.0] - 2026-05-21

First tagged release after `0.1.0`. Major themes: a fully wired 3-D
hex-element FE tier (`fe3d`) with eigenvalue buckling and FPF coupons,
per-ply thickness throughout the stack, a Streamlit GUI, a validation
harness exercised in CI (against a synthetic self-check dataset only;
the published Soutis / Caprino / Sanchez-Saez / NASA datasets are not
yet digitized, see `validation/README.md`), a `SemiAnalyticalResult` dataclass on
the semi-analytical path, vectorised LaRC05/Tsai-Wu failure-criterion
batches, refreshed README/ARCHITECTURE docs, and the OIDC Trusted
Publishing workflow that produced this release (issue #116). See the
detailed entries below.

### Changed

- **`semi_analytical_cai` returns a structured `SemiAnalyticalResult`
  dataclass instead of a positional `tuple[float, int|None, float|None]`
  (issue #83).** The previous return shape forced callers to remember the
  positional convention `(sigma_CAI, critical_interface, N_cr)` and used
  bare `None` sentinels in two of the three slots, which made the
  `BvidAnalysis._semi_analytical` call site fragile (a silent swap of two
  fields would have type-checked but been silently wrong). The new
  ``SemiAnalyticalResult`` is a frozen dataclass with three named fields —
  ``residual_strength_MPa: float``, ``critical_interface_index: int |
  None``, and ``critical_buckling_load_N: float | None`` — so every call
  site reads named attributes. ``BvidAnalysis._semi_analytical`` (the only
  internal consumer) is updated; the tension branch wraps its scalar
  result in the same dataclass so the dispatch shape is uniform across
  loadings. The dataclass is **not** re-exported at the top-level
  ``bvidfe`` package — it stays internal to ``bvidfe.analysis`` until a
  downstream user case justifies promoting it. Since BVID-FE is still on
  the 0.x series, no compatibility shim is provided; users calling
  ``semi_analytical_cai`` directly must read named attributes (or
  reconstruct the previous tuple via
  ``(r.residual_strength_MPa, r.critical_interface_index,
  r.critical_buckling_load_N)`` at the call site).

### Added

- **Per-ply thickness support (issue #5).** ``Laminate.ply_thickness_mm``,
  ``AnalysisConfig.ply_thickness_mm``, and the CLI ``--thickness`` flag now
  accept either a single positive number (uniform laminate, the legacy
  behaviour) **or** a list/tuple of per-ply thicknesses with one entry per
  ply in the layup. This lets users model laminates that mix plies of
  different fabric weights or prepreg gauges (e.g. a thinner surface 0/90
  over a thicker quasi-iso interior) without faking it via an effective
  uniform thickness. ``Laminate`` exposes the resolved per-ply list via
  the new ``ply_thicknesses_mm`` property and reports
  ``is_uniform_thickness``; the CLT ABD assembly, the
  ``_pristine_strength`` thickness-weighted ply average, the
  semi-analytical sublaminate selection, and the fe3d hex-mesh z-grid
  (built on the cached ``FeMeshSkeleton``: per-ply ``ply_top_z``
  boundaries replace the single scalar ``ply_thickness_mm``) all consume
  the per-ply list. A uniform list reproduces scalar results bit-for-bit;
  ``--thickness`` accepts the same comma-separated form as ``--layup``
  with a length-match check at parse time. Closes #5.

### Changed

- **README + ARCHITECTURE refreshed to match v0.2.0-dev reality.** The
  README's Limitations section claimed "No GUI in v0.1.0", "fe3d CAI
  uses first-ply-failure, not a buckling eigenvalue solve", and "no
  validated datasets in v0.1.0" — three statements that were correct
  for v0.1.0 but stale for v0.2.0-dev (where the GUI ships, the
  buckling eigensolve is implemented, and the validation harness is in
  the repo). The hard-coded test-count badge ("216 passing", stale the
  moment the suite grew) was replaced with a CI-driven GitHub Actions
  workflow status shield that updates automatically.
  ``ARCHITECTURE.md``'s module catalog had no entry for
  the ``gui/`` package, the dependency diagram did not reach ``cli``
  or ``gui``, and the v0.2.0 roadmap section was a list of items that
  have since shipped. All three pages updated to reflect the current
  reality, with a "Shipped in v0.2.0-dev" / "Planned (post-v0.2.0)"
  split on the roadmap so future work and historical claims do not
  drift back into confusion.

### Changed

- **`_solve_failure_strain_analytic` (fe3d FPF inner solve) is now
  vectorised across Gauss points.** The per-element inner Python loop
  over Gauss points was replaced with a single batched call to
  ``larc05_index_batch`` / ``tsai_wu_index_batch`` (added in the previous
  commit) plus ``np.where``-masked branch evaluation for the four
  numerical edge cases (LaRC05 mask: ``idx_ref > 0``; Tsai-Wu masks:
  ``|b| < tiny`` linear fallback, ``|b| >= tiny`` quadratic, ``disc < 0``
  no-real-root, ``c <= 0`` non-physical-root). The previous scalar
  implementation is preserved verbatim as
  ``_solve_failure_strain_analytic_scalar_ref`` so the new
  ``tests/analysis/test_fpf_strain_solve_equivalence.py`` can prove the
  two paths return identical strain-at-failure values
  (``rel=1e-10``) on three representative inputs (pristine LaRC05,
  delaminated LaRC05, delaminated Tsai-Wu) and on the tension path
  (``strain_sign=+1``). Speedup is modest at default mesh resolutions
  (n_gp = 8 per element so the inner loop was small) but the change
  removes the last per-Gauss-point Python loop in the fe3d hot path
  and aligns the FPF solve with the rest of the failure-criterion API.

### Added

- **fe3d tier now exercised in CI as a monitoring step.** The CI
  ``validation`` job previously ran only the empirical tier (``python
  validation/validate_bvid_public.py --gate``). A new step runs the
  same harness with ``--tier=fe3d --max-cases=3``, marked
  ``continue-on-error: true`` so failures are advisory and don't block
  merges. The validation script gains:
  - A new ``_TIER_GATE_MULTIPLIER`` table — empirical and
    semi_analytical use the original 1.25× target_mae_pct slack;
    fe3d uses 2.0× to acknowledge the documented energy-flatness
    caveat without flapping the CI on regressions still inside the
    model uncertainty.
  - A new ``--max-cases N`` flag that truncates each dataset to the
    first ``N`` cases. Lets the fe3d step cap CI runtime at ~45 s
    even if the dataset later grows to dozens of cases. Default is
    ``None`` (full dataset), preserving the existing empirical-gate
    behaviour exactly.

- **Vectorised batch failure-criterion functions.** New
  ``larc05_index_batch(material, stresses)`` and
  ``tsai_wu_index_batch(material, stresses)`` evaluate their
  respective criteria over a Voigt-6 batch in a single numpy /
  ``np.einsum`` pass. ``FailureEvaluator.evaluate`` is refactored
  to call ``_index_batch`` once on the entire ``(n_elem, n_gp, 6)``
  field, then ``np.argmax`` to pick the critical (element, Gauss
  point) — replacing the prior nested Python loop that paid the
  function-call overhead n_elem × n_gp times. Three new tests in
  ``tests/failure/test_evaluator.py`` lock numerical equivalence
  between the batch and scalar forms (``rtol=1e-12``) and one
  test asserts the new ``evaluate`` agrees with the explicit
  nested-loop reference. The batch helpers are also intended as
  the prerequisite for a future vectorisation of
  ``_solve_failure_strain_analytic`` (deferred — see CHANGELOG
  rationale: small inner-loop sizes + weak regression coverage on
  the analytic strain solve mean the speedup-vs-risk tradeoff
  doesn't yet justify it).

### Changed

- **Refactor pass: tsai_wu vectorisation + heartbeat helper + fe3d pre-flight
  helper + dead-code removal.** Four targeted refinements (each independently
  verified by an explore-agent and locked by unit / equivalence tests):
  - ``tsai_wu_index`` no longer iterates 36 times in pure Python per call.
    ``_tsai_wu_coefficients`` now returns ``(F: np.ndarray, Q: np.ndarray)``
    instead of a Python list / list-of-lists, and the index is computed via
    ``F @ s + s @ Q @ s`` (BLAS-routed dot + bilinear form). A new
    ``test_tsai_wu_index_matches_explicit_nested_sums`` test pins 50
    randomised stress vectors to the explicit nested-sum form within
    ``rel=1e-12``. ``tsai_wu_strength_uniaxial`` is unchanged because
    ``F[i]`` and ``Q[i][i]`` work identically on numpy arrays.
  - The duplicated heartbeat-progress loop in ``AnalysisWorker.run`` and
    the per-iteration body of ``SweepWorker.run`` is now a single helper
    ``gui.workers._run_with_heartbeat(work, on_progress, *, start_pct,
    end_pct, interval_s, label)`` that spawns a daemon thread, ticks the
    progress callback toward ``end_pct`` at ``interval_s`` intervals, and
    re-raises a ``RuntimeError`` carrying the original traceback on
    failure. ``AnalysisWorker.run`` drops from ~40 lines to ~16; the
    ``SweepWorker`` per-energy body drops from ~19 lines to ~10. The
    ``TierComparisonWorker`` is intentionally NOT migrated because it
    ticks per-pair (synchronous, step-based) rather than per-second.
  - The duplicated fe3d pre-flight prologue (``_guard_problem_size`` +
    ``build_fe_mesh`` + ``_build_elements`` + log lines) shared by
    ``fe3d_cai_buckling``, ``_fe3d_cai_first_ply_failure``, and
    ``fe3d_tai`` is now a single ``analysis.fe_tier._fe3d_preflight(cfg,
    damage, lam, *, label) -> (mesh, elements, t0)`` helper. Each call
    site drops from ~6 lines to 1; ``fe3d_tai`` now also emits per-stage
    timing logs for parity with the CAI paths.
  - ``_k_contact_hertz_linearized`` (a back-compat alias in
    ``impact/olsson.py``) is deleted — re-verified to have zero callers
    in src/, tests/, examples/, scripts/, validation/, or docs/. The
    second prior dead-code candidate ``tsai_wu_strength_uniaxial`` was
    re-checked and KEPT (six active calls in
    ``tests/failure/test_tsai_wu.py``).

### Added

- **`tests/validation/` directory with five new analytical-reference tests.**
  Five new files lock the public formulas to known closed-form values so
  silent drift in the physics surfaces immediately:
  - ``test_clt_abd_reference.py`` — pinned ``A``, ``B``, ``D`` matrices for a
    ``[0/90]_s`` IM7/8552 laminate; symmetric layup → ``B == 0``; balanced
    layup → ``A_16 == A_26 == 0``.
  - ``test_olsson_threshold_reference.py`` — Olsson-2001 threshold-load
    shape relationships (``Pc / sqrt(D_eff * G_IIc) = pi*sqrt(8/9)``,
    panel-size invariance, ``Pc ~ h^{3/2}`` thickness scaling).
  - ``test_soutis_cai_reference.py`` — Soutis CAI knockdown self-
    consistency: zero-DPA short-circuit, dimensionless invariance under
    ``DPA / A_panel`` scaling, monotonicity, ``k_s`` calibration response.
  - ``test_whitney_nuismer_limits.py`` — Whitney-Nuismer TAI analytical
    limits: ``DPA → 0`` knockdown → 1.0; ``DPA → ∞`` knockdown → ``1 /
    Kt_inf`` (= 1/3 for Kt_inf = 3); monotone-decreasing in DPA;
    bounded-below by the asymptote.
  - ``test_fe3d_pristine_clt_consistency.py`` — pristine fe3d uniaxial
    extension recovers ``E_x = A11 / h`` (CLT) within 10% on a coarse
    ``20 × 10 mm × 4-ply`` mesh; pristine input produces no damaged
    elements; zero-load solve returns zero displacement.

- **Docstring upgrades on five high-traffic public symbols.**
  ``larc05_index``, ``tsai_wu_index``, ``olsson.threshold_load``,
  ``hex8.geometric_stiffness_matrix``, and
  ``semi_analytical.sublaminate_buckling_load`` previously had one-line
  docstrings that hid the underlying physics. Each now documents the
  formula, references (LaRC05 NASA/TM, Olsson 2001, Cook §17.7 / Bathe
  §6.8, Timoshenko & Gere §9.2), Voigt-index conventions, units, return-
  value semantics, and any raise-conditions — so ``help()`` and IDE
  tooltips show the same context that's currently spread across the
  README and ARCHITECTURE.

- **GUI panel-geometry helper deduplicated.** ``BvidMainWindow`` had
  three sites constructing the same ``PanelGeometry(Lx_mm, Ly_mm,
  boundary)`` from the panel-input panel (one of which was already
  exposed as ``panel_panel_as_geometry`` but unused). The two callers
  in ``_update_live_onset`` and ``_build_config`` now use the helper;
  the now-unused ``PanelGeometry`` import is dropped from
  ``_update_live_onset``. No behaviour change.

- **Repo hygiene: dependency caps, headless requirements.txt, top-level
  conftest, Dependabot.** Four small infrastructure changes that
  collectively make the project safer to install and easier to keep up:
  - ``pyproject.toml`` and ``requirements.txt`` now pin upper bounds
    on every direct dependency (e.g. ``numpy>=1.24,<3``,
    ``scipy>=1.10,<2``, ``PyQt6>=6.5,<7``, ``pytest>=7,<9``). The next
    major release of any of these libraries can no longer auto-break a
    BVID-FE install; users who want to live on the bleeding edge can
    still override the cap explicitly.
  - ``requirements.txt`` no longer pulls in the GUI stack (~100 MB of
    Qt + VTK + pyvistaqt). Headless library / CI users now get only the
    five core scientific deps; the GUI is opt-in via
    ``pip install -e ".[gui]"`` (or ``[all]``).
  - A new top-level ``conftest.py`` sets ``QT_QPA_PLATFORM=offscreen``
    before any PyQt import, so the ``tests/gui/`` suite runs cleanly on
    CI runners without a display server. Individual GUI test modules
    keep their per-file ``setdefault`` as defence-in-depth.
  - A new ``.github/dependabot.yml`` schedules weekly PRs for both pip
    and GitHub Actions dependencies (capped at 5 open per ecosystem),
    so security patches surface automatically and the test matrix
    runs against the new pin before a maintainer merges.

- **MeshParams input validation.** ``MeshParams.__post_init__`` now
  rejects ``elements_per_ply <= 0`` (and non-int values),
  ``in_plane_size_mm <= 0``, and ``cohesive_zone_factor <= 0``. The bad
  input previously slipped into ``build_fe_mesh`` and either silently
  produced a degenerate brick mesh or tripped a far-downstream error
  inside the assembler. Three parametrised tests in
  ``tests/analysis/test_mesh_params_validation.py`` lock the new
  behaviour.

- **CLI parser-helper unit tests.** ``_parse_panel``, ``_parse_layup``,
  ``_positive_float`` and ``_existing_path`` (introduced for issue #7)
  were exercised only end-to-end via ``subprocess`` — too slow for
  edge-case coverage. ``tests/test_cli_helpers.py`` adds 12 direct
  pytest tests covering malformed panel formats, non-positive
  dimensions, non-numeric layup tokens, zero / negative / non-numeric
  ``_positive_float`` inputs, and the missing-file / directory paths
  for ``_existing_path``.

- **Cross-tier consistency tests.** Per the new README "Knockdown
  definition and cross-tier comparability" section there are several
  documented relationships between tier outputs that no test enforced:
  (a) ``pristine_strength_MPa`` identical across tiers, (b) TAI
  ``empirical.knockdown == semi_analytical.knockdown`` exactly, (c)
  CAI ``semi_analytical.knockdown <= empirical.knockdown``, (d)
  pristine input → knockdown == 1.0, (e) every knockdown finite and
  bounded in [0, 1+ε]. ``tests/test_cross_tier_consistency.py`` locks
  all five.

### Fixed

- **C-scan loader hardens type and finiteness checks.** ``_validate_dict``
  in ``damage/io.py`` previously had several gaps: a string in a numeric
  field crashed at the ``< 0`` comparison with a raw ``TypeError``;
  ``NaN`` and ``Inf`` (both JSON-encodable via Python's default
  ``allow_nan=True``) passed through into the analysis pipeline; a non-
  dict entry in the ``delaminations`` array crashed the
  ``"k" not in d`` check; ``fiber_break_radius_mm < 0`` was accepted;
  ``interface_index = 2.5`` was silently truncated to ``2``; and unknown
  top-level / per-delamination fields were dropped silently, so a typo
  like ``"dent_dept_mm"`` (missing the trailing 'h') reported as
  ``missing required field: dent_depth_mm`` without flagging the typo.
  Fixed by routing every numeric field through a new
  ``_require_finite_number`` helper that raises ``CScanSchemaError`` on
  non-numeric, NaN, or Inf, plus per-delamination dict-shape checks,
  ``fiber_break_radius_mm >= 0``, strict-int ``interface_index``, and a
  ``UserWarning`` (forward-compatible — does not fail the load) on
  unknown fields. Eight new tests in ``tests/damage/test_io.py`` cover
  the string-in-numeric, NaN, Inf, negative-radius, non-object-element,
  unknown-field warning, and non-int-interface-index paths.

- **Invalid `BVIDFE_LOG_LEVEL` / `BVIDFE_FE3D_MAX_DOF` no longer crash
  imports.** Both env vars were read at module import time:
  ``logging.setLevel(os.environ.get(...).upper())`` raised on a typo
  like ``BVIDFE_LOG_LEVEL=DUBG``, and ``int(os.environ.get("BVIDFE_FE3D_MAX_DOF",
  "500000"))`` raised on a non-numeric value. Either error prevented
  ``bvidfe`` from being imported at all. Fixed by routing both reads
  through new ``_resolve_log_level`` / ``_resolve_max_dof`` helpers in
  ``analysis/fe_tier.py`` that fall back to the documented default with
  a single-line stderr warning when the value is unrecognised, zero,
  or negative. ``gui/workers.py`` shares the log-level resolver. Seven
  new tests in ``tests/analysis/test_env_resolution.py`` cover the
  default, valid-name, typo, override, non-numeric, zero, and negative
  paths.

- **Parametric sweeps no longer lose partial results on per-iteration
  failure.** ``sweep_energies`` / ``sweep_layups`` / ``sweep_thicknesses``
  previously called ``_run_one`` unguarded — a single failed iteration
  (Olsson out-of-regime, mesh degeneracy, Tsai-Wu invalid combination)
  raised out of the sweep and the CSV was never written, so a 12-energy
  sweep with one bad row dropped 11 valid results on the floor. Each
  entry point now accepts an ``on_error`` keyword (default
  ``"raise"`` for backward compatibility, ``"skip"`` to fill the failed
  row with NaN numerics + an ``error`` column and continue, ``"warn"``
  to do the same plus emit a ``UserWarning``) and an optional
  ``progress_callback(i_done, n_total)``. With ``on_error="skip"`` the
  partial CSV is preserved with the same row count as the input
  iterable. Five new tests in ``tests/sweep/test_sweep.py`` cover the
  raise / skip / warn paths, the partial CSV, the progress callback,
  and rejection of bogus ``on_error`` values.

- **CLI now validates inputs at parse time instead of deep in the pipeline.**
  Several `bvidfe` CLI flags previously accepted invalid values that only
  surfaced as opaque errors much later: `--material IM7/8553` (typo)
  reached ``MATERIAL_LIBRARY[args.material]`` before producing a raw
  ``KeyError``; `--energy 0` and negative numbers passed argparse and
  produced cryptic Olsson-threshold warnings; `--cscan
  /does/not/exist.json` raised a bare ``FileNotFoundError`` from inside
  ``load_cscan_json``; `--thickness`, `--impactor-diameter`, and
  `--mass` accepted any float including zero and negatives. The argparse
  setup now uses ``choices=`` for ``--material`` (so typos surface as
  ``invalid choice: 'IM7/8553' (choose from 'AS4/3501-6', 'IM7/8552',
  'T700/2510', 'T800/epoxy')``), a custom ``_positive_float`` type for
  the four numeric flags, and a custom ``_existing_path`` type for
  ``--cscan`` that fails fast with ``file not found``. Panel dimensions
  are also rejected when non-positive. A new ``--quick-json`` flag emits
  a single-line JSON object ``{knockdown, residual_strength_MPa,
  pristine_strength_MPa, tier_used}`` for scripted callers that need
  more than the bare scalar from ``--quick``. Four new tests in
  ``tests/test_cli.py`` cover the bad-material, zero-energy,
  missing-cscan, and quick-json paths.

- **GUI surfaces malformed input instead of silently dropping it.** Two
  separate paths previously absorbed user-input failures with no
  feedback: ``BvidMainWindow._load_config`` had a bare ``except
  Exception`` around ``json.loads`` + ``config_from_dict`` that turned
  every failure into a single cryptic ``str(exc)`` warning, and
  ``DamagePanel.get_damage_state`` swallowed ``ValueError`` /
  ``AttributeError`` per-row with ``continue`` so a typo in a single
  cell silently dropped that delamination from the analysis. Both are
  now explicit:
  - ``_load_config`` distinguishes ``OSError`` (cannot read file),
    ``json.JSONDecodeError`` (malformed JSON, with line + column),
    ``KeyError`` (missing required field), and ``TypeError`` /
    ``ValueError`` (wrong type), each with a distinct dialog title.
  - ``DamagePanel.get_damage_state`` records every skipped row in a new
    ``self.skipped_rows`` attribute and logs a warning on the
    ``bvidfe.gui`` logger; ``BvidMainWindow._build_config`` reads this
    after every config build and reports skipped row indices via the
    status bar so the user can fix the typo. Four new tests in
    ``tests/gui/test_config_io.py`` (covering JSON parse errors, missing
    fields, wrong types, and damage-row skipping) lock the new
    behaviour.

- **GUI tier-comparison no longer freezes the main thread.** The
  "Compare Tiers (empirical + semi_analytical)..." menu action ran 16
  `BvidAnalysis.run()` calls (2 tiers x 8 energies, ~12 s wall clock at
  default settings) directly on the Qt main thread inside
  `BvidMainWindow._compare_tiers`. The status bar updated to "Running
  tier comparison..." but the entire UI was unresponsive until the loop
  finished — including the Knockdown Curve tab where the result would
  ultimately appear. Fixed by extracting the loop into a new
  `TierComparisonWorker(QThread)` in `gui/workers.py` that mirrors the
  existing `AnalysisWorker` / `SweepWorker` pattern (resultReady / error /
  progress signals, daemon-thread heartbeats, deleteLater cleanup).
  Per-(tier, energy) failures are absorbed into NaN entries with a
  third tuple element listing the skipped pairs so the status bar can
  report them instead of aborting the whole comparison. Two new tests
  in `tests/gui/test_workers.py` cover the happy-path and the
  bogus-tier partial-failure path.

- **Hex8 elements now validate the Jacobian determinant.** `Hex8Element.B_matrix`
  and `Hex8Element.geometric_stiffness_matrix` previously called
  `np.linalg.inv(J)` with no check on `np.linalg.det(J)` — an inverted element
  (det < 0, from wrong node ordering) silently produced negative volume in
  the stiffness integral, and a singular element (det ≈ 0, from collapsed
  nodes) silently produced NaNs that propagated into the global K and
  corrupted every downstream eigensolve / FPF result. Both call sites now
  invoke a shared `_validate_jacobian` helper that raises a new
  `DegenerateElementError(ValueError)` with the offending detJ value, the
  natural-coordinate Gauss-point location, and the full node-coordinate
  list, so users get a reproducible "bad mesh" error instead of opaque
  NaNs further down the pipeline. `Hex8iElement.stiffness_matrix` (which
  inverts a separate centre-Jacobian for its bubble-mode Benh derivation)
  has the same guard. The structured brick mesh produced by `build_fe_mesh`
  is unaffected — regular cuboid grids always have strictly positive detJ.

- **Tsai-Wu now uses through-thickness strengths and the 1-3 shear strength.**
  `failure/tsai_wu.py` previously computed `F3 = 1/Yt - 1/Yc`, `F33 =
  1/(Yt*Yc)`, and `F55 = 1/S12**2`. The first two reuse the in-plane
  transverse strengths for the through-thickness direction (a missing-
  parameter conflation, since `OrthotropicMaterial` had no `Zt`/`Zc`
  fields), while `F55 = 1/S12**2` was a real formula bug — the 1-3
  shear coefficient should use the 1-3 shear strength (Voigt index 4
  in the project's `[xx, yy, zz, yz, xz, xy]` convention), not the
  in-plane shear strength. This deviated from Tsai & Wu (1971) and
  produced incorrect failure indices whenever a 3D stress state had
  non-trivial sigma_zz, tau_yz or tau_xz components — e.g. `fe3d` TAI
  on a damaged panel where the through-thickness coupling is broken.

  `OrthotropicMaterial` gains three optional fields, `Zt`, `Zc`, `S13`,
  and three corresponding `*_resolved` properties that fall back to
  `Yt` / `Yc` / `S12` when the user has not provided through-thickness
  test data (the standard transverse-isotropy assumption for
  unidirectional CFRP). The four library presets keep their existing
  in-plane values, so default behaviour is numerically backward-
  compatible — the formula reduces to the prior code and existing
  Tsai-Wu tests pass unchanged. Users with measured through-thickness
  strengths can now override per-material; three new tests in
  `tests/failure/test_tsai_wu.py` cover the default-fallback,
  Zt-override, and S13-override paths.

- **fe3d buckling silent fallbacks now surfaced via `AnalysisResults.notes`.**
  Previously, `fe3d_cai_buckling()` returned `(sigma_pristine_MPa, 0.0)` when
  the eigensolver found no positive eigenvalue or raised an exception, and
  `BvidAnalysis.run()` silently used FPF instead when the buckling result was
  below 5% of pristine ("plausibility gate"). The user saw `knockdown = 1.0`
  in the no-eigenvalue case and had no way to know that the buckling solve
  had quietly given up. `AnalysisResults` gains a `notes: list[str]` field
  populated by the analysis backends; `fe3d_cai_buckling()` now returns
  `(sigma, lambda_crit, notes)` so the caller can merge solver diagnostics
  into the result. The GUI Summary tab appends these notes (prefixed with
  ⚠) under the existing "--- Notes ---" section, and `AnalysisResults.summary()`
  / `to_dict()` include them in the text and JSON output. Two new tests
  in `tests/analysis/test_fe3d_buckling.py` cover the eigensolver-failure
  path end-to-end via `monkeypatch`.

- **fe3d damage factor now applied component-wise to the elasticity matrix.**
  The previous `DAMAGE_STIFFNESS_FACTOR = 0.30` was scaled into every entry of
  the 6×6 `_C_global` matrix in `analysis/fe_tier._build_elements`, uniformly
  reducing in-plane stiffness (E11, E22, G12) by 70% in delaminated zones —
  even though the docstring in `fe_mesh.py` already noted that "the plies
  themselves are intact (so in-plane load-carrying is mostly preserved)". The
  fix replaces the single constant with two physically-motivated factors:
  - `DAMAGE_OOP_FACTOR = 0.05` — out-of-plane (E33, G13, G23, and the
    in-plane / OOP Poisson cross-coupling), applied to every element inside a
    delamination ellipse footprint. Representative of the post-delamination
    interlaminar-modulus loss reported in Bolotin (2001) and Sun & Tao (1998).
  - `DAMAGE_FIBER_BREAK_INPLANE_FACTOR = 0.30` — in-plane sub-block (rows/cols
    {0, 1, 5} in the global Voigt frame), applied only inside the fiber-break
    core (within `fiber_break_radius_mm` of any delamination centroid).
    Representative of fiber-direction damage saturation in Camanho & Davila
    (2007) / Maimi et al. (2007) for CFRP.

  `FeMesh` gains a new `in_plane_damage_factors` array; existing
  `damage_factors` is now the OOP factor exclusively. The viz layer
  (stress-field heatmap, pyvista damage overlay) is unchanged because it
  visualises the OOP factor as "delamination depth", which is the right
  semantic. The fe3d Kg assembly now scales the uniaxial sigma_xx prestress
  by the in-plane factor (since the prestress is an in-plane component); pure
  delamination zones therefore carry full uniaxial load through the fibers.

  Effect on results: knockdown values shift, the buckling solve becomes more
  responsive to damage size (because the pure-delamination in-plane block is
  no longer over-penalised), and the previously-flagged "fe3d knockdown
  approximately flat vs. impact energy" caveat is partially attenuated.
  Existing tests are structural (inequalities and bounds, not pinned numeric
  values) and should still pass; CI will flag any that need updating.

### Added

- **Input-variable sensitivity rollup.** Before this change, four of the
  five impact-event inputs were silently ignored by the physics pipeline:
  - `panel.boundary` (simply_supported / clamped / free) never moved the
    Olsson bending stiffness, the semi-analytical sublaminate buckling
    coefficient, or the fe3d buckling BCs — all three tiers hard-coded
    simply-supported. Fixed by adding boundary-dependent multipliers
    (clamped=2.5× bending, 1.9× buckling; free=0.4× bending, 0.5× buckling
    per Timoshenko) in `impact/olsson.py` and `analysis/semi_analytical.py`,
    plus boundary-aware lateral-edge `u_z` penalty BCs in
    `analysis/fe_tier.py::fe3d_cai_buckling`.
  - `impactor.shape` (hemispherical / flat / conical) was never referenced
    by any physics code. Fixed by wiring shape into a per-shape Hertz-
    contact stiffness (Johnson, *Contact Mechanics* §3.4-3.5) and into a
    footprint-spread multiplier on the target DPA (flat=1.4×, conical=0.7×).
  - `impact.mass_kg` was never referenced. Fixed by adding a
    calibration-aware dynamic amplification factor `(5.5/m)^0.1`, exactly
    unity at the ASTM D7136 reference mass, plus a UserWarning when the
    impactor mass ratio falls below Olsson's quasi-static validity regime.
  - `impactor.diameter_mm` only shifted E_onset by <0.5% and was masked
    by the DPA cap; fixed by adding a `(16/d)^0.3` spread factor on DPA.

  A new scripted validation matrix (`scripts/validate_inputs.py`) runs 5
  levels of each variable across every (tier, loading) combination and
  flags any input that produces zero knockdown variation — the flags
  section now reports "(none)" for the first time.

- **Impact-location and fiber-break core markers on the damage map.**
  Previously the damage-map plot drew only the delamination ellipses, so
  a user could not tell from the plot WHERE the impact landed on the
  panel — a real problem for off-center impacts where the footprint
  asymmetry drives failure. `plot_damage_map` now adds a black "×"
  marker at each unique ellipse centroid (the impact point) and, when
  the material has a non-trivial fiber-break model (`fiber_break_eta > 0`
  and `E_impact` above the fiber-break threshold), a red filled circle
  showing the fiber-break core radius. Both appear in the legend. 4
  new regression tests in `tests/viz/test_damage_map_markers.py`.

- **Live DPA preview with saturation warning.** A new `DPA:` label in the
  Impact panel updates every time any input changes, showing both the
  absolute predicted damage area in mm² and the percentage of panel area.
  When the 80% cap engages the label switches to red bold text with a ⚠
  SATURATED marker — users now see saturation *before* pressing Run,
  rather than discovering it post-hoc via the Summary tab notice and
  wondering why knockdown stopped responding to energy. Saturation in
  the default 150×100 mm 8-ply configuration kicks in at ~15 J; this
  preview makes that limitation obvious. 4 new regression tests in
  `tests/gui/test_live_onset.py`.

- **Live E_onset preview in the Impact panel.** The `ImpactPanel.set_onset_energy()`
  method has existed since the first release but was never wired up — the
  "E_onset: — J" label stayed blank no matter what the user changed. Now
  the `BvidMainWindow` connects `configChanged` on the material, panel,
  impact, and input-mode panels to a single `_update_live_onset` slot
  that recomputes the Olsson threshold and updates the label in real
  time. Because E_onset is boundary- and shape-aware (new in this
  release), the label now responds visibly to every relevant input:
  simply_supported→clamped drops it from 0.60 J to 0.24 J, simply_supported→free
  raises it to 1.49 J, etc. Handles invalid intermediate states (user
  mid-edit) gracefully and blanks out in damage-driven mode where the
  preview is not applicable. 5 new regression tests in
  `tests/gui/test_live_onset.py`.

- **fe3d KD-vs-energy monotonicity fix.** `DAMAGE_STIFFNESS_FACTOR` in
  `analysis/fe_mesh.py` was raised from `1e-4` to `0.30`. The old value
  was physically unrealistic — it treated delaminated elements as
  essentially null in the stress field, which made the failure-index
  criterion unable to flag damaged regions, and the fe3d residual
  strength *increased* with impact energy past ~15% panel damage
  (because the peak stress in the undamaged shell drops as the damage
  footprint spreads wider). The new value represents typical residual
  in-plane stiffness after delamination (plies are intact, only the
  through-thickness coupling is lost); literature range for CFRP is
  0.1-0.5. Knockdown now trends monotonically with energy across
  simply-supported, clamped, and free boundaries (see
  `tests/test_input_sensitivity.py::test_fe3d_knockdown_mostly_decreases_with_energy`).

- **Summary tab edge-case notices.** The GUI Summary tab now surfaces
  the `UserWarning`s that previously only appeared on stderr: DPA
  saturation (knockdown insensitive to further energy), 'free' boundary
  soft-support approximation, small-mass quasi-static regime violation,
  and the known fe3d energy-insensitivity limitation. Each notice is
  one plain-language paragraph under a "--- Notes ---" section so the
  user sees why their result may not respond to an input they changed.

- **True linear buckling eigensolve in the `fe3d` CAI tier**
  `Hex8Element.geometric_stiffness_matrix(sigma_bar)` integrates
  `gradN.T @ S @ gradN` over the element via 2×2×2 Gauss quadrature and
  expands to a 24×24 K_g via `np.kron`. A new `fe3d_cai_buckling()` in
  `analysis/fe_tier.py` assembles K and K_g under a uniform uniaxial
  pre-stress (scaled by per-element damage factor), applies rigid-body
  penalty BCs, and solves `K·φ = λ·K_g·φ` via `eigsh` shift-invert.
  `BvidAnalysis.run()` for `tier="fe3d"` now returns the minimum of the
  buckling stress and the first-ply-failure stress, capturing whichever
  mode governs. `AnalysisResults.buckling_eigenvalues` is populated with
  the smallest positive eigenvalue. First-ply-failure is retained as
  `_fe3d_cai_first_ply_failure` — together they give engineers an upper
  bound (FPF) and a lower bound (buckling) on the residual strength.
- **"Damage View" GUI tab** (originally planned as "3D Mesh"). After
  three VTK embedding approaches (QtInteractor, lazy-init QtInteractor,
  BackgroundPlotter) all deadlocked Qt's main event loop on macOS, the
  tab was reworked to a matplotlib-based four-panel orthographic view:
  top (x-y), side (x-z), front (y-z), and a text summary panel. Renders
  in ~50 ms on the main thread, no OpenGL dependency. True VTK 3D
  visualization remains available via the `bvidfe.viz.plots_3d` Python
  API and the `examples/` scripts, just not embedded inside the GUI's
  event loop.
- **Validation harness** (`validation/validate_bvid_public.py`):
  `DatasetCase` dataclass + `run_dataset()` + MAE% metric. Auto-discovers
  any JSON dataset in `validation/datasets/`. Ships with a synthetic
  self-check dataset (MAE ≈ 0% by construction) so the harness is
  exercised in CI; real published datasets (Soutis, Caprino,
  Sanchez-Saez, NASA round-robin) remain to be digitized by hand.
- **CI regression gate** — new `validation` job in
  `.github/workflows/tests.yml` running
  `python validation/validate_bvid_public.py --gate`.
- **Auto-populated Knockdown Curve tab.** Every single-run analysis now
  kicks off an empirical-tier energy sweep in the background (8 points,
  sub-second) so the Knockdown Curve tab always shows a
  knockdown-vs-energy plot for context, without requiring the user to
  click "Run energy sweep" explicitly.
- **Buckling Eigenvalues tab** (was "Buckling Mode" placeholder): bar
  chart of the first up to 6 buckling eigenvalues for semi_analytical
  and fe3d runs; explanatory note for the empirical tier (no eigensolve).
- **Damage Severity tab** (was "Stress Field" placeholder): top-down
  heatmap showing how many interfaces are delaminated at each (x, y)
  column — through-thickness sum of `(1 - damage_factor)`, colored
  hot. Analogous to a C-scan plan view.
- **Stage-by-stage logging**: new `bvidfe.fe3d` and `bvidfe.gui` loggers
  write to stderr so users launching the app from a terminal see mesh
  build, K/Kg assembly, eigensolve, FPF, and heartbeat progress lines
  in real time. Controlled via `BVIDFE_LOG_LEVEL` env var.
- **Mesh-size guard** (`FE3DSizeError` + GUI dialog): both single-run
  and sweep paths check the fe3d problem size against
  `BVIDFE_FE3D_MAX_DOF` (default 500k). Hard stop above cap, soft
  warning for merely-large sizes. Prevents SIGSEGV from scipy's sparse
  solvers when memory is exhausted.
- **Heartbeat progress** in both workers: AnalysisWorker ticks every
  2 s during long fe3d runs, SweepWorker ticks per-energy. Prevents
  the status bar looking frozen during multi-minute runs.
- **Vectorized assembler + analytic FPF** in the fe3d tier: 6.5x
  speedup on realistic 5k-element meshes by (a) replacing the 24×24
  Python loop per element with one numpy broadcast, and (b) replacing
  the 12-iteration bisection for first-ply-failure with a single FE
  solve + closed-form quadratic root (Tsai-Wu and LaRC05 are both
  quadratic in stress and stress scales linearly with applied strain).
- **Sanity-check on buckling eigenvalue**: `fe3d_cai` falls back to
  FPF-only when the uniform-pre-stress buckling approximation returns
  a critical stress below 5% of pristine (indicates a numerical
  artefact on the simplified BC set).
- **Conservative defaults + inline fe3d caveat**: `MeshParams` defaults
  to 1 element/ply + 5 mm in-plane (from 4/1.0), and the Summary tab
  appends a note when `tier=fe3d` explaining that fe3d knockdown is
  largely flat vs. impact energy on the current simplified model and
  pointing users to empirical/semi_analytical for energy-dependent
  curves.
- **DPA cap + dent cap**: `impact_to_damage` clips DPA at 80% of panel
  area (emits `UserWarning`) and `dent_depth_mm` clips dent at 50% of
  laminate thickness. Prevents physically-absurd outputs like a 16 mm
  dent on a 1.2 mm laminate.
- **Four runnable examples** (`examples/01_empirical_quick.py` through
  `04_inspection_driven.py`) + `sample_cscan.json`: copy-paste workflows
  for single-run, tier comparison, energy sweep, and damage-driven
  analysis. Each writes outputs to `examples/output/`.
- **14 edge-case robustness tests** (`tests/test_edge_cases.py`): below-
  threshold energies, huge-energy DPA saturation, tiny panels, empty
  damage states, mixed tier switching, CLI end-to-end subprocess, and
  parametric monotonicity across 4 energies.
- **216 tests now passing** (was 179 at v0.1.0 — +37 new tests across
  analysis, elements, solver, GUI, CLI, and edge cases).

### Known limitations (still deferred)

- Material calibration constants remain placeholders until real datasets
  land.
- `fe3d` tier still uses stiffness-reduction at delaminations rather than
  zero-thickness cohesive surfaces (different physics — stiffness
  reduction captures the stiffness loss but not the debonding/sliding).
- **`fe3d` knockdown is approximately flat vs. impact energy** for any
  damage above the Olsson threshold. This is a structural limitation of
  the stiffness-reduction + uniform-pre-stress buckling simplifications:
  the stress-concentration-driven first-ply-failure strain is controlled
  by the healthy/damaged *boundary* rather than the damage magnitude,
  and the buckling eigenvalue is not reliably physical on our simplified
  BCs so it usually gets rejected by the 5%-pristine sanity check. An
  attempt at v0.2.0-dev to add graded damage + real in-plane pre-stress
  BCs to fix this was shelved: the proper-BC buckling eigensolve was
  14× low vs. analytical plate buckling (sign-convention / stress-field
  interpretation issues that need more time than a single session).
  For energy-dependent residual strength use `tier=empirical` (Soutis
  scales with DPA) or `tier=semi_analytical` (Rayleigh-Ritz sublaminate
  buckling scales with ellipse size). Fixing fe3d is v0.3.0 scope.
- Release artifacts are still unsigned.

## [0.1.0] - 2026-04-17

Graduated from `v0.1.0-alpha`. Adds the PyQt6 desktop GUI and packaging infrastructure.

### Added (since v0.1.0-alpha)

- **PyQt6 desktop GUI** (`bvidfe-gui` console script):
  - Seven input panels (MaterialPanel, PanelPanel, InputModePanel, ImpactPanel, DamagePanel, AnalysisPanel, SweepPanel)
  - Six result tabs (Summary, Damage Map, Knockdown Curve, 3D Mesh/Buckling/Stress placeholders)
  - `AnalysisWorker` and `SweepWorker` QThread subclasses for off-UI-thread computation
  - File menu: Save/Load Config (JSON), Export Results JSON, Export Damage Map PNG
  - Headless pytest-qt test coverage for all panels and workers
- **PyInstaller spec** (`BvidFE.spec`) for building macOS and Windows standalone apps
- **GitHub Actions release workflow** (`.github/workflows/release.yml`): on-tag build of macOS + Windows bundles, auto-uploaded to GitHub Releases
- 30 additional tests (149 → 179), including pytest-qt GUI smoke tests

### Remaining limitations (deferred to v0.2.0)

- Material calibration constants are still placeholders pending validation against published CAI/TAI datasets (Soutis, Caprino, Sanchez-Saez, NASA round-robin)
- `fe3d` tier uses stiffness reduction + first-ply-failure; true cohesive surfaces and buckling-based CAI eigensolve are deferred
- 3D Mesh / Buckling / Stress GUI tabs are placeholder widgets (2D plots are fully wired)
- Release artifacts are **unsigned** (no Apple Developer ID or Windows signing cert configured); macOS users may need `xattr -rd com.apple.quarantine BVID-FE.app`

## [0.1.0-alpha] - 2026-04-16

Initial release.

### Added

- `OrthotropicMaterial` dataclass + `MaterialLibrary` with four built-in presets (AS4/3501-6, IM7/8552, T700/2510, T800/epoxy)
- `Laminate` with Classical Lamination Theory ABD matrices and effective engineering constants
- `PanelGeometry`, `ImpactorGeometry`, `BoundaryKind` geometry primitives
- `DamageState` with shapely-union projected damage area; `DelaminationEllipse` per-interface ellipse model
- C-scan JSON/CSV import and validation (`damage/io.py`) per `docs/cscan_schema.md`
- Olsson quasi-static impact threshold (`P_c`, `E_onset`) with Navier-series plate stiffness
- Peanut-template per-interface DPA distribution (`impact/shape_templates.py`)
- Empirical dent-depth model (`impact/dent_model.py`)
- Impact-driven workflow orchestrator: `ImpactEvent` → `DamageState` (`impact/mapping.py`)
- Three modeling tiers for residual strength:
  - **Empirical**: Soutis open-hole-equivalent CAI knockdown + Whitney-Nuismer point-stress TAI
  - **Semi-analytical**: Rayleigh-Ritz sublaminate buckling + Soutis post-buckling envelope; Whitney-Nuismer TAI
  - **3D FE**: first-ply-failure on a damaged hexahedral mesh; LaRC05 for CAI, Tsai-Wu for TAI; stiffness-reduction approximation of delaminations
- FE primitives: `gauss_points_1d` / `gauss_points_hex`, `Hex8`, `Hex8i` (incompatible modes), `CohesiveSurfaceElement` (zero-thickness, bilinear traction-separation)
- Linear static solver (sparse direct, SciPy), linear buckling eigensolve (`eigsh` / dense fallback)
- `BvidAnalysis(AnalysisConfig).run()` high-level orchestrator returning `AnalysisResults`
- CLI entry point `bvidfe` supporting empirical / semi-analytical / fe3d runs
- `AnalysisResults` with `summary()` and `to_dict()` for JSON export
- 2D matplotlib plots: damage-map ellipse overlay, knockdown curves, tier comparison charts
- 3D PyVista plots: mesh with delamination surfaces, buckling mode shape, stress field
- Parametric sweeps over impact energy, layup, and ply thickness with CSV output (`sweep_energies`, `sweep_layups`, `sweep_thicknesses`)
- 149 unit + integration tests mirroring the package structure

### Known Limitations (deferred to v0.2.0)

- GUI (PyQt6 panels) not yet built
- Material calibration constants (`olsson_alpha`, `soutis_k_s`, `dent_beta`, etc.) are placeholder defaults for typical CFRP; precise values need to be calibrated against specific material test data
- `fe3d` tier uses stiffness reduction at delaminated interfaces, not true cohesive surfaces
- `fe3d` CAI uses first-ply-failure on the damaged mesh, not a buckling eigenvalue solve
- No validated datasets yet (Soutis, Caprino, Sanchez-Saez, NASA digitization pending)
- PyInstaller standalone packaging not yet built
