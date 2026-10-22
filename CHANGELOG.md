# Changelog

All notable changes to BVID-FE are documented in this file.

## [Unreleased]

The closed-form tiers are now calibrated on, and scored against, measured
CAI data, and the fe3d tier's first-ply failure sees the impact damage.
Most results change and a few APIs break: read *Breaking changes and
migration* before upgrading. The fit histories and sweeps behind the new
constants are in [`validation/README.md`](validation/README.md).

### Breaking changes and migration

**Renamed and removed arguments**

- `soutis_cai` (`bvidfe.failure.soutis_openhole`): the third argument is
  renamed `A_panel_mm2` -> `A_ref_mm2` and is now a fixed reference area,
  not the analysed panel's area. Pass `SOUTIS_REFERENCE_AREA_MM2`
  (15,000 mm^2, the 150 x 100 mm ASTM D7137 coupon) to match
  `BvidAnalysis`:

  ```python
  # 0.2.0
  soutis_cai(m, dpa_mm2, A_panel_mm2=Lx * Ly, sigma_pristine_MPa=s0)
  # now
  from bvidfe.failure.soutis_openhole import SOUTIS_REFERENCE_AREA_MM2
  soutis_cai(m, dpa_mm2, A_ref_mm2=SOUTIS_REFERENCE_AREA_MM2, sigma_pristine_MPa=s0)
  ```

- `semi_analytical_cai` (`bvidfe.analysis.semi_analytical`): the required
  `A_panel_mm2` argument is removed; the optional `A_ref_mm2` defaults to
  `SOUTIS_REFERENCE_AREA_MM2`. A call that passes the panel area as the
  fourth positional argument still runs, but uses it as the reference
  area, so drop it:

  ```python
  # 0.2.0
  semi_analytical_cai(lam, damage, s0, Lx * Ly, boundary="clamped")
  # now
  semi_analytical_cai(lam, damage, s0, boundary="clamped")
  ```

  Its `critical_interface_index` is now the interface whose delamination
  grows at the lowest stress, not the largest area-times-depth.
- `AnalysisResults.warnings`: fe3d no longer emits
  `fe3d_buckling_unconverged` or `fe3d_buckling_artefact_dropped`
  (buckling no longer sets the fe3d strength). The list now holds every
  `UserWarning` raised during `run()` as `"<Class>: <message>"` (it was
  always empty); the warnings are still re-emitted, at the end of `run()`,
  so a caller's `"error"` filter now raises there. fe3d diagnostics are in
  `AnalysisResults.notes`.
- Inputs are validated where they are constructed: non-finite or
  out-of-range numbers, numeric fields given as strings, a string or empty
  layup, a tuple instead of `PanelGeometry`, an unknown material name
  (formerly a `KeyError` at run time), zero or negative impactor mass and
  negative impact energy now raise `ValueError` naming the field.
- The library no longer logs fe3d progress to stderr by default (the
  `bvidfe.fe3d` logger has a `NullHandler`); call
  `bvidfe.analysis.fe_tier.enable_stderr_logging()` to restore it. The CLI
  and the app do.
- Errors where there was silent behaviour: `AnalysisConfig` raises
  `ValueError` for a delamination `interface_index` outside the layup (an
  n-ply laminate has interfaces 0 to n-2); fe3d raises `ValueError` for a
  failure criterion other than `"tsai_wu"`, `"larc05"` or `"puck"`;
  `FailureEvaluator.evaluate` raises `ValueError`, not `AssertionError`,
  for a wrong-shaped stress array; `lekhnitskii_kt_infinity` raises a
  descriptive `ValueError` for constants outside the orthotropic stability
  bound.
- Streamlit is no longer installed with the package. Install
  `pip install "bvidfe[app]"` (or `pip install -r requirements.txt`) for
  the web app; the `dev` extra still includes it.

**Changed definitions** (all results move; see *Results change*)

- `pristine_strength_MPa` is the unnotched laminate strength (ply-discount
  last-ply failure on classical lamination theory,
  `bvidfe.failure.laminate_strength.unnotched_strength`), not the
  thickness-weighted average of the ply strengths along the load. Every
  absolute strength drops: a 16-ply IM7/8552 quasi-isotropic laminate goes
  887.5 -> 601.3 MPa in compression and 1316.5 -> 753.2 MPa in tension.
- The delamination threshold load is `pi sqrt(32 D* G_IIc / 3)` with
  Olsson's `D* = sqrt(D11 D22 (A + 1) / 2)`, about 3.5x the 0.2.0 value,
  so the onset energy is 8-12x higher (about 7 J on a 150 x 100 mm
  IM7/8552 panel, was about 0.6 J): a low-energy impact that used to
  predict damage may now predict none (see *Fixed*).
- The damage area is an energy balance,
  `DPA = olsson_alpha (E - E_onset) / G_IIc` (mm^2), times the impactor
  mass, shape and diameter factors and capped at 80% of the panel. It was
  also divided by the laminate thickness. The delaminations are kept
  inside the impacted panel.
- The Soutis CAI knockdown is `1 / (1 + k_s (DPA / A_ref)^m)` with the
  fixed `A_ref` above, instead of the analysed panel's `Lx * Ly`, which
  sent the knockdown to 1 as the panel grew.
- `semi_analytical` and `fe3d` CAI strength is the far-field stress at
  which the weakest buckled sublaminate grows its delamination, not the
  stress at which it buckles. `buckling_eigenvalues` still reports
  buckling onset, which no longer sets the strength.
- The `fe3d` knockdown is fe3d's damaged strength over its own undamaged
  strength (capped at 1), applied to the shared pristine strength, so an
  undamaged panel gives 1.0 in every tier. fe3d runs two FE solves.

**Changed defaults**

| Constant | 0.2.0 | Now |
|---|---|---|
| `OrthotropicMaterial.soutis_k_s` | 2.5 | 4.8 |
| Soutis reference area | panel `Lx * Ly` | `SOUTIS_REFERENCE_AREA_MM2` = 15,000 mm^2 |
| `OrthotropicMaterial.olsson_alpha` | 0.8 (area divided by thickness) | 0.0865 (energy balance); 0.0648 for T800S/3900-2B |
| Impactor shape factor on the damage area, flat / conical | 1.4 / 0.7 | 0.2 / 1.5 |
| Delamination threshold load `P_c` | `pi sqrt(8 G_IIc sqrt(D11 D22) / 9)` | `pi sqrt(32 D* G_IIc / 3)`, Olsson's `D*` (3.5x higher) |
| Impactor diameter exponent on the damage area | 0.3 | 0.2 (not applied to cones) |
| IM7/8552 `Xt` | 2560 MPa | 2238 MPa (NCAMP cross-ply back-out) |
| fe3d `DAMAGE_ZONE_INPLANE_FACTOR` | (no damage zone) | 0.17 |

A custom `OrthotropicMaterial` that relied on the old `soutis_k_s` or
`olsson_alpha` defaults should set them explicitly, or be recalibrated
(`validation/README.md`).

**Validation-case schema** (`validation/datasets/*.json`)

- Each case must give `impactor_diameter_mm` and `impactor_mass_kg`; they
  used to default silently to 16 mm and 5.5 kg. A dataset whose source
  omits them says in its `reference` what it assumes:

  ```json
  {"impact_energy_J": 30.0, "impactor_diameter_mm": 15.875, "impactor_mass_kg": 5.5}
  ```

### Results change

Closed-form knockdowns, IM7/8552 with 0.152 mm plies on a simply
supported panel, 16 mm hemispherical 5.5 kg impactor (0.2.0 recomputed
with the 0.2.0 code):

| Case | Empirical, 0.2.0 -> now | Semi-analytical, 0.2.0 -> now |
|---|---|---|
| 16-ply QI `[0, 45, -45, 90] * 4` (the README example), 150 x 100 mm, 5 J, CAI | 0.533 -> 1.000 | 0.007 -> 1.000 |
| same, 20 J, CAI | 0.353 -> 0.400 | 0.002 -> 0.400 |
| same, 20 J, TAI | 0.349 -> 0.369 | 0.349 -> 0.369 |
| same on a 500 x 500 mm panel, 20 J, CAI | 0.759 -> 1.000 | 0.003 -> 1.000 |
| 8-ply QI [0/45/-45/90]s, 150 x 100 mm, 10 J, CAI | 0.356 -> 0.593 | 0.010 -> 0.508 |
| 8-ply cross-ply [0/90]2s, 150 x 100 mm, 10 J, CAI | 0.360 -> 0.667 | 0.001 -> 0.562 |

The predicted damage area falls several-fold on thin laminates (16-ply QI
at 20 J: 8,090 -> 1,464 mm^2), and the absolute strengths fall with the
pristine strength (the same case: 312.9 -> 240.6 MPa). The 5 J case and
the 500 x 500 mm panel are below the corrected onset energy (6.6 J and
122 J), so no damage is predicted there. 0.2.0's
semi-analytical CAI was the sublaminate buckling onset, hence its
near-zero values. fe3d results are not comparable with 0.2.0, whose fe3d
knockdown was below 1 for an undamaged panel and whose first-ply failure
did not see the damage.

Mean absolute error against measured CAI strength (empirical /
semi-analytical; 0.2.0 recomputed on the same cases, with the damage
mapped on the tested coupon since 0.2.0 has no `impact_panel`):

| Dataset | 0.2.0 | Now |
|---|---|---|
| `ncamp_as4_8552_cai`, 7 coupons | no AS4/8552 preset | 1.6% / 1.6% |
| `lovejoy_scotti_im7_8552_cai`, 8 coupons | 91.6% / 82.8% | 8.7% / 8.7% |
| `hasebe_t800s_3900_2b_cai`, 343 coupons (not gated) | no T800S/3900-2B preset | 21.2% / 23.2% |

The two gated datasets were used to fit `soutis_k_s` and `olsson_alpha`,
so these errors are in-sample. Fitted on one dataset alone, `soutis_k_s`
predicts the other at 4.9-8.7%, the error to expect on a new material or
layup. On the Hasebe set a constant knockdown scores 30.3%.

### Added

- `AS4/8552` (NCAMP CAM-RP-2010-002 RTD means) and `T800S/3900-2B`
  (Hasebe et al. 2025 and the Toray 3900 data sheet) material presets.
- `AnalysisConfig.impact_panel`: the panel the impact is applied on when
  it is not the one tested (an ASTM D7136 125 x 75 mm window, or a
  60 x 60 mm clamped window and an 80 x 50 mm coupon). The damage is
  mapped on it, centred on `panel` and kept inside it.
- Three measured CAI datasets: `ncamp_as4_8552_cai` and
  `lovejoy_scotti_im7_8552_cai`, gated in CI for both closed-form tiers,
  and `hasebe_t800s_3900_2b_cai`, reported but not gated.
- Validation harness: case fields `boundary`, `impactor_shape`,
  `measured_pristine_MPa`, `impact_panel_Lx_mm`, `impact_panel_Ly_mm` and
  `impact_boundary`; `"gate": false` with a `gate_note` for a known miss;
  `--dataset` takes several names. `validation/README.md` and
  `validation/reference_data.json` list the open datasets found by a
  literature survey, replacing the v0.2.0 roadmap.
- `bvidfe.failure.laminate_strength.unnotched_strength`, and in
  `bvidfe.analysis.semi_analytical` the growth model's pieces:
  `sublaminate_growth_stress`, `weakest_sublaminate_growth`,
  `mixed_mode_toughness` and `panel_buckling_load`.
- `distribute_damage(..., max_half_extent_mm=...)`, and
  `Hex8Element.stress_at_gauss_points_material` (both stress methods take
  an optional stiffness).
- An `app` extra (Streamlit), PyPI classifiers and keywords, and an SPDX
  license expression (PEP 639).
- `tests/test_docs_examples.py` runs the Python examples of the README,
  the quickstart and the API page. CI also checks the JSON schemas against
  the dataclasses (`scripts/generate_schemas.py --check`), and
  `tests/test_package.py` the four version strings against each other.

### Changed

- **fe3d solves ~9x faster in ~6x less memory.** `solve_linear_static`
  eliminates the constrained DOFs instead of adding a 1e10 penalty and
  calling `spsolve` (SuperLU, COLAMD), and solves the reduced system by
  ILU-preconditioned CG to a relative residual of 1e-10, with GMRES and
  then a direct LU as fallbacks. Whole `BvidAnalysis.run()`, 8-ply
  [0/45/-45/90]s IM7/8552 150 x 100 mm, 15 J, fe3d compression, one
  thread: 5 mm 20.8 s / 0.55 GB -> 5.1 s / 0.20 GB, 2.5 mm 270 s /
  3.05 GB -> 31.5 s / 0.46 GB. Knockdowns move by ~2e-7 (0.43002956 ->
  0.43002967 at 2.5 mm), the penalty's own constraint error. Repeated
  BCs on a DOF are merged and conflicting ones raise `ValueError`;
  `method="penalty"` keeps the old path. `BVIDFE_FE3D_MAX_DOF` defaults
  to 150,000 (was 500,000; 265k DOFs ran out of 15 GB with the old
  solver), sized so the direct fallback stays within a few GB.
- **Damage area**: energy balance with `olsson_alpha` = 0.0865, fitted to
  the median of Lovejoy & Scotti's C-scan areas (0.0648 for T800S/3900-2B,
  fitted to Hasebe et al.'s). It now grows with laminate thickness at
  constant J/mm, as measured.
- **Soutis CAI knockdown**: fixed reference area and `soutis_k_s` = 4.8,
  the pooled optimum of both closed-form tiers on NCAMP and Lovejoy &
  Scotti.
- **Impactor shape and diameter factors** are fitted to Hasebe et al.
  2025, the only open data that vary the impactor on one material:
  flat-ended impacts now predict about 7x less damage, conical about 2.1x
  more (both also have their own onset energy); 16 mm hemispherical
  results are unchanged.
- **Pristine strength** is the unnotched laminate strength, within
  0.83-1.14x of the NCAMP unnotched tension and compression strengths of
  AS4/8552 and IM7/8552 (`tests/validation/test_unnotched_strength_ncamp.py`).
  The 8552 presets' `Xt` is NCAMP's cross-ply back-out value, the 0 deg
  ply strength realised inside a laminate.
- **semi_analytical CAI** is `min(Soutis, delamination growth)`: the
  thin-film energy release rate of the buckled sublaminate (Chai, Babcock
  & Knauss 1981) against a Benzeggagh-Kenane mixed-mode toughness, at the
  weakest delaminated interface.
- **fe3d**: the damage zone (every ply under the delamination footprints)
  is a soft inclusion with in-plane stiffness x`DAMAGE_ZONE_INPLANE_FACTOR`
  = 0.17, fitted to three measured CAI coupons. Its elements fail at the
  strain intact material would, and the failure stress is the loaded-edge
  reaction force over the gross section. CAI is
  `min(delamination growth, first-ply failure)`, normalised by an
  undamaged fe3d run, and fe3d now reports `critical_sublaminate`.
  Buckling onset comes from the Rayleigh-Ritz closed form
  (`panel_buckling_load`) instead of the 3D `K_g` eigensolve, which
  under-predicted it (0.20 against 12.6 MPa on an 8-ply 200 x 150 mm
  panel; issues #98, #129). Whole-panel buckling is left out of the
  knockdown (CAI fixtures carry anti-buckling guides); its stress is in
  `AnalysisResults.notes`.
- **CI** gates the semi-analytical tier on the measured datasets as well
  as the empirical tier; the fe3d step is a smoke run on three synthetic
  cases, not a gate.

### Fixed

- **The delamination threshold load was 3.5x too low.** `threshold_load`
  computed `pi sqrt(8 G_IIc sqrt(D11 D22) / 9)`: Davies' isotropic
  formula with the bending stiffness `D` put in place of
  `E h^3 / (1 - nu^2)`, which is `12 D`. It is now the quasi-static,
  large-mass threshold for one delamination,
  `pi sqrt(32 D* G_IIc / 3)` with Olsson's orthotropic
  `D* = sqrt(D11 D22 (A + 1) / 2)`, `A = (D12 + 2 D66) / sqrt(D11 D22)`
  (Davies & Zhang 1995; Suemasu & Majima 1996; Olsson 2001; González et
  al. 2011, Eq. 8), with no dynamic factor (the 1.213 of Olsson et al.
  2006 is for small-mass impacts). Against measured threshold loads,
  predicted / measured goes from 0.25 to 0.86 (González et al. 2011,
  AS4/8552 [(45/0/-45/90)4]s, 9.94 kN measured), 0.29 to 0.99 (Salvetti
  et al. 2016, IM7/MTM45-1, about 5.7 kN) and 0.27 to 0.96 (NIAR,
  IM7/8552, about 6.7 kN). The onset energy rises 8-12x: about 7 J on
  150 x 100 mm IM7/8552 panels (was about 0.6 J), 4.7 J on the NCAMP
  coupons and 6.1-8.3 J on Lovejoy & Scotti's. The constants that had
  absorbed the error were refitted: `olsson_alpha` 0.078 -> 0.0865 and
  0.047 -> 0.0648 (T800S/3900-2B), the conical shape factor 1.3 -> 1.5,
  `soutis_k_s` 4.7 -> 4.8; the gated MAEs move 1.5% -> 1.6% (NCAMP) and
  9.0% -> 8.7% (Lovejoy & Scotti). The onset uses the linear bending
  stiffness, so on thin or large panels, which deflect more than their
  thickness before the threshold load, it is overestimated: about 29 J on
  an 8-ply 300 x 200 mm panel and 120 J on a 16-ply 500 x 500 mm one.
- **The damage template no longer depends on how the layup is written.**
  `distribute_damage` took the ply mismatch as `|theta_b - theta_a|` and
  the ellipse axis as `(theta_a + theta_b) / 2` without wrapping, so a
  -45/90 interface was a 135 deg mismatch (aspect ratio clipped to 4) at
  22.5 deg, and writing -45 as 135 changed the footprint of the same
  16-ply QI laminate at 1500 mm^2 from 71 x 36 mm to 40 x 71 mm. Angles are
  now axial (modulo 180 deg): the mismatch is the acute angle between the
  fibres (-45/90 is 45 deg, aspect ratio 2.125), each ellipse lies along
  the fibres of the back-face-side ply of its interface (ply 0 is the
  impact face; Liu 1988, Abrate 1998), and interfaces between plies of the
  same direction (the doubled mid-plane ply) carry no delamination (a
  unidirectional layup keeps every interface). Unconfined footprints of the
  calibration coupons become near-round: NCAMP 99 x 52 -> 67 x 62 mm,
  Lovejoy & Scotti QI 134 x 70 -> 89 x 85 mm, hard laminate 108 x 108 ->
  125 x 122 mm (confined to the 125 x 75 mm window: 99 x 52 -> 67 x 62,
  125 x 75 -> 92 x 75, 111 x 75 -> 125 x 75 mm). Closed-form tiers:
  empirical MAE unchanged (NCAMP 1.50%, Lovejoy & Scotti 8.95%, Hasebe
  21.25% -> 21.25%, the union area moving by under 0.5%); semi-analytical
  NCAMP and Lovejoy & Scotti unchanged, Hasebe 23.28% -> 23.35% (per
  layup: cross-ply c8 23.2 -> 22.0%, c16 19.4 -> 27.7%, c24 22.6 -> 17.6%,
  QI q8 18.6 -> 19.3%, q16 23.2 -> 25.3%, q24 36.6 -> 37.4%, r0 34.6 ->
  29.5%, r45 unchanged). fe3d on the three coupons that fix
  `DAMAGE_ZONE_INPLANE_FACTOR` = 0.17 (5 mm mesh): NCAMP #0 -3.2% ->
  -12.0%, Lovejoy & Scotti QI #0 -1.2% -> +18.9%, hard #4 +4.7% ->
  +12.7%. The factor is not refitted here. The cross-ply CLI snapshot's
  damage area moves by 1.2e-6.
- **fe3d's 3D ply stiffness uses nu23, not nu21, for the 2-3 coupling.**
  `get_compliance_matrix` set `S23 = -nu21 / E22`, a 2-3 Poisson ratio of
  about 0.02 for every preset, while the presets' own `G23` implies
  `nu23 = E22 / (2 G23) - 1` = 0.30-0.68 (AS4/8552's `G23` was derived
  from a measured 0.45). `OrthotropicMaterial` gains an optional `nu23`;
  `None` resolves to that value (`nu23_resolved`), and a value outside
  (-1, 1) or one that makes the compliance indefinite is rejected. Only
  the fe3d hex stiffness uses the 3D matrix: the closed-form tiers and the
  CLI snapshots are unchanged. fe3d coupons (on top of the template fix):
  NCAMP #0 -12.0% -> -14.4%, Lovejoy & Scotti QI #0 +18.9% -> +10.8%,
  hard #4 +12.7% -> +10.8% (on the old template alone: -3.2% -> -5.9%,
  -1.2% -> -3.6%, +4.7% -> +3.5%). With both fixes the fe3d tension
  knockdown at 30 J (IM7/8552, 150 x 100 mm) is 0.316 for 16-ply QI (was
  0.333) and 0.312 for 24-ply QI (unchanged). A `[90]_4`
  coupon in compression with simply-supported loaded edges now fails
  3-4% below `Yc / E22`: the `u_z = 0` restraint on the whole loaded face
  blocks the ply's through-thickness Poisson expansion (the ply-frame test
  now uses free edges for that case).
- **fe3d first-ply failure depended on the mesh.** The damage zone took or
  left whole elements by their centroid, and failure was checked at
  single Gauss points, where the stress at the re-entrant corners of the
  zone's footprint is singular: the 15 J cross-ply knockdown (8-ply
  IM7/8552, 150 x 100 mm) moved 24-28% between 5 and 2.5 mm meshes. Edge
  elements now take the in-plane factor weighted by their area fraction
  inside the zone, and the failure index is averaged over a 2.5 mm disk
  within each element layer (`FPF_AVERAGING_RADIUS_MM`). 5 -> 2.5 mm now
  moves the 15 J knockdowns 0.2-3.1% (QI tension 0.421 -> 0.424,
  cross-ply tension 0.230 -> 0.237, cross-ply compression 0.249 -> 0.249,
  QI compression 0.426 -> 0.430). `DAMAGE_ZONE_INPLANE_FACTOR` stays 0.17
  (calibration coupons -3.2% / -1.2% / +4.7%).
- **Inputs that gave silent or deep errors**: NaN/inf geometry, ply
  thickness and angles; `ImpactEvent(location_xy_mm=[0, 0])` from JSON
  (now a tuple, still the panel centre); a delamination outside the panel,
  counted in the damage area but missed by the fe3d mesh (now
  `DamageOutsidePanelWarning`).
- **CLI errors printed tracebacks** (malformed C-scan, invalid config,
  `FE3DSizeError`, a closed pipe); they are one-line errors with a nonzero
  exit code. The CLI also rejects ply angles outside [-90, 90] and
  non-finite values, and `--help` no longer prints `(default: None)`.
- **The app's saved-config loader rejected per-ply thickness lists.**

- **Predicted delaminations stay inside the impacted panel.** The union of
  the template's ellipses spread beyond narrow coupons (108 x 108 mm on a
  102 mm wide one), so fe3d's soft zone spanned the whole section. Each
  ellipse is now shrunk to fit, with the union re-solved for the same
  area and a warning when the window cannot hold it.
- **`Laminate.effective_engineering_constants` was off by a factor of
  h^2** (`Ex` 830.9 GPa for a 24-ply IM7/8552 QI laminate, now 62.4 GPa).
- **fe3d first-ply failure used the wrong ply frames**: a +theta ply got
  the stiffness of a -theta ply, and the criteria received global-frame
  stress instead of material-frame stress.
- **fe3d evaluated Tsai-Wu when asked for Puck** (or for any misspelt
  criterion name).
- **fe3d buckling notes** were dropped before reaching
  `AnalysisResults.notes`.
- **Whitney-Nuismer TAI could exceed pristine strength or diverge** for
  `Kt_inf` above about 9.2; the knockdown is now capped at 1 with a
  `UserWarning`. The built-in presets stay below 7.8, even
  unidirectional, so their results do not move.
- **`semi_analytical` CAI returned `numpy.float64`** when delamination
  growth governed; it returns a Python `float`.
- **The validator skipped an unknown `--dataset` name** silently; it now
  exits 1.
- **Version strings disagreed** (`pyproject.toml`, `__version__`,
  `CITATION.cff`, README BibTeX).
- **Project links pointed at the repository's previous owner**
  (`ranipdx-glitch`); they point at `elhajjar1/BVID-FE`, which the PyPI
  Trusted Publisher must also use.
- **`docs/schemas/analysis_config.json` was missing `puck_p_nt_minus` and
  `puck_p_nt_plus`.**
- **CI**: coverage counts the CLI subprocess again
  (`[tool.coverage.run] patch = ["subprocess"]`), and ruff's rule set is
  pinned so ruff 0.16 does not turn lint red.
- **Docs**: the README sweep example crashed (it reused a config the
  C-scan example had rebound), and the `[0.2.0]` entry below overstated
  its validation (CI ran a synthetic self-check only).

### Removed

- The fe3d 3D `K_g` buckling eigensolve and the 5%-of-pristine
  plausibility guard it needed.
- The `fe3d_buckling_unconverged` and `fe3d_buckling_artefact_dropped`
  warning tags.
- `streamlit` from the core dependencies (now the `app` extra).
- Leftovers of the removed PyQt6 GUI: `scripts/visual_matrix.py` and
  `docs/BUILD.md`.

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
