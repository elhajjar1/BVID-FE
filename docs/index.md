# BVID-FE

Residual compression (CAI) and tension (TAI) strength of fibre-reinforced
composite laminates with Barely Visible Impact Damage (BVID), from an
impact energy or a C-scan, in three modelling tiers.

## Why this tool?

Low-velocity impacts on composite structures can create internal delaminations
that are invisible to the naked eye yet significantly degrade compression and
tension strength. Engineers certifying composite airframes and pressure vessels
need quick estimates of how much strength is lost and how the knockdown
(residual strength divided by undamaged strength) depends on layup, panel
size, and impact energy. BVID-FE provides three modeling tiers — from
closed-form estimates that run in tens of milliseconds to a 3D
finite-element solution that takes tens of seconds — so the level of
fidelity can be chosen for each stage of the design process.

BVID-FE is the third in a family of defect-specific composite tools, joining
**PorosityFE** (porosity defects) and **WrinkleFE** (fiber waviness). The three
tools share material models, laminate theory, failure criteria, and
documentation conventions.

## Features

- **Two workflow paths** converging on a shared `DamageState`:
    - *Impact-driven*: Olsson quasi-static threshold + energy-balance damage
      area + peanut-template distribution over the interfaces + empirical
      dent model
    - *Inspection-driven*: C-scan JSON import per the documented
      [schema](cscan_schema.md)
- **Three modeling tiers** for residual strength after BVID:
    - *Empirical*: Soutis CAI knockdown + Whitney-Nuismer TAI on an
      equivalent hole for the fibre damage (tens of milliseconds)
    - *Semi-analytical*: buckling-driven delamination growth (Rayleigh-Ritz
      sublaminate buckling + the energy release rate of the post-buckled
      elliptical film (2D blister) + Benzeggagh-Kenane mixed-mode
      toughness, weakest interface governs), capped by Soutis; the
      empirical TAI model (about a second)
    - *3D FE*: first-ply failure on a hexahedral mesh with the damage zone
      modelled as a soft inclusion (fibres intact in tension) and, in
      tension, the closed-form fibre-damage hole as a notch at the impact
      site, normalised by an undamaged fe3d run; LaRC05 for CAI, Tsai-Wu
      for TAI (about 4 s for an 8-ply 150 × 100 mm panel at the default
      5 mm mesh, about 20 s at 2.5 mm; see
      [Tier choice](python_api.md#tier-choice))
- **CAI and TAI loading modes** (Compression-After-Impact and
  Tension-After-Impact)
- **Pristine (unnotched) laminate strength** from classical lamination
  theory with ply-discount last-ply failure, within 0.83–1.14× of NCAMP
  unnotched tension and compression strengths
- **Per-interface ellipse damage model** using `DelaminationEllipse` with
  shapely-union projected damage area
- **Six material presets**: AS4/3501-6, AS4/8552, IM7/8552, T700/2510,
  T800S/3900-2B, T800/epoxy
- **Validation harness** (`validation/`) that scores the tiers against
  four measured CAI datasets (NCAMP AS4/8552 and Lovejoy & Scotti
  IM7/8552, gated in CI; Hasebe et al. T800S/3900-2B and González's
  AS4/8552, reported only) and
  one TAI condition (Körbelin et al. M21/T800S, reported only); see
  [Validation](#validation) below
- **CLI** for single runs, several tiers in one call
  (`--tier empirical,semi_analytical`), and shell pipelines (`--quick`,
  `--quick-json`); see the [CLI reference](python_api.md#cli)
- **Streamlit web app** with sidebar-driven `AnalysisConfig` and result
  tabs — Summary, Damage Map, Knockdown Curve (live empirical sweep around
  the current energy), 3D Damage (Plotly hex mesh coloured by damage
  factor), Buckling (tier-specific buckling indicator: `N_cr` in N/mm for
  `semi_analytical`, `σ_crit` in MPa for `fe3d`), Damage Severity
  (through-thickness sum of per-element `1 − damage_factor`; see
  [Damage severity heatmap](physics_models.md#damage-severity-heatmap)),
  Sweep (CSV-exportable parametric energy sweep), and Export NCR (a
  nonconformance report in Markdown or JSON with a *proposed* disposition
  that a Material Review Board must review). Cached runs, fe3d mesh-size
  guards, JSON/CSV downloads; see [Deployment](deployment.md)
- **Parametric sweeps** over impact energy, layup, or ply thickness with
  CSV output
- **2D plots**: damage map, knockdown curves, tier comparison charts
- **Plotly 3D viz** for the damaged hex mesh (`bvidfe.viz.plotly_3d`)

## Installation

BVID-FE is not on PyPI yet: install the library and the `bvidfe` command
with `pip install git+https://github.com/elhajjar1/BVID-FE.git`, or clone
the repository for the examples, the notebook, the validation harness and
the Streamlit app. [Quickstart → Install](quickstart.md#1-install) has the
commands.

## Validation

Four measured CAI datasets are in the repository's validation harness;
the first two (ASTM D7136/D7137) are gated in CI:

| Dataset | Material | Coupons | Empirical MAE | Semi-analytical MAE |
| --- | --- | --- | --- | --- |
| NCAMP CAM-RP-2010-002 | AS4/8552 [45/0/-45/90]3s, about 30 J | 7 | 1.5% | 1.5% |
| Lovejoy & Scotti (NASA LaRC 2019) | IM7/8552 quasi-isotropic at 54 J and hard laminate at 68 J | 8 | 8.8% | 8.8% |
| Hasebe et al. 2025 (not gated) | T800S/3900-2B, 8 layups, 3 hemispheres, a flat-ended impactor and a cone, 80 × 50 mm JIS K 7089 coupons, 2–21 J | 343 | 21.5% | 23.4% |
| González 2010 thesis, Girona (not gated) | AS4/8552 quasi-isotropic (unclustered and clustered, 32 and 48 plies) and cross-ply, ASTM D7136/D7137, 19–58 J | 19 | 14.3% | 14.4% |

MAE is the mean absolute error between predicted and measured CAI strength.
For NCAMP and Lovejoy & Scotti the predicted strength is the knockdown
times the model's pristine strength, and the two errors partly cancel:
against NCAMP's measured unnotched compression the knockdown MAE is 3.6%
and 13.6% (the hard laminate's pristine strength is 11.9% high, its
knockdown 17.3% low).
These errors are in-sample: the Soutis constant `soutis_k_s` = 4.8 and the
default damage-area constant `olsson_alpha` = 0.0855 were fitted to these two
datasets. Fitted on one dataset alone, `soutis_k_s` predicts the other with
an MAE of 5.0–8.8%, which is the error to expect on a new material or
layup. On the Hasebe coupons the predicted damage area matches the C-scans,
but the tiers beat a constant knockdown (30.3%) only by 7–9 points: the
damage-area-based knockdown misses the strong effect of layup. TAI is not
validated: its one constant is fitted in-sample to a single measured
condition (Körbelin et al. 2020, M21/T800S at 25 J). See
[`validation/README.md`](https://github.com/elhajjar1/BVID-FE/blob/main/validation/README.md)
for the data sources and the calibration.

## Where to next

- [Quickstart](quickstart.md) — run your first analysis in under a minute.
- [Python API](python_api.md) — script BVID-FE workflows beyond the GUI/CLI.
- [C-Scan Schema](cscan_schema.md) — load inspection data from JSON.
- [Physics Models](physics_models.md) — what each tier solves.
- [Deployment](deployment.md) — publish the Streamlit app to a public URL.

## JSON Schemas

Auto-generated JSON Schemas (draft 2020-12) for the public data contracts are
published with this site:

- [`analysis_config.json`](schemas/analysis_config.json) — `AnalysisConfig`
- [`analysis_results.json`](schemas/analysis_results.json) — `AnalysisResults`

The schemas are regenerated from the in-tree dataclasses by
`scripts/generate_schemas.py`.

## Citation

If you use BVID-FE in your research, please cite the project (see the
[CITATION.cff](https://github.com/elhajjar1/BVID-FE/blob/main/CITATION.cff)
at the repository root).

## License

MIT License. See
[LICENSE](https://github.com/elhajjar1/BVID-FE/blob/main/LICENSE) for
details.
