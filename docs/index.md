# BVID-FE

A Python library for predicting residual strength and stiffness of fiber-reinforced composite laminates containing Barely Visible Impact Damage (BVID).

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/elhajjar1/BVID-FE/actions/workflows/tests.yml/badge.svg)](https://github.com/elhajjar1/BVID-FE/actions/workflows/tests.yml)

## Why this tool?

Low-velocity impacts on composite structures can create internal delaminations
that are invisible to the naked eye yet significantly degrade compression and
tension strength. Engineers certifying composite airframes and pressure vessels
need quick estimates of how much strength is lost and how the knockdown
(residual strength divided by undamaged strength) depends on layup, panel
size, and impact energy. BVID-FE provides three modeling tiers — from a
closed-form empirical estimate that runs in tens of milliseconds to a 3D
finite-element solution that takes minutes — so the level of fidelity can be
chosen for each stage of the design process.

BVID-FE is the third in a family of defect-specific composite tools, joining
**PorosityFE** (porosity defects) and **WrinkleFE** (fiber waviness). The three
tools share material models, laminate theory, failure criteria, and
documentation conventions.

## Features

- **Two workflow paths** converging on a shared `DamageState`:
    - *Impact-driven*: Olsson quasi-static threshold + energy-balance damage
      area + peanut-template distribution over the interfaces + empirical
      dent model
    - *Inspection-driven*: C-scan JSON import per the documented schema
      (see [C-Scan Schema](cscan_schema.md))
- **Three modeling tiers** for residual strength after BVID:
    - *Empirical*: Soutis CAI knockdown + Whitney-Nuismer TAI (under a
      second)
    - *Semi-analytical*: buckling-driven delamination growth (Rayleigh-Ritz
      sublaminate buckling + thin-film energy release rate +
      Benzeggagh-Kenane mixed-mode toughness, weakest interface governs),
      capped by Soutis; Whitney-Nuismer for TAI (under a second)
    - *3D FE*: first-ply failure on a hexahedral mesh with the damage zone
      modelled as a soft inclusion, normalised by an undamaged fe3d run;
      LaRC05 for CAI, Tsai-Wu for TAI (minutes)
- **CAI and TAI loading modes** (Compression-After-Impact and
  Tension-After-Impact)
- **Per-interface ellipse damage model** using `DelaminationEllipse` with
  shapely-union projected damage area
- **Pristine (unnotched) laminate strength** from classical lamination
  theory with ply-discount last-ply failure, within 0.83–1.14× of NCAMP
  unnotched tension and compression strengths
- **Six material presets**: AS4/3501-6, AS4/8552, IM7/8552, T700/2510, T800S/3900-2B, T800/epoxy
- **Validation harness** scored against two measured CAI datasets (see
  [Validation](#validation) below)
- **CLI**, **Streamlit web app**, and **Python API**
- **Parametric sweeps** over impact energy, layup, or ply thickness with
  CSV output
- **2D / 3D visualisations** of damage and stress fields

## Installation

```bash
git clone https://github.com/elhajjar1/BVID-FE.git
cd BVID-FE
pip install -e ".[dev]"
pytest -v
```

## Validation

Three measured CAI datasets are in the repository's validation harness;
the first two (ASTM D7136/D7137) are gated in CI:

| Dataset | Material | Coupons | Empirical MAE | Semi-analytical MAE |
| --- | --- | --- | --- | --- |
| NCAMP CAM-RP-2010-002 | AS4/8552 [45/0/-45/90]3s, about 30 J | 7 | 1.5% | 1.5% |
| Lovejoy & Scotti (NASA LaRC 2019) | IM7/8552 quasi-isotropic at 54 J and hard laminate at 68 J | 8 | 9.0% | 9.0% |
| Hasebe et al. 2025 (not gated) | T800S/3900-2B, 8 layups, 3 hemispheres, a flat-ended impactor and a cone, 80 × 50 mm JIS K 7089 coupons, 2–21 J | 343 | 21.3% | 23.3% |

MAE is the mean absolute error between predicted and measured CAI strength.
These errors are in-sample: the Soutis constant `soutis_k_s` = 4.7 and the
default damage-area constant `olsson_alpha` = 0.078 were fitted to these two
datasets. Fitted on one dataset alone, `soutis_k_s` predicts the other with
an MAE of 6.2–9.0%, which is the error to expect on a new material or
layup. On the Hasebe coupons the predicted damage area matches the C-scans,
but the tiers beat a constant knockdown (30.3%) only by 7–9 points: the
damage-area-based knockdown misses the strong effect of layup. There is no measured TAI validation yet. See
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
