# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Repository Overview

This repository contains **BVID-FE**, a Python package (`bvidfe`, src/-layout) for predicting residual strength and stiffness of fiber-reinforced composite laminates containing Barely Visible Impact Damage (BVID). The tool is the third in a family of defect-specific composite tools alongside PorosityFE and WrinkleFE. It provides two workflow paths (impact-driven via Olsson threshold + peanut-template DPA, and inspection-driven via C-scan JSON import), three modeling tiers (empirical, semi-analytical, 3D FE), and both CAI and TAI loading modes. The package uses a `src/` layout, is MIT-licensed, and targets Python 3.9+.

## File Structure

```
BVID-FE/
├── src/
│   └── bvidfe/
│       ├── __init__.py
│       ├── cli.py                    # bvidfe entry point
│       ├── core/
│       │   ├── material.py           # OrthotropicMaterial + MaterialLibrary
│       │   ├── laminate.py           # Laminate, Ply, CLT ABD matrices
│       │   └── geometry.py           # PanelGeometry, ImpactorGeometry, BoundaryKind
│       ├── impact/
│       │   ├── olsson.py             # Quasi-static damage threshold (P_c, E_onset)
│       │   ├── shape_templates.py    # Per-interface peanut DPA templates
│       │   ├── dent_model.py         # Empirical dent-depth model
│       │   └── mapping.py            # ImpactEvent -> DamageState orchestrator
│       ├── damage/
│       │   ├── state.py              # DamageState, DelaminationEllipse
│       │   └── io.py                 # C-scan JSON import and validation
│       ├── elements/
│       │   ├── hex8.py               # Standard Hex8 element
│       │   ├── hex8i.py              # Incompatible-modes Hex8i element
│       │   ├── gauss.py              # Gauss quadrature (1D and hex)
│       │   └── cohesive.py           # Zero-thickness cohesive surface element
│       ├── solver/
│       │   ├── static.py             # Linear static solver
│       │   ├── assembler.py          # Global stiffness assembly
│       │   ├── boundary.py           # Dirichlet BC handler
│       │   └── buckling.py           # Linear buckling eigensolve
│       ├── failure/
│       │   ├── larc05.py             # LaRC05 (Hashin-3D reduction)
│       │   ├── tsai_wu.py            # 3D Tsai-Wu criterion
│       │   ├── soutis_openhole.py    # Soutis CAI + Whitney-Nuismer TAI
│       │   └── evaluator.py          # FailureEvaluator, LaminateFailureReport
│       ├── analysis/
│       │   ├── config.py             # AnalysisConfig dataclass
│       │   ├── bvid.py               # BvidAnalysis orchestrator
│       │   ├── results.py            # AnalysisResults dataclass
│       │   ├── semi_analytical.py    # Semi-analytical tier
│       │   ├── fe_tier.py            # 3D FE tier
│       │   └── fe_mesh.py            # Damaged hex mesh builder
│       ├── viz/
│       │   ├── plots_2d.py           # Damage map, knockdown curves (matplotlib)
│       │   ├── plots_3d.py           # 3D PyVista plots
│       │   └── style.py              # Publication styling
│       └── sweep/
│           └── parametric_sweep.py   # sweep_energies, sweep_layups, sweep_thicknesses
├── tests/
│   ├── core/          test_material.py, test_laminate.py, test_geometry.py
│   ├── impact/        test_olsson.py, test_shape_templates.py, test_dent_model.py, test_mapping.py
│   ├── damage/        test_state.py, test_io.py
│   ├── elements/      test_hex8.py, test_hex8i.py, test_gauss.py, test_cohesive.py
│   ├── solver/        test_static.py, test_assembler.py, test_boundary.py, test_buckling.py
│   ├── failure/       test_larc05.py, test_tsai_wu.py, test_soutis.py, test_evaluator.py
│   ├── analysis/      test_empirical_path.py, test_semi_analytical_path.py,
│   │                  test_fe3d_path.py, test_fe_mesh.py, test_fe_tier.py,
│   │                  test_semi_analytical.py
│   ├── sweep/         test_sweep.py
│   ├── viz/           test_plots_2d.py, test_plots_3d.py
│   ├── test_cli.py
│   └── test_package.py
├── docs/
│   ├── cscan_schema.md
│   └── superpowers/
│       ├── specs/2026-04-16-bvid-fe-design.md
│       └── plans/2026-04-16-bvid-fe-implementation.md
├── pyproject.toml
├── requirements.txt
├── README.md
├── ARCHITECTURE.md
├── CLAUDE.md
├── CITATION.cff
├── LICENSE
├── CHANGELOG.md
└── CONTRIBUTING.md
```

## Running the Tool

```bash
# Single analysis run (CLI)
bvidfe --material IM7/8552 --layup "0,45,-45,90,90,-45,45,0" \
       --thickness 0.152 --panel 150x100 \
       --loading compression --energy 30

# Run all tests
pytest tests/ -v

# Run a specific test module
pytest tests/impact/test_olsson.py -v

# Run tests for a subsystem
pytest tests/failure/ -v
```

## Code Architecture

`BvidAnalysis(AnalysisConfig).run()` in `analysis/bvid.py` is the single entry point for all analysis paths. It:

1. Resolves the material (string preset or `OrthotropicMaterial` instance) using `MaterialLibrary`.
2. Builds a `Laminate` from the layup and ply thickness via CLT.
3. If `config.impact` is provided, calls `impact/mapping.py` to produce a `DamageState` (Olsson threshold → peanut DPA → dent model).
4. If `config.damage` is provided, uses it directly (inspection-driven path).
5. Dispatches to the selected tier (`empirical`, `semi_analytical`, or `fe3d`).
6. Returns an `AnalysisResults` with `residual_strength_MPa`, `knockdown`, `damage`, and optional tier-specific outputs.

### Tier dispatch

- **empirical**: `failure/soutis_openhole.py` — Soutis open-hole-equivalent knockdown (CAI) or Whitney-Nuismer point-stress criterion (TAI). No mesh, no eigensolve.
- **semi_analytical**: `analysis/semi_analytical.py` — Rayleigh-Ritz plate buckling eigenvalue for the damaged sublaminate; Soutis envelope for post-buckling CAI failure. Whitney-Nuismer for TAI.
- **fe3d**: `analysis/fe_tier.py` + `analysis/fe_mesh.py` — builds a structured hexahedral mesh (`elements/`), reduces stiffness at delaminated interfaces, assembles and solves the global system (`solver/`), and evaluates LaRC05 (CAI) or Tsai-Wu (TAI) failure at all Gauss points (`failure/`).

### DamageState

`damage/state.py` is the universal handoff between the impact-mapping stage and the residual-strength engine. `DelaminationEllipse` objects carry `interface_index`, centroid, semi-axes, and orientation. `DamageState.projected_damage_area_mm2` computes the shapely union of all ellipse footprints.

### Impact mapping

`impact/mapping.py` runs three stages: (1) Olsson threshold to check whether the energy exceeds `E_onset`, (2) total DPA from energy above threshold, and (3) per-interface distribution using `shape_templates.py` peanut templates.

## Key Mathematical Relationships

**Olsson quasi-static damage threshold:**
```
P_c = pi * sqrt(8 * G_IIc * D_eff / 9)
E_onset = P_c^2 / (2 * k_cb)
k_cb = 1 / (1/k_bending + 1/k_hertz)
```
where `D_eff = sqrt(D11 * D22)` from the CLT D-matrix.

**Total projected damage area (above threshold):**
```
DPA = alpha * (E_impact - E_onset) / (G_IIc * h)
```
`alpha` (`olsson_alpha`) is a material-calibration constant (rough default for CFRP; needs test calibration).

**Soutis CAI knockdown (open-hole analogy):**
```
sigma_CAI = sigma_0 * (1 - sqrt(DPA / (pi * a^2))) * k_s
```
`a` is the equivalent circular damage radius; `k_s` (`soutis_k_s`) is a stress-concentration correction factor.

**Whitney-Nuismer TAI (point-stress):**
```
sigma_TAI = sigma_UTS * (1 - (a / (a + d_0))^2)
```
`d_0` is the characteristic distance; `a` is the damage-equivalent hole radius.

**Rayleigh-Ritz sublaminate buckling (semi_analytical tier):**
```
N_cr = pi^2 / b^2 * (D11/m^2 * (b/a)^2 + 2*(D12 + 2*D66) + D22*m^2*(a/b)^2)
```
where `a`, `b` are sublaminate dimensions and `m` is the half-wave count.

**3D Tsai-Wu failure:**
```
F_i * sigma_i + F_ij * sigma_i * sigma_j = 1   (i, j = 1..6, Voigt notation)
```

**LaRC05 (minimal Hashin-3D reduction) — fiber compression mode:**
```
(tau_12 / S12)^2 + (sigma_22 / Yc)^2 = 1   (matrix compression)
```

## Testing Strategy

Tests mirror the package structure: each submodule has a corresponding test file in the matching subdirectory under `tests/`. There are 149 tests in total.

- **Unit tests** cover individual functions (e.g., `test_olsson.py` checks `threshold_load` against hand-calculated values; `test_tsai_wu.py` tests uniaxial and biaxial failure envelopes).
- **Integration tests** exercise full tier paths on small meshes: `test_empirical_path.py`, `test_semi_analytical_path.py`, `test_fe3d_path.py` each build a minimal `AnalysisConfig`, call `BvidAnalysis.run()`, and assert that `knockdown` is in `(0, 1]`.
- **End-to-end tests** in `test_cli.py` invoke the `bvidfe` CLI via `subprocess` and check for clean exit and expected fields in stdout.

Run the full suite:
```bash
pytest tests/ -v --tb=short
```

## Modifying the Models

**To add a new material preset:**
Add an entry to `MATERIAL_PRESETS` in `src/bvidfe/core/material.py` using the `OrthotropicMaterial` dataclass. Include all orthotropic constants (MPa, kg/mm^3), strength values, and interlaminar fracture toughnesses. Add the name to `MaterialLibrary.NAMES` and include a source reference in the docstring.

**To add a new modeling tier:**
1. Create `src/bvidfe/analysis/my_tier.py` implementing a function `run_my_tier(config, damage_state, laminate) -> AnalysisResults`.
2. Register the tier string literal in `AnalysisConfig.tier` (update the `Literal` type).
3. Add the dispatch branch in `BvidAnalysis.run()` in `analysis/bvid.py`.
4. Add tests in `tests/analysis/test_my_tier_path.py`.

**To add a new failure criterion:**
Create `src/bvidfe/failure/my_criterion.py` subclassing `FailureCriterion` from `failure/base.py`. Register it in `FailureEvaluator` in `failure/evaluator.py`.

**To adjust calibration constants:**
All empirical calibration constants are defined in the relevant module (e.g., `olsson_alpha` in `impact/olsson.py`, `soutis_k_s` in `failure/soutis_openhole.py`, `dent_beta` in `impact/dent_model.py`). The defaults are reasonable for typical CFRP but should be calibrated against test data for production use.

## Scientific References

1. Olsson, R. (2001). Analytical prediction of large mass impact damage in composite laminates. *Composites Part A*, 32(9), 1207-1215.
2. Soutis, C. (1996). Compressive strength of unidirectional composites: measurement and prediction. *ASTM STP*, 1242, 168-176.
3. Whitney, J.M. & Nuismer, R.J. (1974). Stress fracture criteria for laminated composites containing stress concentrations. *Journal of Composite Materials*, 8(3), 253-265.
4. Tsai, S.W. & Wu, E.M. (1971). A general theory of strength for anisotropic materials. *Journal of Composite Materials*, 5(1), 58-80.
5. Davila, C.G., Camanho, P.P., & Rose, C.A. (2005). Failure criteria for FRP laminates. NASA/TM-2005-213530 (LaRC05).
6. Reddy, J.N. (2004). *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*. CRC Press.
