# BVID-FE Architecture

## Module Dependency Diagram

```
core, damage, solver          (no internal dependencies)
  │
  ├─► impact     (core, damage)
  ├─► elements   (core)
  ├─► failure    (core)
  ▼
analysis         (core, damage, elements, failure, impact, solver)
  │
  ├─► viz        (analysis, core, damage)
  ├─► sweep      (analysis, impact)
  ▼
cli.py, app.py   (user-facing entry points)
```

`core`, `damage` and `solver` import nothing else from `bvidfe` (`solver`
works on any element objects that provide a stiffness matrix). `analysis`
is the only package that pulls the physics layers together.
The two user-facing entry points are the `bvidfe` CLI (`src/bvidfe/cli.py`)
and the Streamlit web app (`app.py` at the repo root). Both run analyses
through `analysis.BvidAnalysis`; they also import `core`, `damage` and
`impact` directly to build inputs, and `app.py` also uses `sweep` and `viz`
for its sweep tab and plots. The shared `Literal` aliases live in
`bvidfe/_types.py`.

## Module Catalog

| Module | Submodule | Role |
|--------|-----------|------|
| `_types.py` | — | Shared `Literal` aliases `TierName`, `LoadingMode`, `CriterionName` and their runtime-validation frozensets (`_TIER_NAMES`, `_LOADING_MODES`, `_CRITERION_NAMES`). Each alias and its frozenset change together. |
| `core/` | `material.py` | `OrthotropicMaterial` dataclass (moduli, strengths, toughnesses, density in kg/mm^3, impact-mapping and knockdown calibration constants such as `soutis_k_s` = 4.7 and `olsson_alpha` = 0.078) + the `MATERIAL_LIBRARY` dict with six presets (AS4/3501-6, AS4/8552, IM7/8552, T700/2510, T800S/3900-2B, T800/epoxy). Ported from WrinkleFE. |
| | `laminate.py` | `Laminate`: Classical Lamination Theory ABD matrices, effective engineering constants, flexural rigidity; uniform or per-ply thicknesses. Ported from WrinkleFE. |
| | `geometry.py` | `PanelGeometry`, `ImpactorGeometry`, `BoundaryKind` (`clamped` / `simply_supported` / `free`), `ImpactorShape` (`hemispherical` / `flat` / `conical`). |
| `impact/` | `olsson.py` | Olsson quasi-static damage-threshold load (`threshold_load`) and onset energy (`onset_energy`) from plate bending + fracture energy balance; boundary- and impactor-shape-aware. |
| | `shape_templates.py` | Layup-dependent per-interface "peanut" templates that distribute the total DPA into per-interface elliptical delaminations (`distribute_damage`). With `max_half_extent_mm` each ellipse is shrunk, aspect ratio kept, to fit a window about the impact point, and the union is re-solved for the same area. |
| | `dent_model.py` | Thickness-normalized empirical dent depth (`dent_depth_mm`) and fiber-break core radius (`fiber_break_radius_mm`) from impact energy. |
| | `mapping.py` | `ImpactEvent` and `impact_to_damage`: `ImpactEvent` → `DamageState`. The DPA is an energy balance, `olsson_alpha (E − E_onset) / G_IIc`, with impactor-mass, shape and diameter corrections and an 80%-of-panel cap. The delaminations are kept inside the impacted panel (`impact_panel`, else `panel`; for an ASTM D7136 coupon the 125 × 75 mm support window), measured from the impact point. |
| `damage/` | `state.py` | `DamageState` dataclass + `DelaminationEllipse`; shapely-union projected damage area. |
| | `io.py` | C-scan JSON import/export and validation (`load_cscan_json`, `save_cscan_json`, `damage_state_from_dict`, `CScanSchemaError`) per `docs/cscan_schema.md`. |
| `elements/` | `hex8.py` | 8-node isoparametric hexahedral element: shape functions, B-matrix, 24x24 stiffness, ply-frame stress (`stress_at_gauss_points_material`); `DegenerateElementError` on non-positive `det(J)`. Used by the `fe3d` tier. |
| | `hex8i.py` | Incompatible-modes hex element (Wilson modes) for improved bending accuracy. Not used by the `fe3d` tier. |
| | `gauss.py` | Gauss-Legendre quadrature points and weights: `gauss_points_1d`, `gauss_points_hex`. |
| | `cohesive.py` | Zero-thickness cohesive surface element with a bilinear traction-separation law (pointwise constitutive model only). Not used by the `fe3d` tier, which models delaminations by stiffness reduction. |
| `solver/` | `static.py` | `solve_linear_static`: sparse-direct linear static solve (SciPy). |
| | `assembler.py` | COO → CSC global stiffness matrix assembly from element contributions. |
| | `boundary.py` | `BoundaryCondition`, `apply_dirichlet_penalty`, `uniaxial_x_bcs`: penalty-method Dirichlet BCs for displacement-controlled compression and tension. |
| | `buckling.py` | Generic linear eigenvalue buckling utility (`scipy.sparse.linalg.eigsh`). Not called from the `fe3d` tier since #129; kept for future plate/shell-element buckling work. |
| `failure/` | `laminate_strength.py` | `unnotched_strength(lam, loading)`: unnotched laminate strength by ply-discount last-ply failure on CLT (max-stress per ply, cracked plies keep 1% of E22 and G12). This is `pristine_strength_MPa` for every tier. |
| | `larc05.py` | Minimal LaRC05 composite failure criterion (Hashin-3D reduction); scalar + batch index. |
| | `tsai_wu.py` | Full 3D Tsai-Wu failure criterion with interaction terms; scalar + batch index. |
| | `puck.py` | Puck action-plane inter-fiber failure + max-stress fiber failure; scalar + batch index. |
| | `soutis_openhole.py` | Soutis CAI knockdown (`soutis_cai`; DPA scaled by the fixed 150 × 100 mm ASTM D7137 coupon area `SOUTIS_REFERENCE_AREA_MM2`), Whitney-Nuismer point-stress TAI (`whitney_nuismer_tai`) and Lekhnitskii `Kt_inf` (`lekhnitskii_kt_infinity`). |
| | `evaluator.py` | `FailureEvaluator`: applies a registered criterion across an (element, Gauss point) stress field; produces `LaminateFailureReport`. Standalone utility; the `fe3d` tier calls the criterion functions directly. |
| `analysis/` | `config.py` | `AnalysisConfig` + `MeshParams` dataclasses (material, layup, panel, loading, tier, impact or damage; mesh resolution for fe3d). `MeshParams` validates positive `elements_per_ply` / `in_plane_size_mm` / `cohesive_zone_factor` at construction. |
| | `bvid.py` | `BvidAnalysis(AnalysisConfig).run()`: main orchestrator. Resolves the material and damage, computes the pristine strength (`unnotched_strength`), dispatches to the tier, and sets `knockdown = residual / pristine`. For fe3d it runs the tier twice (damaged and undamaged) and applies the ratio to the pristine strength. |
| | `results.py` | `AnalysisResults` dataclass (`residual_strength_MPa`, `pristine_strength_MPa`, `knockdown`, `buckling_eigenvalues`, `critical_sublaminate`, `notes`, `warnings`, `field_results`, ...) with `summary()` and `to_dict()`. `FieldResults` is reserved; `field_results` is always `None` today. |
| | `semi_analytical.py` | Semi-analytical tier. CAI: min of the Soutis knockdown and the far-field stress at which the weakest delamination grows by sublaminate buckling (`weakest_sublaminate_growth`, `sublaminate_growth_stress`: thin-film energy release rate against a Benzeggagh-Kenane mixed-mode toughness). Also the Rayleigh-Ritz buckling kernels (`sublaminate_buckling_load`, `panel_buckling_load`). TAI: Whitney-Nuismer. |
| | `fe_tier.py` | 3D FE tier. `fe3d_cai` = min(delamination growth from `weakest_sublaminate_growth`, first-ply failure on the damaged mesh); `fe3d_tai` = first-ply failure (Tsai-Wu). First-ply failure uses one static solve with analytic load scaling per criterion, strain-equivalent stresses in the damage zone (`_failure_stress_material`), and the reaction force on the loaded edge over the gross section. `fe3d_cai_buckling` (Rayleigh-Ritz closed form, #129) gives buckling onset for `buckling_eigenvalues` only. |
| | `fe_mesh.py` | Damaged hexahedral mesh construction from `DamageState`. Carries two per-element factors: `damage_factors` (out-of-plane, `DAMAGE_OOP_FACTOR` inside delamination footprints) and `in_plane_damage_factors` (`DAMAGE_ZONE_INPLANE_FACTOR` = 0.17, fitted to three measured CAI coupons, through the thickness over the damage zone, a soft inclusion; `DAMAGE_FIBER_BREAK_INPLANE_FACTOR` inside the optional fiber-break core). |
| `viz/` | `plots_2d.py` | Damage-map ellipse overlay (with impact location + fiber-break-core markers), knockdown curves, per-tier comparison charts (matplotlib). |
| | `plotly_3d.py` | Plotly 3D figures of the hex mesh coloured by damage factor (`mesh_damage_figure`) or any per-element scalar (`stress_field_figure`), for the Streamlit app. |
| | `style.py` | Publication styling constants (fonts, DPI, color maps). |
| `sweep/` | `parametric_sweep.py` | `sweep_energies`, `sweep_layups`, `sweep_thicknesses` with `on_error` (`raise`/`skip`/`warn`) and optional `progress_callback`; CSV output, partial-result preserving. |
| `cli.py` | — | `bvidfe` entry point; argparse with `choices=` for `--material` / `--loading`, comma-separated `--tier` list, custom positive-float / existing-path types, `--quick` (bare scalar) and `--quick-json` (single-line object) output modes, `--list-materials`. |

Outside the package:

| Path | Role |
|------|------|
| `app.py` | Streamlit web app (`streamlit run app.py`): sidebar inputs, tabs for Summary, Damage Map, Knockdown Curve, 3D Damage, Buckling, Damage Severity, Sweep and Export NCR. |
| `validation/validate_bvid_public.py` | Validation harness: runs each case in `validation/datasets/*.json` through one tier and reports MAE against measured strength; `--gate` fails CI. See `validation/README.md`. |
| `scripts/generate_schemas.py` | Writes `docs/schemas/*.json` from the dataclasses; `--check` fails if they are stale (run in CI). |
| `scripts/self_test.py`, `scripts/benchmark.py`, `scripts/validate_inputs.py` | Install self-test, per-tier timing, and an input-sensitivity sweep. |

## Public API

```python
# Most-used types are re-exported at the top level
from bvidfe import (
    AnalysisConfig, AnalysisResults, BvidAnalysis, MeshParams,
    ImpactEvent, ImpactorGeometry, PanelGeometry,
    MATERIAL_LIBRARY, OrthotropicMaterial,
    DamageState, DelaminationEllipse, load_cscan_json,
)

# The deep import paths work too
from bvidfe.analysis import AnalysisConfig, BvidAnalysis, AnalysisResults
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.damage.io import load_cscan_json, save_cscan_json

# Lower-level access (advanced use)
from bvidfe.impact.olsson import threshold_load, onset_energy
from bvidfe.failure.laminate_strength import unnotched_strength
from bvidfe.failure.soutis_openhole import soutis_cai, whitney_nuismer_tai
from bvidfe.failure.evaluator import FailureEvaluator
from bvidfe.sweep.parametric_sweep import sweep_energies, sweep_layups, sweep_thicknesses
```

## Data-Flow Diagram

```
  ┌─────────────────────┐     ┌──────────────────────┐
  │   ImpactEvent        │     │   C-scan JSON / dict  │
  │  (energy, impactor,  │     │  (delaminations,      │
  │   mass, location)    │     │   dent, fiber break)  │
  └────────┬────────────┘     └──────────┬────────────┘
           │  impact/mapping.py           │  damage/io.py
           │  Olsson threshold            │  load_cscan_json()
           │  shape_templates DPA         │
           │  dent_model                  │
           └───────────┬──────────────────┘
                       ▼
              ┌─────────────────┐
              │   DamageState   │
              │  delaminations  │
              │  dent_depth_mm  │
              │  dpa_mm2        │
              └────────┬────────┘
                       │  analysis/bvid.py  BvidAnalysis.run()
                       │
                       │  pristine = unnotched_strength()  (failure/laminate_strength.py)
          ┌────────────┼────────────────────┐
          ▼            ▼                    ▼
    empirical    semi_analytical          fe3d
    (Soutis CAI, (min of Soutis and   (hex mesh: min of
    WN TAI; ms)   delamination         delamination growth
                  growth; WN TAI)      and first-ply failure;
                                       damaged / undamaged
                                       run; minutes)
          │            │                    │
          └────────────┴────────────────────┘
                       ▼
              ┌─────────────────────┐
              │   AnalysisResults   │
              │  residual_strength  │
              │  pristine_strength  │
              │  knockdown          │
              │  damage (echoed)    │
              │  buckling_eigenvals │
              │  critical_sublam.   │
              │  notes / warnings   │
              └─────────────────────┘
```

## Roadmap

### Shipped in v0.2.0
- **Streamlit web app** (`app.py`), which replaced the PyQt6 desktop GUI (#57): sidebar inputs, result tabs, parametric sweeps, shareable URL state.
- **Per-ply thickness** throughout CLT, the closed-form tiers and the fe3d mesh.
- **Component-wise fe3d damage scaling** (in-plane vs out-of-plane) replacing the prior uniform `DAMAGE_STIFFNESS_FACTOR`.
- **Through-thickness strengths and 1-3 shear strength** (`Zt`, `Zc`, `S13`) on `OrthotropicMaterial` with transverse-isotropy fallback to `Yt` / `Yc` / `S12`.
- **Hex8 Jacobian validation** raising `DegenerateElementError` on non-positive `det(J)`.
- **Robust CLI input validation**: `choices=` materials, positive-float numeric flags, existing-path C-scan.
- **Sweep `on_error` / `progress_callback`** preserving partial results on per-iteration failure.
- **Validation harness** in CI, against a synthetic self-check dataset only.

### Unreleased (since v0.2.0; see CHANGELOG)
- **Rayleigh-Ritz closed-form buckling** in the `fe3d` tier (`fe3d_cai_buckling`, delegated to `semi_analytical.panel_buckling_load`, #129), replacing the 3D K_g eigensolve. It reports buckling onset only.
- **Pristine strength** is the unnotched laminate strength from ply-discount last-ply failure (`failure/laminate_strength.py`), not a thickness-weighted ply-strength average.
- **Delamination-growth CAI** in `semi_analytical` and `fe3d`: compression strength is the stress at which a buckled sublaminate grows its delamination, not the stress at which it buckles.
- **fe3d sees impact damage**: soft-inclusion damage zone, strain-equivalent failure stresses, reaction-force first-ply-failure stress, and knockdown = damaged / undamaged fe3d run.
- **AS4/8552 preset** and the first **measured validation datasets** (`ncamp_as4_8552_cai`, `lovejoy_scotti_im7_8552_cai`); `soutis_k_s` calibrated on them (3.5, then 4.7).
- **T800S/3900-2B preset** and the Hasebe et al. CAI dataset (reported, not gated); impactor shape and diameter factors fitted to it.
- **Energy-balance damage area** (`DPA = olsson_alpha (E − E_onset) / G_IIc`, default `olsson_alpha` 0.078 fitted to Lovejoy & Scotti's C-scans) and a **fixed reference area** in the Soutis knockdown; `soutis_k_s` 4.7.
- **Delaminations confined to the impacted panel** (`distribute_damage(max_half_extent_mm=...)`); the NCAMP and Lovejoy & Scotti cases are impacted on the ASTM D7136 125 × 75 mm window.
- **fe3d damage-zone in-plane factor** 0.30 -> 0.20 -> 0.17, fitted to the first NCAMP, Lovejoy & Scotti quasi-isotropic and Lovejoy & Scotti hard CAI coupons.
- CI checks the generated JSON schemas (`scripts/generate_schemas.py --check`) and gates the semi-analytical tier on the measured datasets.

### Planned
- **True cohesive surfaces**: replace the stiffness-reduction approximation with a bilinear traction-separation law in the `fe3d` tier, which would also give fe3d first-ply failure a characteristic length.
- **3D buckling FE path** with proper load-introduction BCs and plate/shell elements (DKQ or MITC4), so buckling onset can capture damage-zone heterogeneity beyond the closed form (#129).
- **More validated datasets**: the candidates and what each needs are listed in `validation/README.md` ("Next datasets").
- **Calibrated material constants**: `dent_beta` and related parameters against material test data (`soutis_k_s`, the default `olsson_alpha`, the impactor shape and diameter factors and the fe3d `DAMAGE_ZONE_INPLANE_FACTOR` are done; the last fits a quasi-isotropic and a hard laminate within 5%, on three coupons), with held-out data to check them.
- **fe3d on measured data in CI**: CI runs fe3d on three synthetic cases as an ungated smoke run, because a 24-ply NCAMP coupon takes minutes.
