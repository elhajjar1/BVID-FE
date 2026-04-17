# BVID-FE Architecture

## Module Dependency Diagram

```
core ──► impact ──► damage ──► elements ──► solver ──► failure ──► analysis ──► (viz, sweep)
 │                    ▲                                                ▲
 │                    └──────────────────────────────────────────────┘
 └──────────────────────────────── viz ─────────────────────────────┘
```

`core` has no internal dependencies. Each layer depends only on layers to its left.
`damage/state` is consumed by both the `impact` mapping stage and the `analysis` orchestrator.

## Module Catalog

| Module | Submodule | Role |
|--------|-----------|------|
| `core/` | `material.py` | `OrthotropicMaterial` dataclass + `MaterialLibrary` with four built-in presets (AS4/3501-6, IM7/8552, T700/2510, T800/epoxy). Ported from WrinkleFE. |
| | `laminate.py` | `Laminate`, `Ply`, `LoadState`; Classical Lamination Theory ABD matrices and effective engineering constants. Ported from WrinkleFE. |
| | `geometry.py` | `PanelGeometry`, `ImpactorGeometry`, `BoundaryKind` (clamped / simply-supported / free). |
| `impact/` | `olsson.py` | Olsson quasi-static damage-threshold load `P_c` and onset energy `E_onset` from plate bending + fracture energy balance. |
| | `shape_templates.py` | Layup-dependent per-interface "peanut" templates that distribute the total DPA into per-interface elliptical delaminations. |
| | `dent_model.py` | Thickness-normalized empirical dent-depth model producing `dent_depth_mm` from impact energy. |
| | `mapping.py` | Orchestrates the impact-driven workflow: `ImpactEvent` → `DamageState`. |
| `damage/` | `state.py` | `DamageState` dataclass + `DelaminationEllipse`; shapely-union projected damage area computation. |
| | `io.py` | C-scan JSON/CSV import, validation, and manual-entry helpers per `docs/cscan_schema.md`. |
| `elements/` | `hex8.py` | Standard 8-node isoparametric hexahedral element: shape functions, B-matrix, 24x24 stiffness. |
| | `hex8i.py` | Incompatible-modes hex element (Wilson modes) for improved bending accuracy. |
| | `gauss.py` | Gauss-Legendre quadrature points and weights: `gauss_points_1d`, `gauss_points_hex`. |
| | `cohesive.py` | 8-node zero-thickness cohesive surface element with bilinear traction-separation law (stiffness-reduction approximation in v0.1.0). |
| `solver/` | `static.py` | `StaticSolver` — sparse-direct linear static solve (SciPy). |
| | `assembler.py` | COO → CSC global stiffness matrix assembly from element contributions. |
| | `boundary.py` | `BoundaryCondition`, `BoundaryHandler` — penalty-method Dirichlet BCs for compression and tension loading. |
| | `buckling.py` | Linear eigenvalue buckling solve (`scipy.sparse.linalg.eigsh`) for sublaminate-buckling CAI. |
| `failure/` | `larc05.py` | Minimal LaRC05 composite failure criterion (Hashin-3D reduction). |
| | `tsai_wu.py` | Full 3D Tsai-Wu failure criterion with interaction terms. |
| | `soutis_openhole.py` | Soutis open-hole-equivalent CAI model + Whitney-Nuismer TAI (point-stress and average-stress). |
| | `evaluator.py` | `FailureEvaluator` — applies criteria across all elements; produces `LaminateFailureReport`. |
| `analysis/` | `config.py` | `AnalysisConfig` dataclass (material, layup, panel, loading, tier, impact or damage). |
| | `bvid.py` | `BvidAnalysis(AnalysisConfig).run()` — main orchestrator; dispatches to tier. |
| | `results.py` | `AnalysisResults` dataclass with `summary()` and `to_dict()`. |
| | `semi_analytical.py` | Semi-analytical tier implementation (Rayleigh-Ritz sublaminate buckling + Soutis). |
| | `fe_tier.py` | 3D FE tier implementation (mesh build → assemble → solve → failure). |
| | `fe_mesh.py` | Damaged hexahedral mesh construction from `DamageState`. |
| `viz/` | `plots_2d.py` | Damage-map ellipse overlay, knockdown curves, per-tier comparison charts (matplotlib). |
| | `plots_3d.py` | 3D PyVista plots: delamination surface, buckling mode shape, stress contour. |
| | `style.py` | Publication styling constants (fonts, DPI, color maps). |
| `sweep/` | `parametric_sweep.py` | `sweep_energies`, `sweep_layups`, `sweep_thicknesses` — CSV output. |
| `cli.py` | — | `bvidfe` entry point; parses command-line arguments and calls `BvidAnalysis`. |

## Public API

```python
# High-level orchestrator
from bvidfe.analysis import AnalysisConfig, BvidAnalysis, AnalysisResults

# Geometry and event dataclasses
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry
from bvidfe.impact.mapping import ImpactEvent

# Damage state (inspection-driven path)
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.damage.io import load_cscan_json

# Lower-level access (advanced use)
from bvidfe.impact.olsson import threshold_load
from bvidfe.failure.soutis_openhole import soutis_cai_knockdown, whitney_nuismer_tai
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
          ┌────────────┼────────────────────┐
          ▼            ▼                    ▼
    empirical    semi_analytical          fe3d
    (Soutis +    (Rayleigh-Ritz +    (hex mesh +
    WN, ~ms)      Soutis, ~s)         LaRC05/TW,
                                       ~minutes)
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
              │  field_results      │
              └─────────────────────┘
```

## Roadmap: v0.2.0

- **PyQt6 GUI**: material panel, impact/inspection panel, tier selector, results area, sweep tab
- **True cohesive surfaces**: replace stiffness-reduction approximation with full bilinear traction-separation law in the `fe3d` tier
- **Buckling eigenvalue CAI** in the `fe3d` tier (currently deferred; uses first-ply-failure as conservative estimate)
- **Validated datasets**: Soutis (1996), Caprino (1984), Sanchez-Saez (2005), NASA COCOMAT — digitized and integrated into `validation/`
- **Calibrated material constants**: `olsson_alpha`, `soutis_k_s`, `dent_beta`, and related parameters refined against specific material test data
- **PyInstaller macOS app bundle**
