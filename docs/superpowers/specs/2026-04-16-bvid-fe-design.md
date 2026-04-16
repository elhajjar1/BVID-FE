# BVID-FE Design Specification

**Date:** 2026-04-16
**Status:** Draft (pending spec-document-review)
**Authors:** Rani Elhajjar (project owner) + Claude Code (brainstorming partner)

## 1. Purpose and Scope

BVID-FE is a Python package and desktop application for predicting residual strength and stiffness of fiber-reinforced composite laminates containing Barely Visible Impact Damage (BVID). The tool targets composites engineers working on design allowables, certification substantiation, and in-service damage disposition.

BVID-FE is the third in a family of defect-specific composite tools, joining **PorosityFE** (porosity defects) and **WrinkleFE** (fiber waviness). The three tools share material models, laminate theory, failure criteria, deployment pattern (PyPI + PyInstaller desktop app), and validation conventions.

### 1.1 Out of Scope

- Dynamic impact simulation with explicit time integration (use a commercial dynamics code if that fidelity is required).
- Fatigue life prediction (post-impact S-N). Scope for a future v2.
- Environmental conditioning (moisture, thermal cycling).
- Stiffened-panel-level buckling with BVID interacting with bay boundaries. v1 is flat panels.
- In-plane shear-after-impact. Small validation base and low demand.

### 1.2 Primary use cases (explicit)

1. **Design-allowable screening.** Engineer specifies layup, panel size, impact threat (energy, impactor), gets CAI / TAI prediction and a knockdown vs energy curve.
2. **In-service disposition.** Engineer imports a C-scan / dent measurement, gets residual strength estimate for a specific damaged panel.
3. **Parametric study.** Engineer sweeps energy, layup, or thickness and exports CSV results.

## 2. Architecture

### 2.1 Module dependency diagram

```
core ──► impact ──► damage ──► elements ──► solver ──► failure ──► analysis ──► gui
 │                                                                      ▲
 └────────────────────────────────── viz ──────────────────────────────┤
                                                           sweep ──────┘
```

Each layer depends only on layers to its left. No back-edges. `core` has no internal dependencies.

### 2.2 Module catalog

| Module | Role |
|---|---|
| `core/material.py` | `OrthotropicMaterial` dataclass + `MaterialLibrary` with presets AS4/3501-6, IM7/8552, T700/2510, T800/epoxy. Ported from WrinkleFE unchanged. |
| `core/laminate.py` | `Laminate`, `Ply`, `LoadState`; CLT ABD matrices. Ported from WrinkleFE. |
| `core/geometry.py` | `PanelGeometry`, `ImpactorGeometry`, `BoundaryKind` (clamped / simply-supported / free). |
| `impact/olsson.py` | Olsson quasi-static damage-threshold load + onset energy. |
| `impact/shape_templates.py` | Layup-dependent per-interface "peanut" templates distributing DPA into elliptical delaminations. |
| `impact/dent_model.py` | Thickness-normalized empirical dent-depth fits. |
| `impact/mapping.py` | Orchestrates the forward direction: `ImpactEvent → DamageState`. |
| `damage/state.py` | `DamageState` dataclass + `DelaminationEllipse`. |
| `damage/io.py` | NDE import: C-scan JSON/CSV schema, manual entry helpers, validation. |
| `elements/hex8.py`, `hex8i.py`, `gauss.py` | Ported from WrinkleFE. |
| `elements/cohesive.py` | 8-node zero-thickness cohesive surface element (bilinear traction-separation). |
| `solver/static.py`, `assembler.py`, `boundary.py` | Linear static FE assembly/solve/BC handling. Ported from WrinkleFE. |
| `solver/buckling.py` | Linear eigenvalue buckling (SciPy sparse `eigs`) for sublaminate-buckling CAI. |
| `failure/larc05.py`, `tsai_wu.py` | Ported from WrinkleFE. |
| `failure/soutis_openhole.py` | Soutis-style open-hole-equivalent for post-buckling CAI and TAI. |
| `failure/evaluator.py` | Applies criteria across elements; produces `LaminateFailureReport`. |
| `analysis/config.py` | `AnalysisConfig` dataclass. |
| `analysis/bvid.py` | `BvidAnalysis(AnalysisConfig).run()` — the main orchestrator. |
| `viz/plots_2d.py` | Damage-map ellipse overlay, knockdown curves, per-tier comparison charts (matplotlib). |
| `viz/plots_3d.py` | PyVista 3D mesh with delamination surface highlighting, buckling mode shape, stress contour. |
| `viz/style.py` | Publication styling constants. |
| `sweep/parametric_sweep.py` | Sweep driver over energy / layup / thickness; CSV output. |
| `gui/` | PyQt6 desktop app. |

### 2.3 Public API (headline)

```python
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry
from bvidfe.impact.mapping import ImpactEvent

cfg = AnalysisConfig(
    material="IM7/8552",
    layup_deg=[0, 45, -45, 90] * 4,
    ply_thickness_mm=0.152,
    panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
    impact=ImpactEvent(
        energy_J=30,
        impactor=ImpactorGeometry(diameter_mm=16),
        mass_kg=5.5,
    ),
    loading="compression",
    tier="semi_analytical",
)
result = BvidAnalysis(cfg).run()
print(result.summary())   # residual CAI strength, DPA, dent, per-interface delams
```

Lower-level modules are importable for advanced use (e.g., `from bvidfe.impact.olsson import threshold_load`).

## 3. Data Model

### 3.1 Dataclasses

```python
@dataclass
class OrthotropicMaterial:          # ported from WrinkleFE
    name: str
    E11: float; E22: float; nu12: float
    G12: float; G13: float; G23: float
    Xt: float; Xc: float; Yt: float; Yc: float
    S12: float; S23: float
    G_Ic: float; G_IIc: float       # interlaminar fracture toughness (N/mm)
    rho: float

@dataclass
class PanelGeometry:
    Lx_mm: float
    Ly_mm: float
    boundary: Literal["clamped", "simply_supported", "free"] = "simply_supported"

@dataclass
class ImpactorGeometry:
    diameter_mm: float = 16.0
    shape: Literal["hemispherical", "flat", "conical"] = "hemispherical"

@dataclass
class ImpactEvent:
    energy_J: float
    impactor: ImpactorGeometry
    mass_kg: float
    location_xy_mm: tuple[float, float] = (0.0, 0.0)   # default: panel center

@dataclass
class DelaminationEllipse:
    interface_index: int                               # 0 = between plies 0 and 1
    centroid_mm: tuple[float, float]
    major_mm: float
    minor_mm: float
    orientation_deg: float                             # major axis in panel frame
    @property
    def area_mm2(self) -> float:                       # π·a·b
        return math.pi * self.major_mm * self.minor_mm

@dataclass
class DamageState:
    delaminations: list[DelaminationEllipse]
    dent_depth_mm: float
    fiber_break_radius_mm: float = 0.0
    @property
    def projected_damage_area_mm2(self) -> float: ...  # union of ellipse footprints
    @property
    def per_interface_area(self) -> dict[int, float]: ...

@dataclass
class MeshParams:
    elements_per_ply: int = 4
    in_plane_size_mm: float = 1.0
    cohesive_zone_factor: float = 1.0                  # tunes bilinear law

@dataclass
class AnalysisConfig:
    material: str | OrthotropicMaterial
    layup_deg: list[float]
    ply_thickness_mm: float
    panel: PanelGeometry
    loading: Literal["compression", "tension"]         # CAI or TAI
    tier: Literal["empirical", "semi_analytical", "fe3d"]
    impact: ImpactEvent | None = None                  # exactly ONE of these two
    damage: DamageState | None = None
    mesh: MeshParams | None = None                     # required only for fe3d

    def __post_init__(self):
        assert (self.impact is None) ^ (self.damage is None), \
            "Provide exactly one of impact or damage"
        if self.tier == "fe3d" and self.mesh is None:
            self.mesh = MeshParams()

@dataclass
class AnalysisResults:
    residual_strength_MPa: float
    pristine_strength_MPa: float
    knockdown: float                                   # residual / pristine
    damage: DamageState                                # always populated
    dpa_mm2: float
    buckling_eigenvalues: list[float] | None = None    # semi_analytical + fe3d
    critical_sublaminate: int | None = None            # interface index
    field_results: FieldResults | None = None          # fe3d only
    tier_used: str
    config_snapshot: dict
    def summary(self) -> str: ...
    def to_dict(self) -> dict: ...                     # for JSON export
```

### 3.2 Design notes on the data model

- **`DamageState` is the universal handoff** between the impact-mapping stage and the residual-strength engine. Both workflow paths converge on a `DamageState` before hitting the solver tier.
- **`AnalysisResults.damage` is always populated.** If the user provided a `DamageState` as input, it is echoed back. If they provided an `ImpactEvent`, the impact mapping writes its output there. Downstream tooling (comparison studies, visualization) can then treat the two paths identically.
- **`DelaminationEllipse.interface_index`** lets the C-scan importer and the 3D FE mesh builder index plies consistently (interface `i` = between ply `i` and `i+1`).
- **Tier-specific outputs** (buckling eigenvalues, field results) are optional — empty for tiers that don't compute them.
- **Immutability:** all dataclasses are `frozen=False` to allow fluent GUI editing, but `config_snapshot` in `AnalysisResults` is a deep-copied dict to preserve provenance.

### 3.3 C-scan JSON schema (summary)

Documented in `docs/cscan_schema.md`:

```json
{
  "schema_version": "1.0",
  "panel": {"Lx_mm": 150, "Ly_mm": 100},
  "layup_deg": [0, 45, -45, 90, 90, -45, 45, 0],
  "dent_depth_mm": 0.45,
  "fiber_break_radius_mm": 3.0,
  "delaminations": [
    {"interface_index": 3,
     "centroid_mm": [75, 50],
     "major_mm": 28, "minor_mm": 18, "orientation_deg": 45}
  ]
}
```

The importer validates: non-negative dimensions, interface indices within `[0, n_plies-1]`, ellipse footprints within panel bounds.

## 4. Physics Models

### 4.1 Impact mapping (`ImpactEvent → DamageState`)

Executed only when `AnalysisConfig.impact` is provided. Three stages.

**4.1.1 Olsson quasi-static damage-threshold load** (`impact/olsson.py`).

Closed-form threshold load from plate bending + delamination fracture energy balance:

```
P_c = π · √(8 · G_IIc · D_eff / 9)
E_onset = P_c² / (2 · k_cb)
k_cb = 1 / (1/k_bending + 1/k_contact_Hertz)
```

where `D_eff` is the effective flexural rigidity from CLT ABD (geometric mean of D11, D22), `k_bending` is a Navier series for the panel bending stiffness at the impact location, and `k_contact_Hertz` is the Hertzian contact stiffness of the impactor on the laminate surface. If `E_impact ≤ E_onset`, no damage is produced (`DamageState` is empty).

**4.1.2 Total damage projected area**:

```
DPA = α · (E_impact - E_onset) / (G_IIc · h)
```

with `h` the laminate thickness and `α` a material-family calibration constant (~0.8 for CFRP, from Olsson 2001 and subsequent fits). Calibration constants live in `core/material.py` alongside each preset.

**4.1.3 Per-interface distribution + dent + fiber break** (`impact/shape_templates.py`, `dent_model.py`).

- DPA is distributed across interfaces using a layup-dependent "peanut" template: ellipses are largest near the back face and oriented along the back-face ply direction. Ellipse aspect ratio per interface `AR_i = f(Δθ_i)` where `Δθ_i` is the ply-angle mismatch across interface `i`. Larger mismatch → more elongated lobe. The template conserves total DPA after accounting for overlap (polygon union).
- Dent depth: thickness-normalized empirical fit `d/h = β · ((E - E_onset) / (G_Ic · h²))^γ`, with `β` and `γ` material-specific constants.
- Fiber-break core radius: `r_fb = η · √(max(0, E - E_fb_threshold))`. Often `η = 0` or `E_fb_threshold` is very high, yielding no fiber-break core (which matches most BVID tests).

### 4.2 Residual-strength engine — three tiers × two loading modes

The tier is selected by `AnalysisConfig.tier`.

**4.2.1 Empirical tier** (`tier="empirical"`).

- **CAI** — Soutis-style knockdown curve:
  ```
  σ_CAI / σ_0 = 1 / (1 + k_s · (DPA / A_panel)^m)
  ```
  `k_s`, `m` material-specific. Calibrated against CAI test data in the validation datasets.
- **TAI** — Whitney-Nuismer point-stress on the equivalent open hole:
  ```
  σ_TAI / σ_0 = f_WN(DPA, d_0)
  ```
  with characteristic distance `d_0` as a per-material constant.

**4.2.2 Semi-analytical tier** (`tier="semi_analytical"`).

- **CAI** — sublaminate buckling + post-buckling propagation:
  1. Identify the critical delaminated interface: largest through-thickness asymmetry × largest ellipse.
  2. Compute sublaminate buckling load via orthotropic Rayleigh-Ritz over the elliptical delamination footprint (5-term series is adequate for engineering accuracy). Simply-supported boundary on the ellipse perimeter.
  3. Post-buckling propagation via ERR-based advance (Griffith-style) coupled with Soutis's residual-strength envelope. Returns `σ_CAI`.
- **TAI** — Soutis notch model: treat DPA as an equivalent open hole of diameter `2·√(DPA/π)`, apply cohesive-zone / average-stress criterion with in-situ tensile strength.

**4.2.3 3D FE tier** (`tier="fe3d"`).

- Structured hex mesh: resolution from `MeshParams` (default ~4 elements per ply through thickness, ~1 mm in-plane). One element set per ply with orthotropic stiffness rotated by ply angle.
- **Zero-thickness cohesive surface elements** inserted at the delaminated interfaces only:
  - Inside each `DelaminationEllipse` footprint: tractions held at zero (pre-cracked).
  - Outside ellipses: bilinear traction-separation law with strengths derived from material `G_Ic`, `G_IIc`.
- **Fiber-break core:** elements inside `fiber_break_radius_mm` get isotropic near-zero stiffness (~1 MPa), same treatment PorosityFE uses for void elements.
- **CAI path:** assemble linear stiffness K from the damaged mesh, solve generalized eigenproblem `K · φ = λ · K_g · φ` for the lowest mode via `scipy.sparse.linalg.eigs`, scale mode to unit amplitude, evaluate LaRC05 at post-buckled stress state. Residual strength = applied load × eigenvalue at first-ply failure.
- **TAI path:** static solve under prescribed uniaxial tension displacement; Tsai-Wu evaluated at every Gauss point; fail when max Tsai-Wu index = 1. No buckling step.

All tiers return the same `AnalysisResults` schema. `field_results` is only populated for `fe3d`.

## 5. Workflow / End-to-End Data Flow

```
                      ┌────────────────────────────────┐
                      │        AnalysisConfig          │
                      │  (material, layup, panel, ...) │
                      └──────────────┬─────────────────┘
                                     │
                       has impact?   │   has damage?
                  ┌──────────────────┴──────────────────┐
                  ▼                                     ▼
      ┌────────────────────┐                 ┌──────────────────────┐
      │  impact.mapping    │                 │  damage.io (optional)│
      │   Olsson + shape   │                 │   C-scan / manual    │
      │   templates + dent │                 │   → DamageState      │
      └──────────┬─────────┘                 └──────────┬───────────┘
                 │  DamageState                         │  DamageState
                 └──────────────────┬───────────────────┘
                                    ▼
                        ┌───────────────────────┐
                        │   tier dispatcher     │
                        └──────────┬────────────┘
                                   │
                  ┌────── empirical ──┼── semi_analytical ──┬── fe3d ──┐
                  ▼                   ▼                     ▼          ▼
              Soutis CAI        Sublaminate             Linear      Pre-cracked
              + WN TAI          Rayleigh-Ritz           buckling    cohesive +
                                buckling + Soutis       eigensolve  static/eig
                                                                    + LaRC05
                  └───────────────────┬─────────────────────────────────┘
                                      ▼
                        ┌────────────────────────────┐
                        │  σ_residual, knockdown,    │
                        │  buckling modes, fields    │
                        └──────────────┬─────────────┘
                                       ▼
                               ┌─────────────────┐
                               │ AnalysisResults │
                               └─────────────────┘
```

Orchestration in `analysis/bvid.py`:

```python
class BvidAnalysis:
    def __init__(self, config: AnalysisConfig):
        self.config = config
    def run(self) -> AnalysisResults:
        damage = self._resolve_damage_state()     # impact→damage OR echo input
        solver = self._select_solver()            # empirical/semi/fe3d
        residual = solver.solve(damage, self.config.loading)
        return AnalysisResults(
            residual_strength_MPa=residual.sigma,
            pristine_strength_MPa=residual.sigma_pristine,
            knockdown=residual.sigma / residual.sigma_pristine,
            damage=damage,
            dpa_mm2=damage.projected_damage_area_mm2,
            buckling_eigenvalues=residual.buckling_eigs,
            critical_sublaminate=residual.critical_interface,
            field_results=residual.field_results,
            tier_used=self.config.tier,
            config_snapshot=deepcopy(asdict(self.config)),
        )
```

**Key invariant:** The residual-strength engine is fully testable in isolation using synthetic `DamageState` objects. This means the impact-mapping code and the residual-strength code can be developed, tested, and validated independently.

## 6. GUI

PyQt6 desktop application. Main window: left dock (input panels), tabbed central area (results), menu bar (file I/O).

### 6.1 Input panels (left dock)

| Panel | Contents |
|---|---|
| `MaterialPanel` | Material library dropdown + editable orthotropic fields; layup (comma-separated angles) + ply thickness. |
| `PanelPanel` | Panel `Lx`, `Ly`; boundary condition dropdown. |
| `InputModePanel` | Radio: *Impact-driven* vs *Damage-driven*. Switches active panel below. |
| `ImpactPanel` | Impact energy (J), impactor diameter, impactor mass, location. Shows computed `E_onset` live. |
| `DamagePanel` | Damage-entry table (one row per delaminated interface: index, major, minor, orientation, centroid); dent-depth field; fiber-break-radius field; **Import C-scan…** button. |
| `AnalysisPanel` | Tier dropdown (empirical / semi_analytical / fe3d); loading dropdown (compression / tension); mesh params (only when fe3d selected); **Run** button. |
| `SweepPanel` | Parametric sweep tab (energy range, layup variations, output CSV path). |

### 6.2 Results tabs (central area)

| Tab | Contents |
|---|---|
| Summary | Residual strength, pristine, knockdown, DPA, dent, tier, run time. |
| Damage Map | Top-down ellipse overlay (all interfaces, color-coded by depth) + side-view dent profile. |
| Knockdown Curve | Residual strength or KD vs impact energy (or vs DPA) — re-runs the selected tier at several energies. |
| Mesh (fe3d only) | PyVista 3D mesh with delamination surface highlights. |
| Buckling Mode (fe3d / semi_analytical) | Eigenvector contour. |
| Stress Field (fe3d only) | σ_xx / σ_zz contour. |

### 6.3 Threading

`AnalysisWorker(QThread)` and `SweepWorker(QThread)` run compute off the UI thread. Workers emit `resultReady(AnalysisResults)` / `progress(int)` / `error(str)` signals. Pattern ported from WrinkleFE.

### 6.4 File I/O

- Save/Load `AnalysisConfig` as JSON via File menu.
- Export results: PNG per tab, JSON of full `AnalysisResults.to_dict()`, CSV of sweep outputs.
- C-scan import: JSON per `docs/cscan_schema.md`.

## 7. Testing Strategy

Target: ~200-400 unit tests (PorosityFE has 165+, WrinkleFE has 414). Organized in `tests/` mirroring package layout.

| Test file | Coverage |
|---|---|
| `test_core_material.py` | Library loading, dataclass defaults, CLT ABD regression. |
| `test_core_laminate.py` | CLT golden values (ported from WrinkleFE). |
| `test_core_geometry.py` | Panel / impactor / boundary types. |
| `test_impact_olsson.py` | Threshold formula, onset energy, edge cases (G_IIc → 0, E < E_onset). |
| `test_impact_shape_templates.py` | Ellipse AR vs Δθ; DPA conservation under union; per-interface count equals `n_plies - 1`. |
| `test_impact_dent_model.py` | Monotonic in energy; zero below threshold; thickness scaling. |
| `test_impact_mapping.py` | End-to-end: impact → damage state with expected DPA / dent / delam count. |
| `test_damage_state.py` | `projected_damage_area_mm2` for overlapping ellipses; serialization round-trip. |
| `test_damage_io.py` | C-scan JSON round-trip; schema validation errors. |
| `test_elements_cohesive.py` | Zero-thickness element stiffness, bilinear traction-separation, ERR on analytical DCB. |
| `test_solver_buckling.py` | Orthotropic plate buckling benchmarks (compare to classical closed-form). |
| `test_failure_soutis.py` | Soutis open-hole golden values. |
| `test_analysis_pipeline.py` | End-to-end both paths; impact-driven and damage-driven yield same `AnalysisResults` when inputs are consistent. |
| `test_gui_smoke.py` | pytest-qt headless: panels construct, signals fire, worker runs mock. |

All tests run under `pytest -v`. CI gates on green tests.

## 8. Validation

Mirrors WrinkleFE / PorosityFE validation structure.

```
validation/
├── validate_bvid_public.py             # regenerates all validation figures
├── reference_data.json                 # all digitized test data, single source of truth
├── datasets/
│   ├── soutis_as4_3501-6.json
│   ├── caprino_as4_epoxy.json
│   ├── sanchez-saez_im7_8552.json
│   └── nasa_cai_round_robin.json
└── figures/
    ├── validation_cai_all.png
    └── validation_per_dataset_*.png
```

**Target pass table (published in README):**

| Dataset | Loading | Cases | Target MAE | Reference |
|---|---|---|---|---|
| Soutis AS4/3501-6 | CAI | ~15 | <12% | Soutis & Curtis (1996) |
| Caprino AS4/epoxy | TAI + CAI | ~10 | <15% | Caprino (1984) |
| Sanchez-Saez IM7/8552 | CAI | ~12 | <12% | Sanchez-Saez et al. (2005) |
| NASA round-robin | CAI | ~10 | <15% | NASA/TM-2007 |

`python validation/validate_bvid_public.py` is run as a CI job to prevent regressions.

## 9. Packaging and Distribution

### 9.1 Repository layout

```
BVID-FE/
├── bvidfe/                  # package (per Section 2)
├── tests/                   # pytest suite
├── validation/              # datasets + script + figures
├── docs/
│   ├── superpowers/
│   │   ├── specs/           # this document
│   │   └── plans/           # implementation plan (next step)
│   └── cscan_schema.md
├── figures/                 # README screenshots
├── joss/                    # JOSS paper (deferred to v1.1)
├── screenshots/             # GUI screenshots for README
├── pyproject.toml
├── BvidFE.spec              # PyInstaller spec (macOS + Windows)
├── requirements.txt
├── CLAUDE.md
├── README.md
├── ARCHITECTURE.md
├── CITATION.cff
├── CONTRIBUTING.md
├── CHANGELOG.md
└── LICENSE                  # MIT
```

### 9.2 Dependencies

Runtime (conservative, matches siblings):
`numpy`, `scipy`, `matplotlib`, `PyQt6`, `pyvista`, `pyvistaqt`.
Dev: `pytest`, `pytest-qt`, `pyinstaller`.

### 9.3 CLI entry points (in `pyproject.toml`)

- `bvidfe-gui` → `bvidfe.gui.app:main`
- `bvidfe` → `bvidfe.cli:main` (batch/scripted runs, sweeps, validation driver)

### 9.4 Distribution

- Source: `pip install -e ".[all]"` from GitHub.
- macOS + Windows standalone app via `pyinstaller BvidFE.spec`.
- GitHub Actions CI: pytest matrix over macOS / Ubuntu / Windows, Python 3.9-3.12.
- Tagged releases produce GitHub Release artifacts (match WrinkleFE v1.1.0 convention).

### 9.5 Documentation

- **README.md:** screenshots, validation table, quick-start (GUI + Python API), citation block.
- **ARCHITECTURE.md:** module dependency diagram, data flow.
- **CLAUDE.md:** contributor / future-Claude guidance.
- **docs/cscan_schema.md:** C-scan input format.
- **JOSS paper:** deferred to v1.1.

### 9.6 License and citation

MIT License, matching siblings. CITATION.cff and BibTeX stub in README, to be completed when a DOI is minted.

## 10. Open Items (to be resolved during implementation)

- Exact calibration values `α, β, γ, η, k_s, m, d_0` per material preset — populated in the implementation plan after the first validation dataset is ingested.
- Whether to include an optional "quick-sanity" tier below empirical (dent-depth-only lookup); left out of v1 unless needed.
- Whether `fe3d` CAI should use geometric nonlinearity (arc-length) for post-buckling — v1 uses linear eigenvalue + Koiter-style first-ply-failure heuristic; nonlinear arc-length deferred to v1.1 if validation MAE demands it.

## 11. References

Physics and calibration sources that the models draw on.

- Olsson, R. (2001). "Analytical prediction of large mass impact damage in composite laminates." *Composites Part A*, 32(9).
- Olsson, R. (2010). "Analytical model for delamination growth during small mass impact on plates." *International Journal of Solids and Structures*, 47(21).
- Soutis, C. & Curtis, P.T. (1996). "Prediction of the post-impact compressive strength of CFRP laminated composites." *Composites Science and Technology*, 56(6).
- Soutis, C., Curtis, P.T. & Fleck, N.A. (1993). "Compressive failure of notched carbon fibre composites." *Proc. R. Soc. Lond. A*, 440.
- Caprino, G. (1984). "Residual strength prediction of impacted CFRP laminates." *Journal of Composite Materials*, 18(6).
- Sanchez-Saez, S., Barbero, E., Zaera, R. & Navarro, C. (2005). "Compression after impact of thin composite laminates." *Composites Science and Technology*, 65(13).
- Whitney, J.M. & Nuismer, R.J. (1974). "Stress fracture criteria for laminated composites containing stress concentrations." *Journal of Composite Materials*, 8(3).
- Pinho, S.T., Davila, C.G., Camanho, P.P., Iannucci, L., Robinson, P. (2005). "Failure models and criteria for FRP under in-plane or three-dimensional stress states including shear non-linearity." NASA/TM-2005-213530. (LaRC05 criterion.)
- Tsai, S.W. & Wu, E.M. (1971). "A general theory of strength for anisotropic materials." *Journal of Composite Materials*, 5(1).
- Elhajjar, R. (2025). "Fat-tailed failure strength distributions and manufacturing defects in advanced composites." *Scientific Reports*, 15:25977.
