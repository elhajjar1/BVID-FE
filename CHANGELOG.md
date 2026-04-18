# Changelog

All notable changes to BVID-FE are documented in this file.

## [0.2.0-dev] - unreleased

In-progress work toward v0.2.0. No tag yet.

### Added

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
