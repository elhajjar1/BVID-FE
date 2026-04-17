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
- **Working 3D Mesh GUI tab** via `pyvistaqt.QtInteractor`. The old
  `v0.2.0` placeholder is replaced with a real interactive 3D viewer that
  renders the FE mesh coloured by per-element damage factor. Works for
  all tiers (empirical / semi-analytical / fe3d) — for tiers that don't
  otherwise build a mesh, the viewer constructs one with default
  `MeshParams` just for visualization.
- **Validation harness** (`validation/validate_bvid_public.py`):
  `DatasetCase` dataclass + `run_dataset()` + MAE% metric. Auto-discovers
  any JSON dataset in `validation/datasets/`. Ships with a synthetic
  self-check dataset (MAE ≈ 0% by construction) so the harness is
  exercised in CI; real published datasets (Soutis, Caprino,
  Sanchez-Saez, NASA round-robin) remain to be digitized by hand.
- **CI regression gate** — new `validation` job in
  `.github/workflows/tests.yml` running
  `python validation/validate_bvid_public.py --gate`.
- Seven new tests on the geometric-stiffness / buckling path; three on
  the 3D Mesh tab; three on the validator smoke path. **192 tests now
  passing** (up from 179 at v0.1.0).

### Known limitations (still deferred)

- Material calibration constants remain placeholders until real datasets
  land.
- `fe3d` tier still uses stiffness-reduction at delaminations rather than
  zero-thickness cohesive surfaces (different physics — stiffness
  reduction captures the stiffness loss but not the debonding/sliding).
- Buckling Mode / Stress Field GUI tabs are still placeholders.
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
