# BVID-FE Improvement Plan

*Prepared 2026-07-16 from a full audit of the codebase, tests, CI, docs, packaging,
validation assets, and the open issue/PR backlog.*

This document is a prioritized, phased plan. Each item carries file/line anchors,
a rationale, and acceptance criteria so it can be turned into issues/PRs directly.
Where an item overlaps an existing GitHub issue, the issue number is cited.

---

## 1. Current-state snapshot

| Dimension | State |
|---|---|
| Version | `pyproject.toml` says `0.2.0` — but see [A4](#a4): four conflicting version strings in-repo |
| Tests | 425 passing in ~18 s, 0 failures |
| Coverage | 93% statement (CI gate at 90%); weakest: `viz/plotly_3d.py` 22%, `solver/buckling.py` 76%, `impact/shape_templates.py` 84% |
| Static typing | `py.typed` shipped, but **no mypy/pyright anywhere** (local or CI) |
| Lint/format | ruff (default rules only) + black, enforced in CI and pre-commit |
| Open issues | 10 (#82, #89, #110, #114, #115, #117, #118, #129, #132, #133) |
| Open PRs | 6 Dependabot bumps + #75 (per-ply thickness, stale) |
| App | `app.py` — 1,608-line Streamlit monolith, import-smoke-tested only |
| Validation | CI gate runs **only** against a synthetic self-check dataset (tautological, ~0% MAE by construction) |

Overall: the library core is in good shape — well-tested, well-guarded numerics in
most hot paths, disciplined warnings/notes plumbing. The dominant problems are
(1) a handful of real correctness/robustness gaps, (2) heavy documentation rot
from the abandoned PyQt GUI era, (3) three hand-rolled copies of config
serialization, and (4) infrastructure that advertises guarantees it doesn't
enforce (typed API with no type checker, a validation "gate" that can't fail
meaningfully, schemas that can drift silently).

---

## 2. Workstream A — Correctness & numerical robustness (P0)

These are defects, not enhancements. All are small; land them first as a `v0.2.1`
patch series.

### A1. fe3d silently computes Tsai-Wu when asked for Puck ⚠ highest priority
`_solve_failure_strain_analytic` (`src/bvidfe/analysis/fe_tier.py:258-267`, scalar
ref `:332-338`) branches `if criterion == "larc05": ... else: <Tsai-Wu>`. The
parameter is typed `CriterionName`, which includes `"puck"` — so passing
`"puck"` silently evaluates Tsai-Wu with no error or warning. Puck (added in
#121) is therefore unreachable from the fe3d tier and dangerous to request.
**Fix:** validate `criterion` at the top of the function; raise `ValueError` for
unsupported names until Puck's analytic-scaling branch exists (its plane-search
form is not quadratic in the load multiplier, so it needs either a dedicated
branch or a documented exclusion).
**Accept:** `pytest` case proving `"puck"` raises (or is correctly evaluated);
type of `criterion` narrowed to the actually-supported subset.

### A2. Whitney–Nuismer denominator can go non-physical
`whitney_nuismer_tai` (`src/bvidfe/failure/soutis_openhole.py:88-90`):
`denom = 2 + ξ² + 3ξ⁴ − (Kt_inf−3)(5ξ⁶−7ξ⁸)`. For strongly orthotropic layups
(production passes the Lekhnitskii Kt from `bvid.py:177` / `semi_analytical.py:423`;
a near-UD laminate gives Kt ≈ 7), the correction term can push `kd = 2/denom`
above 1 (residual > pristine, at Kt ≳ 10) and, in the extreme, through a pole
(denom → 0). Also `lekhnitskii_kt_infinity` (`:33`) takes an unguarded
`math.sqrt` that can receive a negative argument for unusual engineering
constants.
**Fix:** clamp `kd` to `(0, 1]`, raise/warn on non-positive denominator, guard
the sqrt argument.
**Accept:** property test over a grid of layups asserting `0 < kd ≤ 1`.

### A3. `solver/buckling.py` violates its own contract
Docstring promises "smallest **positive** eigenvalues" but both paths return the
raw ascending sort with no sign filtering (`src/bvidfe/solver/buckling.py:51-54`,
`:63-65`) — geometric stiffness is indefinite, so spurious negative eigenvalues
can be returned as buckling factors. The dense fallback (`:68-73`) calls
`scipy.linalg.eigh(K, Kg)` which *requires* PD `Kg` and will raise for exactly
the indefinite matrices the module documents. The sparse `eigsh(sigma=0)`
shift-invert path has no convergence-failure handling. The module is currently
dead in the analysis path (kept deliberately for #132), but #132 builds on it —
fix the contract before the plate/shell work starts.
**Accept:** tests covering positive-filtering, an indefinite-Kg dense case, and
ARPACK non-convergence handling; sparse path gains its first test.

### A4. Version chaos — four different strings
`pyproject.toml:7` = `0.2.0`; `src/bvidfe/__init__.py:8` = `0.2.0.dev0`;
`CITATION.cff:5` = `0.1.0`; `README.md` bibtex = `0.1.0-alpha`. The publish
workflow verifies the tag against `pyproject.toml` only — a `v0.2.0` wheel ships
with `bvidfe.__version__ == "0.2.0.dev0"`, and `bvidfe --version` disagrees with
`pip show`.
**Fix:** single-source via `importlib.metadata.version("bvidfe")` in
`__init__.py` (or a release-time check); extend `publish.yml` tag verification to
`__init__` and `CITATION.cff`; update the `release` skill checklist.
**Accept:** one authoritative version; CI check that all in-repo version strings
agree.

### A5. `assert` used for input validation
`FailureEvaluator._index` validates the stress-field shape with `assert`
(`src/bvidfe/failure/evaluator.py:76`) — stripped under `python -O`, turning a
clear error into downstream garbage. Replace with `ValueError`.

### A6. `STRICT_HEX8I_CONDENSATION` is not actually an environment variable
It's a module-level bool (`src/bvidfe/elements/hex8i.py:34`) settable only by
monkeypatching, despite env-var phrasing in commit messages/docs. Worse, in the
default non-strict mode a singular `Kaa` makes the element **silently degrade to
the locking-prone plain Hex8** (`:104-145`) with only a log line — no
result-level signal, so a mesh can mix element formulations invisibly. Related:
`Kaa` is explicitly inverted rather than solved, and conditioning is only
checked *after* an exception (`:105-109`).
**Fix:** read the flag from the environment like `BVIDFE_FE3D_MAX_DOF` does
(`fe_tier.py:63`); propagate a `hex8i_condensation_fallback` tag into
`AnalysisResults.warnings`; use `np.linalg.solve`/`lstsq` with a pre-check on
conditioning.
**Accept:** env-resolution test (mirroring `tests/analysis/test_env_resolution.py`),
warning-tag test, Hex8i patch test (currently missing — Hex8 has one, Hex8i doesn't).

### A7. Penalty BC conditioning and docstring mismatch
`apply_bcs` adds a fixed `1e10` to the diagonal (`src/bvidfe/solver/boundary.py:55`)
— unscaled relative to `diag(K)`, so conditioning and enforcement accuracy vary
with material/mesh scaling; the module docstring (`:8`) describes a different
(multiplicative) scheme than the code implements. Also `_nodes_on_plane` uses an
absolute `1e-9` mm tolerance (`:60-63`) and `boundary="free"` (`:110-119`) can
leave rotational rigid-body modes unconstrained with no singularity check.
**Fix:** scale penalty as `α · max(diag(K))`, fix docstring, make the plane
tolerance relative to element size, add a rigid-body/conditioning check for
`free`.

### A8. Material-card guards
`get_stiffness_matrix` inverts compliance with no conditioning guard
(`src/bvidfe/core/material.py:147`) — unlike `Laminate`, which guards at
`cond > 1e10`. The `nu12` validation bound (`-1 < ν < 0.5`, `:95`) is an
isotropic rule; the orthotropic stability bound is `ν12 < sqrt(E11/E22)`.
**Fix:** mirror the laminate guard; implement the correct thermodynamic bound.

### A9. Impact-location sentinel collision
`ImpactEvent.location_xy_mm` defaults to `(0.0, 0.0)`, which `mapping.py:126-127`
reinterprets as "panel center" — a genuine corner-origin impact at (0,0) is
unrepresentable. `olsson.onset_energy` uses a different sentinel (`None`) for
the same concept.
**Fix:** default to `None` everywhere; `(0,0)` becomes a real location. Breaking
change — do it while the user base is small, with a CHANGELOG note.

### A10. Remaining small guards
- Puck mode-A envelope factor can silently go negative
  (`src/bvidfe/failure/puck.py:168`) — warn or clamp.
- `knockdown = ... if sigma_0 > 0 else 0.0` (`src/bvidfe/analysis/bvid.py:132`)
  silently reports 0.0 for a degenerate pristine strength — raise instead.
- `save_cscan_json` writes non-atomically with no error handling
  (`src/bvidfe/damage/io.py:149-150`) — temp-file-and-rename.
- `laminate._reduced_stiffness` denominator `1 − ν12·ν21` unguarded
  (`src/bvidfe/core/laminate.py:73`); float `==` thickness-uniformity test
  (`:229`).

**Workstream A effort:** ~2–4 focused PRs, each small and independently testable.

---

## 3. Workstream B — Dead code & documentation rot (P1)

A PyQt6 desktop GUI was removed in favor of the Streamlit app, but its
fingerprints are everywhere and several docs actively mislead users.

### B1. Broken and misleading files (delete or rewrite)
- `scripts/visual_matrix.py` — imports `PyQt6` and `bvidfe.gui.tabs.*`
  (`:27-34`), **neither of which exists**. The script cannot run. Delete, or
  rewrite against the Streamlit app if visual QA is still wanted.
- `docs/BUILD.md` — documents a PyInstaller `BVID-FE.exe`/`.app` pipeline,
  `bvidfe-gui` console script, `build_windows.bat`, `BvidFE.spec`, and
  `build-artifacts`/`release` workflows — **none exist**. Delete.
- `docs/python_api.md:139-151` — instructs users to
  `from bvidfe.viz.plots_3d import mesh_to_pyvista, plot_mesh_with_damage`;
  module and functions don't exist (actual module: `viz/plotly_3d.py`).
  Copy-pasting the published docs fails immediately.
- `ARCHITECTURE.md` — module catalog lists an entire `gui/` package
  (`main_window.py`, panels, tabs, workers, `config_io.py`), a PyVista
  `plots_3d.py`, and a `bvidfe-gui` entry point; the "Shipped in v0.2.0-dev"
  roadmap lists the PyQt GUI and PyInstaller packaging as shipped. Rewrite to
  match reality (Streamlit `app.py` is the GUI layer).
- `CONTRIBUTING.md` — `pip install -e ".[all]"` (`:25`): no `all` extra exists,
  so the very next step (`pytest`) fails on a clean checkout. `:89` tells
  contributors to subclass `FailureCriterion` from `failure/base.py`, which
  doesn't exist (contradicts the `add-failure-criterion` skill). `:76` hardcodes
  "149+ tests" (actual: 425).
- `CITATION.cff` — wrong version (`0.1.0`) and its `repository-code` disagrees
  with every other URL in the repo; README bibtex likewise stale.

### B2. Stale Qt-era references in living code
`viz/plots_2d.py:1-12,27-36` (design rationale written around a removed
`FigureCanvasQTAgg` embedding), `sweep/parametric_sweep.py:23` ("the GUI
`SweepWorker`"), `examples/README.md:25` and `examples/05_...py:1-3` ("the
GUI's File → Compare Tiers action"), `plotly_3d.py:17-20` (points to a
`mesh_3d_orthographic_figure` that exists nowhere).

### B3. Dead code to remove or wire up (decide per item)
- `stress_field_figure` (`viz/plotly_3d.py:164-190`) — no live caller; its data
  source `FieldResults` is documented as "currently always None"
  (`analysis/results.py:15-24`). Either populate `FieldResults` from the fe3d
  solve and add the promised "Stress Field" tab, or remove both.
- `_solve_failure_strain_analytic_scalar_ref` (`fe_tier.py:305-361`) — exists
  only for one equivalence test; move into the test module.
- `tsai_wu_strength_uniaxial` (`failure/tsai_wu.py:127-143`) — unreferenced.
- Scalar duplicates of vectorized code: `_aspect_ratio`/`_orientation_deg`/
  `_relative_size` (`impact/shape_templates.py:31-42`) and `_point_in_ellipse`
  (`analysis/fe_mesh.py:112-120`) — unused twins of the `_template_arrays` /
  `_points_in_ellipse` vectorized forms; drift hazards.
- `_CRITERION_NAMES` frozenset (`_types.py:42,51`) — third copy of the criterion
  name set; never consulted (see C2).
- `threshold_load(lam, pan, imp)` ignores `pan` and `imp`
  (`impact/olsson.py:210,232-234`).
- `plot_tier_comparison` (`viz/plots_2d.py:164-186`) — tested but unused by app
  or examples; wire into the app/example 05 (see E5) rather than delete.

### B4. Deployment-doc duplication
`DEPLOYMENT_STREAMLIT.md` (root) and `docs/deployment.md` are near-verbatim
copies. Keep the docs-site version; reduce the root file to a pointer.

### B5. CHANGELOG accuracy
`CHANGELOG.md` `[0.2.0]` claims "validated public datasets exercised in CI" —
contradicted by `validation/README.md` (only a tautological synthetic self-check
exists) and README's own Limitations section. Correct the entry; overclaiming
validation status is a credibility risk for a certification-adjacent tool.

### B6. GitHub hygiene
- **#129 is still open** although its body says "Closing — addressed by PR #130/#131." Close it.
- Six Dependabot PRs pending (#79, #80, #135, #136, #138, #139) — merge or batch.
- **PR #75** (per-ply thickness) is stale; `tests/test_per_ply_thickness.py`
  already exists on `main`, suggesting the feature landed another way. Triage &
  close or rebase.

**Workstream B effort:** 1–2 PRs of deletions/doc rewrites; near-zero risk, high
trust payoff. The ARCHITECTURE.md rewrite should land with C-workstream naming
decisions to avoid rewriting twice.

---

## 4. Workstream C — API & configuration hardening (P1)

### C1. `AnalysisConfig.to_dict()` / `from_dict()` in the library
Today config serialization is hand-rolled **three times**: `app.py:989-1039`
(write), `app.py:421-474` (read), `cli.py:248-275` (build). The library already
has `AnalysisResults.to_dict` and `damage_state_from_dict` but nothing for
config. Every new config field silently breaks the app's save/load and share-URL
paths with no test.
**Fix:** implement round-trip (de)serialization next to the dataclass, versioned
like the C-scan schema; migrate app and CLI onto it; regenerate
`docs/schemas/analysis_config.json` (the `regenerate-schemas` skill covers this).
**Accept:** round-trip property test (`from_dict(to_dict(cfg)) == cfg`); app and
CLI use the library path; a new-field regression test.

### C2. Make the failure criterion a first-class config knob (finishes #82's intent)
The fe3d criterion is hardcoded (`larc05` CAI at `fe_tier.py:375`, `tsai_wu` TAI
at `:501`); `AnalysisConfig` has no `criterion` field, so Puck and the
`FailureEvaluator` registry are unreachable from `BvidAnalysis.run()`. Meanwhile
the valid-criterion set lives in three places that must be hand-synced: the
`CriterionName` Literal, `_CRITERION_NAMES` frozenset (`_types.py`), and the two
parallel registries in `evaluator.py:20-36` (batch + scalar, "any extension
touches both" per its own comment).
**Fix:** add `cai_criterion`/`tai_criterion` (or one `criterion`) to
`AnalysisConfig` with validation at construction; derive the Literal-validated
set from a single registry that carries both scalar and batch callables (and,
later, the analytic-scaling strategy needed by A1). Update the
`add-failure-criterion` skill to the single-registration pattern.
**Accept:** an end-to-end fe3d run with each supported criterion; one place to
register a new criterion.

### C3. Deduplicate resolution/validation logic
- `_resolve_material` duplicated in `bvid.py:34-37` and `fe_tier.py:176-181`;
  both raise a bare `KeyError` for unknown materials while every other config
  field gets a friendly error. Consolidate + `ValueError` listing valid presets.
- "Pick the thinner sublaminate" logic copied three times
  (`semi_analytical.py:183-196`, `:388-393`, `fe_tier.py:443-446`) — all three
  must agree for cross-tier buckling consistency; extract one helper.
- `geometry.py:29,47` re-lists the `BoundaryKind`/`ImpactorShape` Literal values
  in tuple membership checks — derive from `typing.get_args`.
- `config.py:106-118` re-implements `laminate._normalise_ply_thicknesses`.

### C4. Encapsulation break: `elem._C_global` mutation
`fe_tier.py:171` reaches into `Hex8Element` internals to apply damage masks.
Add a public constructor parameter or `set_damage_mask()` on the element; this
also unblocks the ply-batched assembly in D2.

### C5. Calibration constants: centralize and document provenance
Uncalibrated magic constants are scattered across the impact chain:
`_SHAPE_DPA_FACTOR` (flat 1.4 / conical 0.7), diameter/mass exponents
(`mapping.py:43-73`), boundary bending factors 2.5/0.4 (`olsson.py:54-58`),
hardcoded 30° cone angle (`olsson.py:62`), steel-impactor contact constants
(`olsson.py:189-196`), 50%-thickness dent cap (`dent_model.py:30`), peanut
template constants (`shape_templates.py:32-57`), and the shared per-material
defaults (`material.py:46-64` — all four presets inherit identical
`olsson_alpha`/`soutis_k_s`/etc.). The README already flags calibration as a
limitation; the code should make it *visible and overridable*: gather these into
a documented calibration dataclass (per-material where physical), reference
their literature sources, and surface "default-calibration in use" in
`AnalysisResults.notes` for certification users. Feeds G6 (validation datasets).
Also note: `fiber_break_eta = 0.0` for every preset means the entire
fiber-break-core path is inert unless hand-overridden — either provide a
calibrated default or document it as opt-in.

### C6. Smaller API items
- CLI single-tier vs multi-tier silently switches JSON object ↔ array
  (`cli.py:304-307`) — always emit an envelope, or gate the array behind the
  multi-tier flag explicitly.
- Import-time side effects in `fe_tier.py:79-98` (logger handler install, env
  parse) — move to first-use.
- `MeshParams.cohesive_zone_factor` validated but never consumed
  (`config.py:47-74`) — deprecate until G1 actually consumes it.
- `_types.py` docstring promises a migration that should be finished (#82) and
  the issue closed.

---

## 5. Workstream D — Performance & scalability (P2)

The fe3d tier is bounded today by pure-Python assembly plus one-shot sparse LU,
guarded by a DOF cap (`FE3D_MAX_DOF = 500_000`, `fe_tier.py:98`) that proxies
for memory. Ordered by payoff:

### D1. Reuse the factorization and the element set
Every solve re-builds all element objects and re-factorizes from scratch:
`_build_elements` (`fe_tier.py:132-172`) and `solve_linear_static`
(`solver/static.py:38-45`) run once for CAI-FPF and again for TAI in the same
`fe3d_cai` call, and `sweep_energies` repeats it per energy point even though
only damage factors change. The `FeMeshSkeleton` cache (`fe_mesh.py:362-375`)
caches geometry only, holds exactly one entry, and is not thread-safe.
**Fix:** cache the assembled pristine `K` (and `scipy.sparse.linalg.splu`
factor) keyed by mesh signature; apply damage as a low-rank/element-subset
update where possible, or at minimum reuse the element list across channels;
make the skeleton cache a small LRU.
**Expected:** ~2× on a single fe3d run (CAI+TAI), much larger on sweeps.

### D2. Vectorize the assembly hot path
Per-element Python loops dominate: `assemble_coo` loops elements in Python
(`assembler.py:54-67`) calling `Hex8Element.stiffness_matrix`, which Python-loops
8 Gauss points (`hex8.py:207-209`); stress recovery FPF loops per element again
(`fe_tier.py:254-298`). On the regular grid, all elements of a ply share
`B`/`detJ` — `Ke` can be computed **once per (ply, damage-state) class** and
broadcast, and the Gauss loop is an `einsum` over a stacked `(n_gp, 6, 24)`
tensor the `GeometryTable` already holds. Same treatment for
`stress_at_gauss_points` and the (currently dead, but #132-relevant) `K_g`.
Also: `gauss.py:39-40` rebuilds the constant 8-point rule per call — memoize;
`hex8i` skips the geometry-table fast path entirely and computes each Jacobian
twice (`hex8i.py:78,96`).
**Expected:** order-of-magnitude assembly speedup; raises the practical mesh cap.

### D3. Memory-aware guard instead of a DOF proxy
COO buffers alone reach ~1.5 GB near the cap (`assembler.py:48-52`), and the cap
exists because a BLAS OOM is an uncatchable SIGSEGV (`fe_tier.py:91-97`).
Estimate bytes (COO + CSC + LU fill estimate) in `_guard_problem_size` and
compare against available memory (`psutil` optional dependency), keeping the DOF
cap as a fallback.

### D4. Sweep throughput
All sweeps are sequential loops (`parametric_sweep.py:236,260,284`).
`ProcessPoolExecutor` with the existing `on_error` semantics gives near-linear
speedup for fe3d sweeps. Add the missing sweep dimensions users hand-roll today
(impactor diameter, panel size, cartesian grids — exactly what
`examples/05` open-codes).

### D5. Micro-items
- Tsai-Wu batch criterion evaluated twice per element to fit the quadratic
  (`fe_tier.py:270-272`) — expose linear/quadratic parts from the criterion.
- `DamageState.projected_damage_area_mm2` runs a shapely `unary_union` on every
  access and is called inside the Brent iteration loop
  (`damage/state.py:62-68`, `shape_templates.py:102`) — cache on the instance.
- `fe_mesh` node generation is a Python list comprehension over all nodes
  (`fe_mesh.py:215-218`) — `meshgrid`; keep `element_dof_maps` as one `(n_elem,
  24)` array instead of a list of row arrays (`:253`) so the assembler can stop
  re-`asarray`ing each row.

**Sequencing note:** benchmark first (`scripts/benchmark.py` exists but is run
nowhere — wire it to a CI job or at least a `make bench` target so D-items have
before/after numbers).

---

## 6. Workstream E — User-facing layers (P2)

### E1. Break up the `app.py` monolith
1,608 lines containing pure logic (`_recommend_disposition`, `_build_ncr`,
`_parse_layup`, `_config_from_dict`, ~250 lines of NCR/MRB document generation)
that is untestable except via an import smoke test
(`tests/test_streamlit_app.py:35`). Extract the pure logic into an importable
package (e.g. `src/bvidfe/webapp/`), keep `app.py` as a thin page shell, and
unit-test the disposition banding and NCR assembly directly.
Specific defects to fix during extraction:
- Severity heatmap re-derives mesh dimensions and element ordering by hand
  (`app.py:1308-1312`) — duplicated from `estimate_fe_mesh_size`
  (`fe_mesh.py:75-88`) and the skeleton's internal C-ordering; if the library
  ever reorders elements the heatmap silently renders garbage. Move the
  projection into `viz/` with library-owned dimensions.
- Cache-hit toast is a wall-clock guess (`< 0.05 s` ⇒ "cache hit",
  `app.py:1106-1114`) — misleading; derive a real signal (e.g. a counter inside
  the cached function).
- `st.session_state["last_result"]` retains full `FeMesh` arrays alongside the
  `st.cache` copy (`app.py:1042-1109`) — double retention on the 1 GB Streamlit
  Cloud tier.
- `ply_thickness_mm` is the one primary input excluded from preset/share-URL
  state (`app.py:184-193,768-778`).

### E2. Close the CLI capability gap
The CLI cannot express: boundary condition (always `simply_supported`), impactor
shape (always hemispherical), impact location, mesh params (`elements_per_ply`,
`in_plane_size_mm` — fe3d always runs defaults with no DOF pre-check, unlike the
app), custom materials, file output (stdout only), or **any sweep** despite
`sweep_energies/layups/thicknesses` existing. Add `--boundary`,
`--impactor-shape`, `--impact-location`, `--elements-per-ply`,
`--in-plane-size`, `--output`, and a `bvidfe sweep` subcommand emitting the same
CSV schema as the library. C1's `from_dict` also enables `--config run.json`.

### E3. Unify visualization
Two unshared styling stacks: `style.py` palette is matplotlib-only;
`plotly_3d.py` hardcodes `"Viridis"`/`"Reds"`/`"Plasma"`; the app's severity
heatmap hardcodes `"hot_r"`. Define one palette/theme consumed by both; wire the
tested-but-unused `plot_tier_comparison` into the app; add the missing buckling
mode-shape visualization once #132 provides data; resolve `stress_field_figure`
per B3.

### E4. Examples & C-scan I/O polish
- `examples/05_tier_comparison_sweep.py:56-74` hand-rolls the energy loop
  instead of calling `sweep_energies` — the canonical example bypasses the
  canonical function and diverges from its CSV schema.
- `damage/io.py:66-69` rejects any `schema_version != "1.0"` exact-string — no
  minor-version tolerance policy despite the stated intent that old files keep
  loading. Define the policy now (accept `1.x` with warnings) before a `1.1`
  producer exists.
- Validation knows types but not physics: `interface_index` is unbounded above
  (`io.py:111-118`; interface 999 loads cleanly) and ellipses are never checked
  against panel extents. Add a `DamageState.validate_against(panel, n_plies)`
  step at analysis time, where the context exists.

---

## 7. Workstream F — Testing, CI, tooling (P1–P2)

### F1. Static typing — the biggest infra gap
The package ships `py.typed` (advertising typed APIs to downstream mypy users)
but **no type checker runs anywhere**. Adopt mypy (or pyright) with a lenient
baseline, ratchet strictness per-module, add the CI step next to ruff/black.
Fixes discovered along the way will overlap A1/C2 (Literal misuse is exactly
what a checker catches — issue #82's motivation).

### F2. Broaden lint & coverage configuration
- Ruff runs default rules only (`pyproject.toml:47-49`): enable `I` (isort),
  `B` (bugbear), `UP` (pyupgrade), `NPY`, `RUF`; `B011` alone would have flagged A5.
- No `[tool.coverage]` config exists: add `branch = true`, `exclude_lines`, and
  move the 90% floor from a CLI flag in `tests.yml:24` into versioned config.
- Register pytest markers (none exist): `slow` for fe3d cases, enabling a fast
  local `-m "not slow"` loop and unblocking #89's h-refinement test (which
  should land as designed: parametrized mesh sizes, plateau assertion).

### F3. Fill the specific test gaps
- `analysis/bvid.py` — the knockdown/`_pristine_strength` single source of truth
  has no direct unit test (only path tests); pin the thickness-weighted formula.
- `viz/plotly_3d.py` at 22% coverage; `solver/buckling.py` sparse path untested
  (A3); no Hex8i patch test (A6); `AnalysisResults.summary()/to_dict()`
  round-trip untested; no negative test that the validation gate *fails* on bad
  MAE (`validate_bvid_public.py:183-189` is unverified logic).
- Add `hypothesis` property tests for cross-tier invariants the README states
  as guarantees: `semi_analytical CAI ≤ empirical CAI`, knockdown monotone
  non-increasing in DPA, `0 < kd ≤ 1` (pairs with A2).

### F4. CI workflow improvements
- `cache: pip` on all three workflows (9 matrix legs × full scipy stack today).
- Unify action versions (`docs.yml` pins checkout@v4/setup-python@v5 vs v6
  elsewhere) — then merge the 4 pending Dependabot action PRs.
- Run `mkdocs build --strict` on PRs (docs currently only build on push to
  `main`, so link rot lands before it's caught — exactly how B1's dead API docs
  survived).
- Wire `scripts/generate_schemas.py --check` into CI (the script itself says
  it's intended for CI but was never wired), so `docs/schemas/*.json` can't
  drift from the dataclasses.
- Add a real Streamlit smoke test (`streamlit run --headless` + HTTP probe, or
  `streamlit.testing.v1.AppTest`) — both deployment docs list this as a TODO.
- Upload coverage to Codecov (or equivalent) for PR diffs; the artifact-only
  upload today means nobody sees coverage moves.
- Supply-chain: `pip-audit` on PRs + weekly (issue #117) and a CycloneDX SBOM on
  releases (issue #118) — both cheap, both already specified in the issues.

### F5. Packaging
- Split heavy UI deps: `streamlit` + `plotly` are hard deps of the *library*
  (`pyproject.toml:19-20`); move to an `[app]` extra (README/deploy docs adjust;
  `requirements.txt` remains the Streamlit Cloud manifest). Breaking for
  `pip install bvidfe` + `streamlit run app.py` users — coordinate with a minor
  release and CHANGELOG guidance.
- `requirements.txt` duplicates `pyproject` deps with no sync check — add a CI
  assertion or generate one from the other.

---

## 8. Workstream G — Physics roadmap (v0.3.x, aligns with existing issues)

Ordered by (user value ÷ effort), consistent with `ARCHITECTURE.md` "Planned"
and the open backlog:

1. **Observability quick wins (#133, #110).** Two strings for #133
   (closed-form-delegation note + `fe3d_buckling_governs` tag — directly
   requested, trivially cheap); per-stage solver logging for #110 (pairs
   naturally with D1/D2 instrumentation and gives Streamlit a live progress
   panel instead of a frozen spinner).
2. **Stiffness knockdown reporting (#114).** The damaged mesh already carries
   the per-element factors; homogenizing effective `E11/E22/G12` is mostly
   plumbing + a CLT baseline. High certification value (airworthiness modulus
   floors) for moderate effort.
3. **Validation datasets (G6).** Digitize Soutis (1996) and Caprino first —
   they anchor the two calibration constants users must trust most
   (`soutis_k_s`, `olsson_alpha`). Replace the tautological synthetic self-check
   as the only CI gate; keep fe3d advisory until (4) lands, then tighten.
   This is the single highest-credibility investment available: today the
   "validation gate" cannot fail for physics reasons, and the CHANGELOG already
   (incorrectly) claims dataset validation.
4. **Cohesive surfaces in fe3d (G1).** `elements/cohesive.py` is fully
   unit-tested but orphaned — imported nowhere in the analysis path; its tangent
   ignores the softening state (`cohesive.py:110-114`), so integration needs a
   consistent tangent + damage-state persistence first. This is the acknowledged
   fix for fe3d's energy-insensitivity (README Limitations; the biggest physics
   caveat in the tool). Large item; schedule as the v0.3.0 headline with A/D
   prerequisites done (D1/D2 make the nonlinear iterations affordable).
5. **Plate/shell buckling FE (#132).** Restores damage-zone heterogeneity that
   the closed-form delegation can't express. Depends on A3 (fix `buckling.py`'s
   contract first) and benefits from D2's vectorized assembly. The acceptance
   criteria in the issue (converge to Rayleigh-Ritz on pristine SSSS within 5%,
   measurable divergence on asymmetric damage, h-refinement regression per #89)
   are well-posed — keep them.
6. **Biaxial loading (#115).** Largest scope; the issue itself suggests a 2–3 PR
   split. Schedule after C2 (criterion plumbing) since the empirical tier's
   biaxial path runs through the criterion envelopes.
7. **LaRC05 full plane-search kinking** (currently a minimal Hashin-3D
   reduction, `larc05.py:1-8`) and per-material calibration of the C5 constants
   — as coupon data becomes available.

---

## 9. Suggested sequencing

| Phase | Contents | Outcome |
|---|---|---|
| **1. Patch series (v0.2.1)** | Workstream A (correctness), B1/B2/B5/B6 (rot deletions, CHANGELOG fix, issue hygiene), F4 cheap CI wins (pip cache, schema check, docs-on-PR, pip-audit #117) | Trustworthy baseline; docs stop lying; CI guards what it claims to guard |
| **2. API hardening (v0.2.2)** | C1 (config serialization), C2 (criterion knob, finish #82), C3/C4, E2 (CLI flags + sweep), F1 (mypy baseline), F2/F3 (markers, coverage config, missing tests, #89), #133/#110 quick wins | Library is the single source of truth; app/CLI are thin clients; typed API is actually checked |
| **3. Performance (v0.2.3)** | D1–D5 with benchmark baselines, E1 (app modularization), E3 (viz unification), F5 (dependency split) | fe3d usable at 2–4× current mesh density; sweeps parallel; app testable |
| **4. Physics (v0.3.0)** | G3 validation datasets → tightened gates, #114 stiffness knockdown, G1 cohesive surfaces, then #132 and #115 | The roadmap items the README already promises, on de-risked foundations |

Rationale for the order: correctness before features (A), credibility before
growth (B/validation), and the API consolidation (C1/C2) before the app/CLI work
that would otherwise triple-implement it again. Performance (D) precedes
cohesive surfaces (G1) because nonlinear cohesive iterations multiply the cost
of every inefficiency in the current assembly path.

---

## 10. Item → issue cross-reference

| Existing issue | Covered by |
|---|---|
| #82 stringly-typed API | C2, F1 |
| #89 h-refinement test | F2, F3 |
| #110 solver progress logging | G(1), D1 instrumentation |
| #114 stiffness knockdown | G(2) |
| #115 biaxial loading | G(6) |
| #117 pip-audit | F4 |
| #118 SBOM | F4 |
| #129 (self-described closed) | B6 — close it |
| #132 plate/shell buckling FE | G(5), prerequisite A3 |
| #133 fe3d buckling-channel tags | G(1) |

New issues should be opened for: A1 (Puck mis-dispatch), A2 (WN guard), A4
(version single-sourcing), A6 (Hex8i strict flag + fallback surfacing), B1
(doc-rot sweep), C1 (config serialization), C5 (calibration provenance), D1/D2
(fe3d performance), E1 (app modularization), E2 (CLI gaps), F1 (mypy), and the
validation-dataset digitization (G3).
