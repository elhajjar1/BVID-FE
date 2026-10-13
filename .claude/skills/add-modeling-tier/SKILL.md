---
name: add-modeling-tier
description: Add a new BVID-FE residual-strength modeling tier alongside empirical / semi_analytical / fe3d. Use when the user wants to add a new tier, fidelity level, or analysis backend, or extend the AnalysisConfig.tier dispatch.
---

# Add a new modeling tier

The three existing tiers (`empirical`, `semi_analytical`, `fe3d`) are
wired together by **four** core sync points, plus three hard-coded tier
lists in the user-facing tools, that must all change in one atomic edit. The dispatch in `BvidAnalysis.run` raises
`NotImplementedError` for unknown tiers, but the surrounding
`AnalysisConfig` validation, the type alias, and the runtime-validation
frozenset all catch the omission earlier and at different layers — get
any one of them wrong and the failure mode is confusing.

## Pre-flight

Confirm with the user:

1. The tier name (the string users pass to `AnalysisConfig(tier=...)`).
   Use `snake_case`, lowercase. Match the existing pattern
   (`empirical`, `semi_analytical`, `fe3d`).
2. Whether it supports **both** loading modes (`compression` and
   `tension`) or just one. `empirical` and `semi_analytical` support
   both via separate functions; `fe3d` splits compression into a
   delamination-growth channel + first-ply-failure channel and a tension
   channel.
3. What `BvidAnalysis.run` should populate on the returned
   `AnalysisResults`:
   - `buckling_eigenvalues` (set when the tier produces them, else `None`)
   - `critical_sublaminate` (set when a worst-interface index is known)
   - `field_results` (set when the tier produces an element/Gauss-point
     stress field)
   - `notes` / `warnings` (free-form provenance and warning tags)
4. How the tier's strength relates to the shared pristine strength.
   `run()` sets `pristine_strength_MPa` from `_pristine_strength` (the
   unnotched laminate strength, `failure/laminate_strength.py`) for every
   tier and `knockdown = residual / pristine`, so the residual must be on
   that scale. If the tier's own failure stresses are a different
   quantity (as fe3d's are), normalize by the tier's own undamaged run
   and multiply by `sigma_0`, as `BvidAnalysis._fe3d` does.

## The four sync points

| File | What changes |
|---|---|
| `src/bvidfe/analysis/<tier>.py` | New module implementing the tier's compression and/or tension paths |
| `src/bvidfe/_types.py` | Add the name to `TierName` Literal **and** `_TIER_NAMES` frozenset |
| `src/bvidfe/analysis/bvid.py` | New `elif self.config.tier == "<name>":` branch in `BvidAnalysis.run` |
| `tests/analysis/test_<tier>_path.py` | End-to-end test that exercises `BvidAnalysis(config).run()` with `tier="<name>"` |

The user-facing tools keep their own tier lists, which do not import
`TierName`:

| File | What changes |
|---|---|
| `src/bvidfe/cli.py` | `_VALID_TIERS` (validates `--tier`) |
| `app.py` | `TIERS` (Streamlit tier picker) |
| `validation/validate_bvid_public.py` | `--tier` `choices` and `_TIER_GATE_MULTIPLIER` (gate slack on `target_mae_pct`; unknown tiers fall back to 1.25) |

`scripts/validate_inputs.py` also has a `TIERS` list for its sensitivity
sweep; extend it if the new tier should be swept.

`AnalysisConfig` itself does not currently hold a `Literal` for tier —
it stores `tier: TierName` via the alias in `_types.py`. So extending
`TierName` is what extends the public typed surface.

## Steps

### 1. Implement the tier module

Create `src/bvidfe/analysis/<name>.py`. Reuse:

- `bvidfe.core.laminate.Laminate` — already constructed by `run()`
  before dispatch, passed in.
- `bvidfe.damage.state.DamageState` — already resolved (either user-
  supplied or from `impact_to_damage`).
- `_pristine_strength(lam, loading)` in `analysis/bvid.py` (the unnotched
  laminate strength, `failure/laminate_strength.unnotched_strength`) —
  already computed and passed in as `sigma_0`.

Mirror the closest existing tier:

- `analysis/semi_analytical.py` (returns a `SemiAnalyticalResult`
  dataclass with `residual_strength_MPa`, `critical_interface_index`,
  `critical_buckling_load_N`)
- `analysis/fe_tier.py` (free functions: `fe3d_cai`,
  `_fe3d_cai_first_ply_failure`, `fe3d_tai`, and `fe3d_cai_buckling` for
  buckling onset), wrapped by `BvidAnalysis._fe3d`, which runs them on the
  damaged and the undamaged mesh

### 2. Extend the type alias

In `src/bvidfe/_types.py`:

```python
TierName = Literal["empirical", "semi_analytical", "fe3d", "<name>"]
_TIER_NAMES: frozenset[str] = frozenset({"empirical", "semi_analytical", "fe3d", "<name>"})
```

The module comment explicitly says to keep these in lockstep. Both, or
neither.

### 3. Add the dispatch branch in `BvidAnalysis.run`

In `src/bvidfe/analysis/bvid.py`, add the `elif self.config.tier == "<name>":`
branch **above** the trailing `else: raise NotImplementedError(...)`.
Pattern by what the tier produces — example shapes:

```python
elif self.config.tier == "<name>":
    if self.config.loading == "compression":
        sigma, eigs, tier_notes = <name>_cai(...)
        notes.extend(tier_notes)
        buckling_eigs = eigs
    else:
        sigma = <name>_tai(...)
        buckling_eigs = None
    critical_interface = None
    field_results = None
```

The branch is responsible for populating every field that
`AnalysisResults(...)` reads at the bottom of `run()` — `sigma`,
`buckling_eigs`, `critical_interface`, `field_results`. Default to
`None` when the tier doesn't produce a given piece.

### 4. End-to-end tests

Create `tests/analysis/test_<name>_path.py`. At minimum:

- One CAI run (compression) returning a sane `residual_strength_MPa`
  in `(0, sigma_0]`.
- One TAI run (tension) if the tier supports it.
- An assertion that `results.tier_used == "<name>"`.
- An assertion that `results.config_snapshot["tier"] == "<name>"`
  (provenance round-trip).

Mirror the structure of `test_semi_analytical_path.py` or
`test_fe3d_path.py` — both run through `BvidAnalysis(config).run()`
rather than calling the tier functions directly, which is what locks
the dispatch wiring.

Also check the tests that loop over tiers (grep `tests/` for
`"semi_analytical"`): `test_per_ply_thickness.py`,
`test_cross_tier_consistency.py`, `test_edge_cases.py` and
`test_validation_smoke.py` loop over the closed-form tiers
(`empirical`, `semi_analytical`). Add the new tier where it is cheap
enough and the asserted property should hold for it.

### 5. JSON Schema and docs

The `AnalysisConfig` JSON Schema picks up the new tier through the
`Literal` automatically — regenerate with the `regenerate-schemas`
skill so `docs/schemas/analysis_config.json` stays in sync.

If the tier deserves narrative coverage, add a section to
`docs/physics_models.md` mirroring the existing tier sections, and a row
to the `analysis/` part of the module table in `ARCHITECTURE.md`.

### 6. Validate

```bash
ruff check src tests
black src tests
pytest tests/analysis -v
python scripts/generate_schemas.py     # regenerate (or use the skill)
python validation/validate_bvid_public.py --tier <name> \
    --dataset ncamp_as4_8552_cai lovejoy_scotti_im7_8552_cai
```

The last command scores the tier against the measured CAI datasets.
`synthetic_selfcheck` is generated from the empirical tier, so only the
measured datasets are a fair check of a new tier. CI does not run the new
tier until a step is added to `.github/workflows/tests.yml`; ask the user
before adding one.

## Things this skill must never do

- Never add the tier to `TierName` without also adding to
  `_TIER_NAMES`, or vice versa.
- Never add a new dispatch branch without an end-to-end test that
  goes through `BvidAnalysis(config).run()` — calling the tier
  function directly does not exercise the wiring.
- Never raise `NotImplementedError` for partial loading-mode support
  silently; populate `notes`/`warnings` and return a sane fallback, the
  way fe3d records in `notes` the raw stresses and the excluded panel
  buckling stress behind its knockdown. (The old `fe3d_buckling_fallback`
  warning tag is retired; no backend emits a `warnings` tag today.)
- Never widen `AnalysisResults` to add a new field for this tier
  without first checking whether an existing field
  (`notes`, `warnings`, `field_results`) already covers the need.
