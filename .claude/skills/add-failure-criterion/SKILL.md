---
name: add-failure-criterion
description: Add a new failure criterion (e.g. Hashin, Christensen, max-stress) to the BVID-FE failure module. Use when the user wants to add, register, or wire up a failure index alongside the existing Tsai-Wu / LaRC05 / Puck criteria.
---

# Add a new failure criterion

Failure criteria plug into `FailureEvaluator` via a pair of registries.
A criterion is "added" only when it shows up in **four** places that
must move in lockstep:

1. The implementation module under `src/bvidfe/failure/`
2. `_CRITERION_REGISTRY` (batch) in `src/bvidfe/failure/evaluator.py`
3. `_CRITERION_SCALAR_REGISTRY` (scalar) in the same file
4. `CriterionName` Literal **and** `_CRITERION_NAMES` frozenset in
   `src/bvidfe/_types.py`

Forgetting any one of these causes a confusing runtime error far
from the actual omission. Treat it as one atomic change.

The fe3d tier does **not** go through `FailureEvaluator`: its
first-ply-failure solve in `src/bvidfe/analysis/fe_tier.py` imports the
criterion functions directly and has one load-scaling branch per
criterion (`_critical_strain_multiplier`, plus the scalar reference
`_solve_failure_strain_analytic_scalar_ref`). Its guard only checks the
name against `_CRITERION_NAMES`, so a new name with no branch there
**silently falls through to the Tsai-Wu `else` branch**. Either add the
branch (step 4) or tell the user the criterion is not available in fe3d.
`BvidAnalysis` does not expose a criterion choice: fe3d uses LaRC05 for
CAI and Tsai-Wu for TAI, so a new criterion is reached through
`FailureEvaluator` or the `criterion=` argument of
`_fe3d_cai_first_ply_failure`.

## Pre-flight

Confirm with the user:

1. The criterion name (the registry key, e.g. `"hashin"`). Use
   `snake_case`, lowercase. This is also the string users pass to
   `FailureEvaluator(criterion=...)`.
2. The literature source for the formulation.
3. Whether the criterion is direction-dependent or has any
   in-plane-only assumptions worth documenting.

## Required public interface

The new module must expose **two** functions with these exact
signatures, mirroring `larc05.py` / `puck.py` / `tsai_wu.py`:

```python
def <name>_index(material: OrthotropicMaterial, stress: np.ndarray) -> float:
    """Per-point failure index. stress shape: (6,) in Voigt order."""

def <name>_index_batch(material: OrthotropicMaterial, stresses: np.ndarray) -> np.ndarray:
    """Vectorised. stresses shape: (..., 6). Returns (...): leading axes kept."""
```

`FailureEvaluator` passes `(n_elem, n_gp, 6)`; the fe3d solve passes
`(n_gp, 6)` per element. `tests/failure/test_evaluator.py::
test_index_batch_preserves_leading_axes` checks the shape contract for the
existing criteria.

The batch form is **not optional** — `FailureEvaluator.evaluate` calls
it directly. The two functions must agree numerically.
`tests/failure/test_evaluator.py::test_evaluate_matches_scalar_loop`
checks this, but only for the criteria listed in its loop (Tsai-Wu and
LaRC05 today), so add the new criterion to that loop (step 5).

Voigt stress order in this codebase: `[σ11, σ22, σ33, σ23, σ13, σ12]`,
in the **ply material frame** (fe3d rotates each ply's stress into it with
`Hex8Element.stress_at_gauss_points_material` before calling a
criterion). Use the `*_resolved` properties for `Zt`, `Zc`, `S13`, which
fall back to `Yt`, `Yc`, `S12`. Confirm against any other criterion file
before assuming.

## Steps

### 1. Implement the criterion

Create `src/bvidfe/failure/<name>.py`. Use the closest existing
criterion as a structural template:

- `tsai_wu.py` — quadratic interaction, simplest reference
- `larc05.py` — fibre/matrix/kinking modes, decomposed
- `puck.py` — action-plane search, most complex

Vectorise the batch form properly (numpy operations on the full
`(n_elem, n_gp, 6)` array, not a Python loop) — the evaluator
delegates to it inside the hot path.

### 2. Register in `failure/evaluator.py`

Add the import and the two registry entries:

```python
from bvidfe.failure.<name> import <name>_index, <name>_index_batch

_CRITERION_REGISTRY[CriterionName] = {
    ...,
    "<name>": <name>_index_batch,
}

_CRITERION_SCALAR_REGISTRY[CriterionName] = {
    ...,
    "<name>": <name>_index,
}
```

Add the entry to **both** registries. The evaluator caches both
function pointers in `__init__`.

### 3. Extend the type alias in `_types.py`

In `src/bvidfe/_types.py`, add the new string to:

```python
CriterionName = Literal["tsai_wu", "larc05", "puck", "<name>"]

_CRITERION_NAMES: frozenset[str] = frozenset({"tsai_wu", "larc05", "puck", "<name>"})
```

The module comment explicitly calls out that the Literal and the
frozenset must move together — keep them in lockstep.

### 4. fe3d load-scaling branch (if fe3d should support it)

The fe3d first-ply-failure solve runs one static solve at a reference
strain and scales it analytically: stress is linear in the load
multiplier `c`, so the critical `c` follows from how the index scales.
Add an `elif criterion == "<name>":` branch, above the Tsai-Wu `else`, in
both `_critical_strain_multiplier` (vectorised) and
`_solve_failure_strain_analytic_scalar_ref` (scalar) in
`src/bvidfe/analysis/fe_tier.py`, using the existing ones as models:

- degree-2 homogeneous (LaRC05): `c = 1 / sqrt(idx_ref)`
- degree-1 homogeneous (Puck): `c = 1 / idx_ref`
- linear + quadratic terms (Tsai-Wu): solve `a c + b c^2 = 1` from two
  samples

Document in a comment why the scaling holds (no sign-selected branch may
flip as `c` grows). Then add the name to the `criterion` parametrize
lists in `tests/analysis/test_fpf_strain_solve_equivalence.py`.

### 5. Tests in `tests/failure/`

Create `tests/failure/test_<name>.py`. At minimum:

- A pristine-stress sanity check (index < 1 for a clearly safe state)
- A failure-stress sanity check (index ≥ 1 at the documented strength)
- A scalar-vs-batch equivalence test on a random stress field
- A check that `FailureEvaluator(material, criterion="<name>")`
  constructs and `.evaluate(...)` runs

Then update `tests/failure/test_evaluator.py`:
`test_criterion_registry_contains_known_keys` pins the exact registry key
set and fails until the new name is added there; add the new
`("<name>", <name>_index)` pair to the loop in
`test_evaluate_matches_scalar_loop`. Grep `tests/` for `"puck"` to find
any other list of all criteria.

### 6. Validate

```bash
ruff check src tests
black src tests
pytest tests/failure tests/analysis/test_fpf_strain_solve_equivalence.py -v
```

The full failure suite must stay green; the scalar-vs-batch test in
`test_evaluator.py` is the canonical guard against batch/scalar drift.

## Things this skill must never do

- Never add the criterion to only one of the two registries.
- Never extend `CriterionName` without also extending `_CRITERION_NAMES`.
- Never assume Voigt order — confirm it against an existing criterion
  file in this repo (it is `[σ11, σ22, σ33, σ23, σ13, σ12]`).
- Never skip the scalar-vs-batch test — it is what catches the most
  common implementation bug.
- Never leave a new name to fall through to the Tsai-Wu branch in
  `fe_tier.py` silently: add its scaling branch or say it is not
  supported in fe3d.
