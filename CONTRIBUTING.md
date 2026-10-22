# Contributing to BVID-FE

Thank you for your interest in contributing to BVID-FE!

## Reporting Bugs

Open an [issue](https://github.com/elhajjar1/BVID-FE/issues) with:
- Python version and OS
- Steps to reproduce the problem
- Expected vs actual behavior
- Error messages or tracebacks

## Suggesting Features

Open an issue describing:
- The use case (what composite damage analysis problem you are solving)
- Expected behavior and outputs
- Any relevant references or equations

## Development Setup

```bash
git clone https://github.com/elhajjar1/BVID-FE.git
cd BVID-FE
pip install -e ".[dev]"
pytest tests/ -v
```

## Running Tests

```bash
# Full suite
pytest tests/ -v

# Single subsystem
pytest tests/impact/ -v
pytest tests/failure/ -v
pytest tests/analysis/ -v

# What CI's test job runs (coverage floor 90%)
pytest -v --cov=bvidfe --cov-report=term --cov-fail-under=90
```

Coverage also measures subprocesses (`[tool.coverage.run] patch =
["subprocess"]` in `pyproject.toml`, needs coverage >= 7.10), so the CLI
tests that run `python -m bvidfe.cli` count towards it.

## Code Style

- Format with **black** and lint with **ruff** (pre-commit hooks enforce both)
- Add docstrings to all public functions and classes
- Include units in variable names and docstrings (MPa, mm, rad, J)
- Use SI units throughout (millimetres, MPa, kg)

## Pre-commit setup

Install pre-commit hooks once after cloning so black, ruff, and basic file
hygiene checks run automatically before each commit:

```bash
pip install pre-commit
pre-commit install
```

### Before committing

If you haven't installed pre-commit, run these three commands locally before
pushing — CI runs `black --check` and **will fail** if any file is not
already formatted, so you must auto-format with `black` (not just verify
with `black --check`):

```bash
ruff check src tests
black src tests        # auto-format in place; do NOT use --check here
pytest -q
```

## Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Write tests for any new functionality before implementing it
4. Ensure the full suite passes (`pytest tests/ -v`), and run the
   validation harness if you changed a model (see below)
5. Submit a pull request with a clear description of what was changed and why

## Examples

Runnable examples live in `examples/`. The script set (`01_empirical_quick.py` ... `05_tier_comparison_sweep.py`) covers the headless workflows; `examples/quickstart.ipynb` is an annotated Jupyter notebook that mirrors `01_empirical_quick.py` with LaTeX derivations of the Olsson, Soutis, and Whitney-Nuismer models plus inline Plotly visualisations. **Strip notebook outputs before committing** with `nbstripout examples/quickstart.ipynb` (`pip install nbstripout`); commits with output cells will be flagged in review.

The ```` ```python ```` blocks in `README.md`, `docs/quickstart.md` and `docs/python_api.md` are executed by `tests/test_docs_examples.py`, in order and sharing one namespace per file, so a later block can use a name an earlier one defined. Quote printed numbers rounded so small model changes do not stale them. Mark a block that cannot run on its own (a fragment, or one that needs fe3d or Streamlit) with `<!-- docs-test: skip (reason) -->` on the line before its fence.

## Adding New Materials

Add a new entry to `MATERIAL_LIBRARY` in `src/bvidfe/core/material.py` using the `OrthotropicMaterial` dataclass, with `name` equal to the dict key. Include all orthotropic stiffness constants (MPa), strength values (MPa), interlaminar fracture toughnesses (N/mm), and density (kg/mm^3, e.g. `1.58e-6`). Cite a literature source for every constant in a comment above the entry (see the `AS4/8552` entry). Leave the calibration constants (`soutis_k_s`, `olsson_alpha`, ...) at their defaults unless you have material-specific test data. Add tests in `tests/core/test_material.py` and add the name to the preset lists in `test_material_library_has_six_presets` and `tests/test_cli_version.py`. The step-by-step procedure is in `.claude/skills/add-material/SKILL.md`.

## Adding New Failure Criteria

Create `src/bvidfe/failure/<name>.py` with a scalar `<name>_index(material, stress)` and a vectorised `<name>_index_batch(material, stresses)` (ply-frame Voigt stress `[s11, s22, s33, t23, t13, t12]`). Register both in `_CRITERION_REGISTRY` and `_CRITERION_SCALAR_REGISTRY` in `failure/evaluator.py`, and add the name to `CriterionName` and `_CRITERION_NAMES` in `src/bvidfe/_types.py`. If the fe3d tier should be able to use it, add a load-scaling branch in `_critical_strain_multiplier` (and the scalar reference) in `analysis/fe_tier.py`; otherwise the name falls through to the Tsai-Wu branch there. Add tests in `tests/failure/`. Details: `.claude/skills/add-failure-criterion/SKILL.md`.

## Adding New Modeling Tiers

1. Implement the tier in `src/bvidfe/analysis/my_tier.py`
2. Add the tier string to `TierName` and `_TIER_NAMES` in `src/bvidfe/_types.py` (`AnalysisConfig.tier` is typed `TierName`)
3. Add the dispatch branch in `BvidAnalysis.run()` (`analysis/bvid.py`); the residual must be on the same scale as `pristine_strength_MPa`
4. Add the name to the hard-coded tier lists in `cli.py` (`_VALID_TIERS`), `app.py` (`TIERS`) and `validation/validate_bvid_public.py` (`--tier` choices, `_TIER_GATE_MULTIPLIER`)
5. Add end-to-end tests in `tests/analysis/test_my_tier_path.py` and regenerate the JSON schemas (`python scripts/generate_schemas.py`)

Details: `.claude/skills/add-modeling-tier/SKILL.md`.

## Adding a Validation Dataset

Validation datasets live in `validation/datasets/<name>.json`; the harness
(`validation/validate_bvid_public.py`) discovers every `*.json` there.

1. Write the file with a top-level `name`, `reference` (full citation and
   any assumptions), `target_mae_pct` and `cases`. Each case needs
   `material`, `layup_deg`, `ply_thickness_mm`, `panel_Lx_mm`,
   `panel_Ly_mm`, `impactor_diameter_mm`, `impactor_mass_kg` (no defaults:
   if the source omits them, state the assumed values in `reference`),
   `impact_energy_J` and `measured_strength_MPa`. Optional: `loading`
   (default `compression`), `boundary`, `impactor_shape`,
   `measured_dent_mm`, `measured_dpa_mm2`, `note`.
2. Run it: `python validation/validate_bvid_public.py --dataset <name>`
   (and `--tier semi_analytical`).
3. If the model is known to miss the dataset, set `"gate": false` with a
   `gate_note` explaining why: its error is still printed but cannot fail
   `--gate`.
4. Do not tune calibration constants against the dataset you validate on
   without checking them on held-out data; see `validation/README.md`.
5. For a measured dataset, also add its name to the semi-analytical step
   in `.github/workflows/tests.yml` and to the parametrize list of
   `test_measured_cai_gates_pass_for_the_closed_form_tiers` in
   `tests/test_validation_smoke.py`.

## Continuous Integration

`.github/workflows/tests.yml` runs on every push and pull request:

- **test** (Ubuntu, macOS and Windows; Python 3.10, 3.11, 3.12):
  `ruff check src tests`, `black --check src tests`, and
  `pytest` with a 90% coverage floor.
- **validation** (Ubuntu, Python 3.11, after **test** passes):
  - `python scripts/generate_schemas.py --check`: the files in
    `docs/schemas/` must match the dataclasses. Regenerate them with
    `python scripts/generate_schemas.py` in the commit that changes a field.
  - Empirical tier on every dataset, `--gate` (hard gate).
  - Semi-analytical tier on the measured datasets
    (`ncamp_as4_8552_cai`, `lovejoy_scotti_im7_8552_cai`,
    `hasebe_t800s_3900_2b_cai`), `--gate` (hard gate; the Hasebe set is
    reported only, `"gate": false`). `synthetic_selfcheck` is generated from the empirical tier, so
    other tiers cannot match it.
  - fe3d tier on the first 3 cases of `synthetic_selfcheck` only, as a
    smoke run without `--gate` (`continue-on-error`).

`--gate` fails when a dataset's MAE exceeds its `target_mae_pct` times
1.25 (empirical, semi-analytical) or 2.0 (fe3d), unless the dataset sets
`"gate": false`.

`.github/workflows/docs.yml` builds the MkDocs site with `mkdocs build
--strict` and deploys it on every push to `main`.

## Claude Code on the web

This repo ships a `SessionStart` hook at `.claude/hooks/session_start.sh` (wired
up via `.claude/settings.json`) that runs `black --check src tests` and
`ruff check src tests` at the start of every Claude Code session and reports
any drift. It is informational only (always exits 0) and prevents agents from
delegating commits while formatting is broken. If `black`/`ruff` are not
installed yet, the hook prints a friendly note and exits cleanly.

## Releasing

BVID-FE publishes to PyPI via **OIDC Trusted Publishing** — there is no
long-lived API token stored in the repo. The release flow is driven by
`.github/workflows/publish.yml`, which triggers on any tag matching `v*`.

### One-time PyPI setup (maintainer only)

Trusted Publishing must be configured **once** on PyPI for this repo. See
[PyPI's Trusted Publishing docs](https://docs.pypi.org/trusted-publishers/)
for the canonical instructions. Register a pending or existing publisher
with:

- **Repository owner:** `elhajjar1`
- **Repository name:** `BVID-FE`
- **Workflow filename:** `publish.yml`
- **Environment name:** `release`

A matching GitHub Actions environment named `release` must also exist on
the repo (Settings → Environments → New environment → `release`). The
workflow declares `environment: release`, so it will queue indefinitely
until that environment is created — which is the intended fail-safe.

### Cut-a-release flow

1. Bump the version (e.g. `0.2.0` → `0.2.1` for a patch, `0.3.0` for a
   feature release) in all four places: `[project].version` in
   `pyproject.toml`, `__version__` in `src/bvidfe/__init__.py`,
   `version` and `date-released` in `CITATION.cff`, and the BibTeX
   `version = {...}` line in `README.md`. `tests/test_package.py` fails if
   they disagree.
2. Update `CHANGELOG.md`: move the `## [Unreleased]` notes into a
   new `## [X.Y.Z] - YYYY-MM-DD` section (today's date) and leave a
   fresh empty `## [Unreleased]` above it.
3. Commit those changes on `main` and wait for CI to pass; the publish
   workflow does not run the tests.
4. Tag the commit with `git tag vX.Y.Z` (note the leading `v`) and
   push with `git push origin vX.Y.Z`.
5. The `Publish to PyPI` workflow takes over: it verifies the tag
   matches `pyproject.toml`, builds an sdist + wheel, runs
   `twine check --strict`, and uploads to PyPI through OIDC. No manual token
   handling.

If the tag/version check fails, delete the tag (`git push --delete
origin vX.Y.Z` and `git tag -d vX.Y.Z`), fix the version mismatch,
and re-tag.

## Code of Conduct

Be respectful. Focus feedback on code and ideas, not people. Contributions at any experience level are welcome.
