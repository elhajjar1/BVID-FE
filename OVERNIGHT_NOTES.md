# Overnight session notes — 2026-04-17 → 2026-04-18

## TL;DR

- **217 tests passing** (up from 197 at start of night)
- `dist/BVID-FE.app` rebuilt with all improvements (quarantine cleared)
- Launch from terminal for visible progress logs:
  ```bash
  ./dist/BVID-FE.app/Contents/MacOS/BVID-FE
  ```

## What shipped overnight

| Commit | What |
|---|---|
| `0c2477e` | Overnight improvements: auto-knockdown-curve after every single run, two placeholder tabs replaced with real matplotlib content, 14 edge-case tests |
| `85f0d02` | Vectorized `Hex8Element.shape_derivatives` and `B_matrix` (removed 8-node Python loops) |
| `ce4c9d4` | New **File → Compare Tiers** action: overlays empirical + semi_analytical knockdown curves |
| `4566dac` | `examples/05_tier_comparison_sweep.py` — Python-API equivalent of the Compare Tiers action |
| `88763bb` | `bvidfe --version` CLI flag + GUI **Help** menu with About dialog and tier-interpretation guide |
| `377c8ab` | `--version` flag smoke test (217 tests) |

## New GUI tabs (replaced placeholders)

- **Buckling Eigenvalues** — bar chart of the first up to 6 buckling eigenvalues from `semi_analytical` or `fe3d` runs. Static note for empirical.
- **Damage Severity** — top-down heatmap of the panel, colored by number of delaminated interfaces stacked at each (x, y) column. Analogous to a C-scan plan view.

## Knockdown Curve now auto-populates

After every single-run analysis, a background empirical-tier energy sweep fires (8 points, sub-second) and fills the Knockdown Curve tab. No need to click "Run energy sweep" explicitly.

## Compare Tiers (new menu item)

**File → Compare Tiers (empirical + semi_analytical)…** runs both tiers on the current config across 8 energies and overlays them on the Knockdown Curve tab. Useful for quickly understanding whether the buckling contribution is meaningful for your panel.

## Help menu (new)

**Help → About BVID-FE…** shows version, license, repo link.
**Help → How to interpret the tiers…** short comparison table of when to use each tier.

## fe3d buckling refactor — attempted, shelved

Tried to fix the flat-knockdown problem by wiring proper in-plane pre-stress BCs + Cauchy-stress recovery into `fe3d_cai_buckling`. The refactor compiled and all existing tests still passed, but the pristine buckling result was 14× below the analytical plate-buckling formula. Likely a sign-convention issue with how `K_g` interacts with compressive pre-stress, but a proper fix needs more time than a single session. Reverted to the known-good (conservative, FPF-dominant) buckling path.

A graded-damage implementation was also tried (sample element footprints for fractional ellipse overlap → smooth damage_factor between 1.0 and 1e-4). Reverted because it exposed the same buckling sign issue and perversely made "damaged" sometimes report higher strength than "pristine".

Both experiments are in git history (reverted cleanly) if you want to re-explore with fresh eyes.

## fe3d limitation — documented, not fixed

- `fe3d` knockdown stays approximately flat vs. impact energy above the Olsson threshold (stress concentration at healthy/damaged boundary dominates FPF strain; buckling path gets rejected by 5%-pristine sanity check).
- Summary tab now shows an inline note pointing users to `empirical` or `semi_analytical` for energy-dependent curves.
- CHANGELOG's Known Limitations and README's Limitations both explicitly call this out.
- v0.3.0 fix requires proper in-plane pre-stress BCs + cohesive surface elements.

## Performance

- fe3d single run on 150×100 panel with default mesh: ~10 s per run
- Empirical sweep (8 energies): <1 s
- Semi_analytical sweep (8 energies): ~1 s
- Compare Tiers action: ~2 s

## If you want to keep iterating

1. **Calibrate material constants** against one digitized dataset (Soutis / Caprino / Sanchez-Saez). Paste me `(energy, measured CAI strength)` pairs and I'll tune `olsson_alpha / soutis_k_s / dent_beta` to minimize MAE.
2. **Real fe3d buckling** — the refactor I shelved is the right idea; it needs a proper sign/BC audit and an analytical plate-buckling reference case to calibrate against.
3. **Cohesive surfaces** instead of stiffness reduction — would fix both the flat knockdown and the unphysical stress concentrations. Substantial v0.3.0 scope.
4. **Icon + code signing** for the .app bundle so Gatekeeper doesn't prompt on download.

Nothing in this session broke anything. Full test suite passes; the worst-case rollback is `git reset --hard` to any commit from before last night.
