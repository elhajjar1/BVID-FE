# BVID-FE Examples

Short, self-contained scripts that exercise the main BVID-FE workflows.
Each one runs standalone from the repository root:

```bash
python examples/01_empirical_quick.py
python examples/02_tier_comparison.py
python examples/03_energy_sweep.py
python examples/04_inspection_driven.py
python examples/05_tier_comparison_sweep.py
```

Outputs (PNGs, CSVs) are written to `examples/output/`, which is
gitignored; delete them when you are done.

## What each example shows

| Script | Demonstrates |
|---|---|
| `01_empirical_quick.py` | Minimal API call: 30-J impact on 8-ply IM7/8552 → CAI knockdown (about 0.31, DPA about 3,160 mm²) |
| `02_tier_comparison.py` | All three tiers (empirical / semi_analytical / fe3d) on the same 15-J impact, with a bar-chart summary (knockdown about 0.40, 0.40 and 0.26; at this damage area, about 1,530 mm², fe3d's 5 mm mesh is within about 5% of a 2 mm one, see below) |
| `03_energy_sweep.py` | `sweep_energies()` → DataFrame → CSV + knockdown-vs-energy plot (empirical, 5–40 J: about 0.55 at 5 J, just above the onset energy of about 0.8 J, then 0.45 at 10 J → 0.28 at 40 J) |
| `04_inspection_driven.py` | Loading a C-scan JSON (`sample_cscan.json`, 7 delaminations, about 1,040 mm²) and running the damage-driven path (no impact event); semi_analytical knockdown about 0.44 |
| `05_tier_comparison_sweep.py` | Empirical vs. semi_analytical over a 12-point energy sweep — overlaid knockdown curves and a CSV with both series. |
| `quickstart.ipynb` | Jupyter walkthrough of `01_empirical_quick.py` with LaTeX derivations (Olsson, Soutis, Whitney-Nuismer) and inline Plotly damage map / knockdown sweep / 3D damage mesh. Outputs are stripped on commit — run `jupyter notebook examples/quickstart.ipynb` to populate them. |

## Runtime expectations (default 5 mm mesh, one thread)

- empirical / semi_analytical: well under a second per run
- fe3d: about 5 s per run on the 8-ply 150 x 100 mm panel (it solves
  the damaged and the undamaged panel), so `02_tier_comparison.py` takes a
  few seconds. Thicker laminates and finer meshes take longer: about 35 s
  for the same panel at 2.5 mm, about 40 s for a 24-ply ASTM coupon at
  5 mm. For damage areas below about 1,000 mm² (a zone only a few 5 mm
  elements across) use `MeshParams(in_plane_size_mm=2.5)` or finer: at
  5 mm the fe3d knockdown can be about 25% high there. fe3d TAI is not
  recommended (far more conservative than the closed-form TAI; see the
  README's Limitations).
