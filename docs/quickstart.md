# Quickstart

The fastest path from a fresh clone to a residual-strength estimate.

## 1. Install

```bash
git clone https://github.com/elhajjar1/BVID-FE.git
cd BVID-FE
pip install -e ".[dev]"
```

To also build this documentation site locally, install the `docs` extras:

```bash
pip install -e ".[docs]"
mkdocs serve     # http://127.0.0.1:8000
```

## 2. Try it in Jupyter

The fastest tour of the physics is
[`examples/quickstart.ipynb`](https://github.com/elhajjar1/BVID-FE/blob/main/examples/quickstart.ipynb)
— a runnable notebook that walks through:

- the Olsson quasi-static damage threshold,
- the peanut-template DPA distribution across ply interfaces,
- the Soutis CAI knockdown formula,
- the Whitney-Nuismer TAI point-stress criterion,

with LaTeX derivations and inline Plotly damage map / knockdown sweep / 3D
damage mesh visualisations.

```bash
pip install jupyter
jupyter notebook examples/quickstart.ipynb
```

## 3. Run from Python

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
        energy_J=20,
        impactor=ImpactorGeometry(diameter_mm=16),
        mass_kg=5.5,
    ),
    loading="compression",
    tier="semi_analytical",
)
result = BvidAnalysis(cfg).run()
print(result.summary())
```

This prints:

```text
BVID Analysis Results
  tier_used              : semi_analytical
  pristine_strength_MPa  : 601.3
  residual_strength_MPa  : 224.3
  knockdown              : 0.373
  dpa_mm2                : 1918.5
  dent_depth_mm          : 1.216
  n_delaminations        : 15
```

`knockdown` is `residual_strength_MPa / pristine_strength_MPa`, and
`dpa_mm2` is the projected damage area. The pristine strength is the
unnotched laminate strength from classical lamination theory. In this case
the Soutis bound (224.3 MPa) is below the delamination growth stress, so
the semi-analytical result equals the empirical one.

## 4. Run from the command line

```bash
bvidfe --material IM7/8552 \
       --layup "0,45,-45,90,90,-45,45,0" \
       --thickness 0.152 \
       --panel 150x100 \
       --loading compression \
       --energy 10
```

All lengths are in **millimetres**. `--panel` is `<Lx>x<Ly>` (literal `x`
separator, no spaces — e.g. `150x100` is `Lx = 150` mm, `Ly = 100` mm).
`--thickness` is the per-ply thickness in mm (one value, or one per ply);
`--layup` is a comma-separated list of ply angles in degrees; `--energy` is
the impact energy in joules. The command prints the full results JSON from
the default `empirical` tier. Add `--tier empirical,semi_analytical` to
compare tiers, `--quick-json` for a one-line result, or `--list-materials`
to see the six presets. `bvidfe --help` lists every option.

## 5. Streamlit web app

```bash
streamlit run app.py
```

Opens the BVID-FE web UI at <http://localhost:8501>. Configure the material,
layup, panel, impact (or C-scan), pick a tier, click **Run analysis**. Results
appear across the Summary / Damage Map / Knockdown Curve / 3D Damage /
Buckling / Damage Severity / Sweep / Export NCR tabs. Export NCR writes a
nonconformance report (Markdown or JSON) with a *proposed* disposition; a
Material Review Board must review and approve it.

See [Deployment](deployment.md) for publishing this app to a public URL on
Streamlit Community Cloud.

## 6. Compare with measured data

The validation harness runs a tier on every dataset in
`validation/datasets/` (three measured CAI datasets and a synthetic
self-check) and prints each prediction next to the test result:

```bash
python validation/validate_bvid_public.py --tier empirical
```

The empirical tier's mean absolute error is 1.5% on the NCAMP AS4/8552
coupons and 9.0% on the Lovejoy & Scotti IM7/8552 coupons. Both datasets
were used to fit `soutis_k_s` (and Lovejoy & Scotti's C-scan areas to fit
the default `olsson_alpha`), so expect errors of about 6–9% on a material
or layup outside them. See
[`validation/README.md`](https://github.com/elhajjar1/BVID-FE/blob/main/validation/README.md).

## More examples

The repository's [`examples/`](https://github.com/elhajjar1/BVID-FE/tree/main/examples)
folder contains five end-to-end scripts:

- `01_empirical_quick.py` — minimal one-shot empirical analysis.
- `02_tier_comparison.py` — all three tiers on the same damage state.
- `03_energy_sweep.py` — parametric sweep + plot.
- `04_inspection_driven.py` — C-scan JSON loading.
- `05_tier_comparison_sweep.py` — empirical vs. semi-analytical across
  impact energies.
