# BVID-FE

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/elhajjar1/BVID-FE/actions/workflows/tests.yml/badge.svg)](https://github.com/elhajjar1/BVID-FE/actions/workflows/tests.yml)
[![Docs](https://img.shields.io/badge/docs-mkdocs--material-blue)](https://elhajjar1.github.io/BVID-FE/)

**Documentation: <https://elhajjar1.github.io/BVID-FE/>** (Quickstart, Python API, C-scan schema, physics models, deployment and the [JSON Schemas](https://elhajjar1.github.io/BVID-FE/schemas/analysis_config.json)).

BVID-FE predicts the residual compression (CAI) and tension (TAI) strength of fibre-reinforced composite laminates with Barely Visible Impact Damage (BVID). Give it a material, layup, panel and impact energy, or a C-scan of the damage, and it returns the knockdown: the residual strength divided by the strength of the undamaged laminate. Three modelling tiers trade runtime for fidelity, from closed-form estimates in tens of milliseconds to a 3D finite-element solution in seconds, through a Python API, a CLI and a Streamlit app. It is the third of a family of defect-specific composite tools, with **PorosityFE** (porosity) and **WrinkleFE** (fibre waviness).

## Which tier?

| Tier | Runtime | Validation status | Use it for |
|---|---|---|---|
| `empirical` | Tens of milliseconds | CAI calibrated on two measured datasets (NCAMP AS4/8552 and Lovejoy & Scotti IM7/8552; in-sample MAE 1.5% and 8.8%, about 5–9% cross-dataset). 22% on Hasebe et al.'s T800S/3900-2B coupons (reported, not gated). TAI not validated: its one constant is fitted to a single measured condition (Körbelin et al. 2020, 3.5% MAE in-sample). | Screening, energy sweeps, design trade studies |
| `semi_analytical` | Tens of milliseconds | Same CAI datasets and errors (the Soutis bound governs every one of their coupons); 24% on Hasebe et al. TAI identical to `empirical` (not validated). | Checking whether delamination growth governs (low energy, cross-ply laminates) |
| `fe3d` | About 4 s for an 8-ply 150 × 100 mm panel at the default 5 mm mesh (two FE solves), about 20 s at 2.5 mm; longer for thick laminates | **Qualitative.** Its damage-zone stiffness is fitted in-sample to three measured CAI coupons (Lovejoy & Scotti within 1.2%, NCAMP 21% low); not run on the measured datasets in CI. **TAI** keeps the zone's fibres intact (delaminations cut none) and takes the broken fibres from the closed-form equivalent hole, so it falls with energy 2–8% below the closed-form TAI: 3.6% MAE on the one measured TAI dataset, partly in-sample (the hole size is fitted to it). Use a 2.5 mm mesh or finer for damage areas below about 1,000 mm², and 1.25 mm for a zone under about 20 mm across in a coupon only 30–40 mm wide | Stress-field context and damage through the thickness; not for allowables |
| Biaxial (`transverse_stress_ratio ≠ 0`), `empirical` only so far | As the tier | **Experimental, not validated**: no biaxial residual-strength data exist to validate it. The pristine envelope is CLT ply discount along the load ray, the empirical CAI knockdown is assumed independent of the ratio, and TAI adds the Lekhnitskii lateral-stress term at the hole. `semi_analytical` and `fe3d` raise `NotImplementedError` for a nonzero ratio until follow-up PRs of #115 | How a lateral stress might move the residual strength; not for allowables |

All three report the knockdown against the same pristine strength; see
[Knockdown definition](https://elhajjar1.github.io/BVID-FE/physics_models/#knockdown-definition-and-cross-tier-comparability)
and [Validation](https://github.com/elhajjar1/BVID-FE#validation).

## Installation

BVID-FE is not on PyPI yet. Install it from GitHub:

```bash
pip install git+https://github.com/elhajjar1/BVID-FE.git                # library and the `bvidfe` command
pip install "bvidfe[app] @ git+https://github.com/elhajjar1/BVID-FE.git"  # also installs Streamlit, for the web app
```

Once 0.3.0 is on PyPI, `pip install bvidfe` (or `pip install "bvidfe[app]"`)
will do the same. From 0.3.0, each [GitHub release](https://github.com/elhajjar1/BVID-FE/releases)
carries a CycloneDX SBOM (`bvidfe-X.Y.Z.cdx.json`) of what `pip install bvidfe` installs.

The package contains the library and the CLI. The Streamlit app
(`app.py`), the example scripts and data (`examples/`, including the
`sample_cscan.json` used in the [Quickstart](https://elhajjar1.github.io/BVID-FE/quickstart/)), the Jupyter notebook and the validation
harness and datasets are in the repository, so clone it to use those:

```bash
git clone https://github.com/elhajjar1/BVID-FE.git
cd BVID-FE
pip install -e ".[dev]"       # tests, linters and Streamlit
pytest -v
```

## Quick start

```bash
bvidfe --material IM7/8552 --layup "0,45,-45,90,90,-45,45,0" \
       --thickness 0.152 --panel 150x100 --loading compression --energy 10 \
       --tier empirical,semi_analytical --quick-json
```

prints one JSON line per tier (full precision; rounded here):

```text
{"knockdown": 0.448, "residual_strength_MPa": 269.1, "pristine_strength_MPa": 601.3, "tier_used": "empirical"}
{"knockdown": 0.448, "residual_strength_MPa": 269.1, "pristine_strength_MPa": 601.3, "tier_used": "semi_analytical"}
```

Lengths are in mm (`--thickness` is per ply) and `--energy` in J; `bvidfe --help` lists every option. From Python, on a 16-ply panel at 20 J:

```python
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent

config = AnalysisConfig(
    material="IM7/8552",
    layup_deg=[0, 45, -45, 90] * 4,
    ply_thickness_mm=0.152,
    panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
    impact=ImpactEvent(energy_J=20, impactor=ImpactorGeometry(diameter_mm=16), mass_kg=5.5),
    loading="compression",
    tier="semi_analytical",
)
result = BvidAnalysis(config).run()
print(f"knockdown {result.knockdown:.3f}, DPA {result.dpa_mm2:.0f} mm²")  # knockdown 0.375, DPA 1811 mm²
```

The [Quickstart](https://elhajjar1.github.io/BVID-FE/quickstart/) covers C-scan import, parametric sweeps, the Streamlit app (`streamlit run app.py`) and the [Jupyter notebook](https://github.com/elhajjar1/BVID-FE/blob/main/examples/quickstart.ipynb).

## Validation

`validation/validate_bvid_public.py` scores a tier against measured residual strengths (mean absolute error, MAE):

| Dataset | Material and layup | Coupons | Measured strength (MPa) | Empirical MAE | Semi-analytical MAE |
|---|---|---|---|---|---|
| `ncamp_as4_8552_cai` (NCAMP CAM-RP-2010-002) | AS4/8552 [45/0/-45/90]3s, about 30 J | 7 | 171–182 | 1.5% (max 3.0%) | 1.5% (max 3.0%) |
| `lovejoy_scotti_im7_8552_cai` (NASA LaRC 2019) | IM7/8552 [45/0/-45/90]4s at 54 J; hard [45/0/90/0/-45/0/45/0/-45/0]2s at 68 J | 8 | 133–147; 236–249 | 8.8% (max 17.4%) | 8.8% (max 17.4%) |
| `hasebe_t800s_3900_2b_cai` (Hasebe et al. 2025), not gated | T800S/3900-2B, 8 layups (8–24 plies), 12.7–31.8 mm hemispheres, a 15.9 mm flat-ended impactor and a 60° cone at 2–21 J; 80×50 mm JIS K 7089 coupons | 343 | 81–592 | 21.5% (max 143.1%) | 24.0% (max 143.1%) |
| `korbelin_m21_t800s_tai` (Körbelin et al. 2020), TAI, not gated | M21/T800S [45/0/-45/90]2s at 25 J, 250×60 mm tension coupons (T800S/3900-2B card) | 5 | 884–1014 (tension) | 3.5% (max 9.8%), in-sample | 3.5% (max 9.8%) |

The NCAMP and Lovejoy & Scotti strengths are the predicted knockdown times the model's pristine strength, so the two errors can cancel. Scored as a knockdown against NCAMP's measured unnotched compression, the knockdown MAE is 3.6% (NCAMP) and 13.6% (Lovejoy & Scotti): the hard 50/40/10 laminate's pristine strength is 11.9% high and offsets a knockdown 17.3% too low. Hasebe and Körbelin are already scored on their measured pristine strength. The validator prints both numbers.

These errors are in-sample: `soutis_k_s` and the default `olsson_alpha` were fitted to the two gated datasets and `tai_notch_fraction` to the TAI one; fitted on one CAI dataset alone, `soutis_k_s` predicts the other with an MAE of 5.0–8.8%, so expect about 5–9% on a new material or layup. Sources, calibration and the fe3d results are in [`validation/README.md`](https://github.com/elhajjar1/BVID-FE/blob/main/validation/README.md).

## Limitations

- The calibration constants are fitted to two CAI datasets on the 8552 matrix (AS4/8552 and IM7/8552); the AS4/3501-6, T700/2510 and T800/epoxy presets have no measured residual-strength comparison. Calibrate against coupon data for the material before using results in certification.
- **TAI is not validated**: `tai_notch_fraction` is fitted to one measured condition (M21/T800S at 25 J), so its energy, layup and thickness trends are untested.
- **Biaxial loading is experimental and not validated** (`AnalysisConfig.transverse_stress_ratio`, `--transverse-stress-ratio`; `empirical` only): no biaxial CAI/TAI data exist, the CAI knockdown is assumed independent of the lateral stress, and under biaxial tension the pristine strength ends when every ply has cracked, which is conservative. See [Physics Models → Biaxial loading](https://elhajjar1.github.io/BVID-FE/physics_models/#biaxial-loading).
- On Hasebe et al.'s coupons the damage area matches the C-scans, but the damage-area-based knockdown misses the effect of layup.
- `fe3d` is qualitative: a component-wise stiffness reduction instead of cohesive surfaces, with a damage-zone factor (0.15) fitted in-sample to three CAI coupons. Its "LaRC05" is a quadratic Hashin-type index, not full LaRC05.
- The delamination growth stress is a 1D thin-film estimate, and the Soutis bound governs every coupon of the gated datasets, so the growth channel is effectively unvalidated.

The full list, with numbers, is in [Physics Models → Limitations](https://elhajjar1.github.io/BVID-FE/physics_models/#limitations).

## Citation

If you use BVID-FE in your research, please cite:

```bibtex
@software{elhajjar2026bvidfe,
  author    = {Elhajjar, Rani},
  title     = {{BVID-FE}: Barely Visible Impact Damage residual-strength analysis
               for composite laminates},
  year      = {2026},
  version   = {0.3.0},
  publisher = {GitHub},
  url       = {https://github.com/elhajjar1/BVID-FE},
  note      = {University of Wisconsin-Milwaukee}
}
```

## Contributing and license

See [CONTRIBUTING.md](https://github.com/elhajjar1/BVID-FE/blob/main/CONTRIBUTING.md) for reporting bugs, suggesting features and submitting pull requests. MIT License; see [LICENSE](https://github.com/elhajjar1/BVID-FE/blob/main/LICENSE).
