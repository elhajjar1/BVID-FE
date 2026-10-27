# BVID-FE

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/elhajjar1/BVID-FE/actions/workflows/tests.yml/badge.svg)](https://github.com/elhajjar1/BVID-FE/actions/workflows/tests.yml)
[![Docs](https://img.shields.io/badge/docs-mkdocs--material-blue)](https://elhajjar1.github.io/BVID-FE/)

> **Documentation:** the full MkDocs site is published at
> <https://elhajjar1.github.io/BVID-FE/> and includes a Quickstart, the
> Python API reference, the C-scan JSON schema, the physics models, and the
> Streamlit deployment guide. Auto-generated JSON Schemas for `AnalysisConfig`
> and `AnalysisResults` live at
> [`/schemas/`](https://elhajjar1.github.io/BVID-FE/schemas/analysis_config.json).

A Python library for predicting residual strength and stiffness of fiber-reinforced composite laminates containing Barely Visible Impact Damage (BVID).

## Why This Tool?

Low-velocity impacts on composite structures can create internal delaminations that are invisible to the naked eye yet significantly degrade compression and tension strength. Engineers certifying composite airframes and pressure vessels need quick estimates of how much strength is lost and how the knockdown (residual strength divided by undamaged strength) depends on layup, panel size, and impact energy. BVID-FE provides three modeling tiers — from closed-form estimates that run in tens of milliseconds to a 3D finite element solution that takes tens of seconds — so the level of fidelity can be chosen for each stage of the design process.

BVID-FE is the third in a family of defect-specific composite tools, joining **PorosityFE** (porosity defects) and **WrinkleFE** (fiber waviness). The three tools share material models, laminate theory, failure criteria, and documentation conventions.

## Features

- **Two workflow paths** converging on a shared `DamageState`:
  - *Impact-driven*: Olsson quasi-static threshold + energy-balance damage area + peanut-template distribution over the interfaces + empirical dent model
  - *Inspection-driven*: C-scan JSON import per the documented [schema](https://elhajjar1.github.io/BVID-FE/cscan_schema/)
- **Three modeling tiers** for residual strength after BVID:
  - *Empirical*: Soutis CAI knockdown + Whitney-Nuismer TAI on an equivalent hole for the fibre damage (tens of milliseconds)
  - *Semi-analytical*: buckling-driven delamination growth (Rayleigh-Ritz sublaminate buckling + thin-film energy release rate + Benzeggagh-Kenane mixed-mode toughness, weakest interface governs), capped by Soutis; the empirical TAI model (tens of milliseconds)
  - *3D FE*: first-ply failure on a hexahedral mesh with the damage zone modelled as a soft inclusion, normalised by an undamaged fe3d run; LaRC05 for CAI, Tsai-Wu for TAI (tens of seconds; see [Which tier?](https://github.com/elhajjar1/BVID-FE#which-tier))
- **CAI and TAI loading modes** (Compression-After-Impact and Tension-After-Impact)
- **Pristine (unnotched) laminate strength** from classical lamination theory with ply-discount last-ply failure, checked against NCAMP unnotched test data
- **Per-interface ellipse damage model** using `DelaminationEllipse` with shapely-union projected damage area
- **Six material presets**: AS4/3501-6, AS4/8552, IM7/8552, T700/2510, T800S/3900-2B, T800/epoxy
- **Validation harness** (`validation/`) that scores the tiers against three measured CAI datasets (NCAMP AS4/8552 and Lovejoy & Scotti IM7/8552, gated in CI; Hasebe et al. T800S/3900-2B, reported only) and one TAI condition (Körbelin et al. M21/T800S, reported only); see [Validation](https://github.com/elhajjar1/BVID-FE#validation)
- **CLI** for single runs, several tiers in one call (`--tier empirical,semi_analytical`), and shell pipelines (`--quick`, `--quick-json`)
- **Streamlit web app** with sidebar-driven `AnalysisConfig` and result tabs — Summary, Damage Map, Knockdown Curve (live empirical sweep around the current energy), 3D Damage (Plotly hex mesh coloured by damage factor), Buckling (tier-specific buckling indicator: `N_cr` in N/mm for `semi_analytical`, `σ_crit` in MPa for `fe3d`), Damage Severity (through-thickness sum of per-element `1 − damage_factor`; see [Damage severity heatmap](https://github.com/elhajjar1/BVID-FE#damage-severity-heatmap)), Sweep (CSV-exportable parametric energy sweep), and Export NCR (a nonconformance report in Markdown or JSON with a *proposed* disposition that a Material Review Board must review). Cached runs, fe3d mesh-size guards, JSON/CSV downloads
- **Parametric sweeps** over impact energy, layup, or ply thickness with CSV output
- **2D plots**: damage map, knockdown curves, tier comparison charts
- **Plotly 3D viz** for the damaged hex mesh (`bvidfe.viz.plotly_3d`)

## Which tier?

| Tier | Runtime | Validation status | Use it for |
|---|---|---|---|
| `empirical` | Tens of milliseconds | CAI calibrated on two measured datasets (NCAMP AS4/8552 and Lovejoy & Scotti IM7/8552; in-sample MAE 1.5% and 8.8%, about 5–9% cross-dataset). 22% on Hasebe et al.'s T800S/3900-2B coupons (reported, not gated). TAI not validated: its one constant is fitted to a single measured condition (Körbelin et al. 2020, 3.5% MAE in-sample). | Screening, energy sweeps, design trade studies |
| `semi_analytical` | Tens of milliseconds | Same CAI datasets and errors (the Soutis bound governs every one of their coupons); 23% on Hasebe et al. TAI identical to `empirical` (not validated). | Checking whether delamination growth governs (low energy, cross-ply laminates) |
| `fe3d` | About 5 s for an 8-ply 150 × 100 mm panel at the default 5 mm mesh (two FE solves), about 35 s at 2.5 mm; longer for thick laminates | **Qualitative.** Its damage-zone stiffness is fitted in-sample to three measured CAI coupons (Lovejoy & Scotti within 1.2%, NCAMP 21% low); not run on the measured datasets in CI. **TAI not recommended:** 62% low on the one measured TAI dataset. Use a 2.5 mm mesh or finer for damage areas below about 1,000 mm² | Stress-field context and damage through the thickness; not for allowables |

All three report the knockdown against the same pristine strength; see
[Knockdown definition](https://github.com/elhajjar1/BVID-FE#knockdown-definition-and-cross-tier-comparability)
and [Validation](https://github.com/elhajjar1/BVID-FE#validation).

## Installation

From PyPI:

```bash
pip install bvidfe            # library and the `bvidfe` command
pip install "bvidfe[app]"     # also installs Streamlit, for the web app
```

The package contains the library and the CLI. The Streamlit app
(`app.py`), the example scripts and data (`examples/`, including the
`sample_cscan.json` used below), the Jupyter notebook and the validation
harness and datasets are in the repository, so clone it to use those:

```bash
git clone https://github.com/elhajjar1/BVID-FE.git
cd BVID-FE
pip install -e ".[dev]"       # tests, linters and Streamlit
pytest -v
```

## Quick Start

### Try it in Jupyter

The fastest tour of the physics is [`examples/quickstart.ipynb`](https://github.com/elhajjar1/BVID-FE/blob/main/examples/quickstart.ipynb) — a runnable notebook that walks through the Olsson damage threshold, peanut-template DPA distribution, Soutis CAI, and Whitney-Nuismer TAI with LaTeX derivations and inline Plotly damage map / knockdown sweep / 3D damage mesh visualisations. Open it with `jupyter notebook examples/quickstart.ipynb` (from a clone, after `pip install -e . jupyter`).

### Streamlit web app

```bash
streamlit run app.py
```

Opens the BVID-FE web UI in your browser: configure material + layup, panel size, impact event or C-scan inspection, pick a tier, click **Run analysis**. Results appear across the Summary / Damage Map / Knockdown Curve / 3D Damage / Buckling / Damage Severity / Sweep / Export NCR tabs. Use the JSON / CSV download buttons to export results. See [DEPLOYMENT_STREAMLIT.md](https://github.com/elhajjar1/BVID-FE/blob/main/DEPLOYMENT_STREAMLIT.md) for deploying a public URL on Streamlit Community Cloud.

### Command-line interface

```bash
bvidfe --material IM7/8552 \
       --layup "0,45,-45,90,90,-45,45,0" \
       --thickness 0.152 \
       --panel 150x100 \
       --loading compression \
       --energy 10
```

All lengths are in **millimetres**. `--panel` is `<Lx>x<Ly>` (in-plane
dimensions, x then y, literal `x` separator, no spaces — e.g. `150x100`
is Lx = 150 mm, Ly = 100 mm). `--thickness` is the per-ply thickness in
mm (one value, or a comma-separated list with one value per ply), `--layup`
is a comma-separated list of ply angles in degrees, and `--energy` is the
impact energy in joules. The default tier is `empirical` and the output is
the full results JSON. Other options: `--tier` (one tier or a
comma-separated list), `--cscan scan.json` in place of `--energy`,
`--impactor-diameter` and `--mass` (defaults 16 mm and 5.5 kg), `--quick` or
`--quick-json` for a one-line result, and `--list-materials`. For example:

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

### Python API — impact-driven path

```python
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry
from bvidfe.impact.mapping import ImpactEvent

impact_cfg = AnalysisConfig(
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
result = BvidAnalysis(impact_cfg).run()
print(result.summary())
```

```text
BVID Analysis Results
  tier_used              : semi_analytical
  pristine_strength_MPa  : 601.3
  residual_strength_MPa  : 225.4
  knockdown              : 0.375
  dpa_mm2                : 1810.9
  dent_depth_mm          : 0.219
  n_delaminations        : 15
```

DPA is the projected damage area in mm². Here the Soutis bound (225.4 MPa)
is below the delamination growth stress (306.9 MPa at interface 11), so the semi-analytical result
equals the empirical one. No damage is predicted below the onset energy,
about 3.3 J for this panel and impactor. The per-interface ellipses are in
`result.damage.delaminations`, the governing interface in
`result.critical_sublaminate`, and the runtime notes in `result.notes`.

### Inspection-driven path (C-scan import)

```python
from bvidfe.damage.io import load_cscan_json
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry

damage_state = load_cscan_json("examples/sample_cscan.json")

cscan_cfg = AnalysisConfig(
    material="AS4/3501-6",
    layup_deg=[0, 45, -45, 90] * 4,
    ply_thickness_mm=0.125,
    panel=PanelGeometry(Lx_mm=200, Ly_mm=150),
    damage=damage_state,
    loading="compression",
    tier="empirical",
)
cscan_result = BvidAnalysis(cscan_cfg).run()
print(f"Knockdown: {cscan_result.knockdown:.3f}")
# Knockdown: 0.442
```

### Parametric sweep

`sweep_energies` reruns an impact-driven config (here `impact_cfg` from
above) at each energy and returns a pandas DataFrame:

```python
from bvidfe.sweep.parametric_sweep import sweep_energies

results_df = sweep_energies(
    impact_cfg,
    energies_J=[10, 20, 30, 40, 50],
    csv_path="knockdown_vs_energy.csv",
)
print(results_df[["energy_J", "knockdown", "dpa_mm2"]].round(3))
```

## Physics Models

### Empirical tier

For CAI, the Soutis knockdown (after Soutis & Curtis 1996) scales the pristine strength by `1 / (1 + k_s·(DPA/A_ref)^m)`, where DPA is the projected damage area and `A_ref` = 15,000 mm² is the area of the 150 × 100 mm ASTM D7137 coupon (`SOUTIS_REFERENCE_AREA_MM2`). The reference area is fixed because impact damage is a local notch: it used to be the analysed panel's area, which sent the knockdown to 1 as the panel grew. The material constants `soutis_k_s` and `soutis_m` default to 4.8 and 0.5 for every preset; `k_s` = 4.8 is fitted to the two gated measured CAI datasets (see [Validation](https://github.com/elhajjar1/BVID-FE#validation)). For TAI, the fibre damage is treated as an equivalent circular hole whose diameter is `tai_notch_fraction` (0.05) times the width of the delamination footprint across the load (delaminations cut no fibres, so the hole is far smaller than the footprint). The Whitney-Nuismer point-stress criterion gives its notched strength, with the orthotropic stress concentration `Kt_inf` from Lekhnitskii, a material characteristic distance `wn_d0_mm` (1 mm) and Tan's (1988) orthotropic finite-width correction for the panel width. `tai_notch_fraction` is fitted to one measured TAI condition (in-sample). On 150 × 100 mm IM7/8552 panels the TAI knockdown is 0.82 / 0.77 / 0.68 at 10 / 15 / 30 J on 8-ply QI (CAI 0.45 / 0.40 / 0.31) and 0.84 / 0.80 / 0.73 on [0/90]2s (CAI 0.45 / 0.39 / 0.31). Both models are closed-form.

### Semi-analytical tier

At each delaminated interface the thinner sublaminate over the largest delamination is treated as a plate over the delamination's enclosing rectangle, and a Rayleigh-Ritz closed form gives its buckling load. Buckling alone does not fail the laminate (a one-ply sublaminate over a BVID-sized delamination buckles at a few MPa, and the one at the governing interface typically at 70-370 MPa on IM7/8552 150×100 mm panels at 5-30 J); the CAI strength is the far-field stress at which the buckled sublaminate grows its delamination, from the thin-film energy release rate `G = h/(2E_f)·(σ−σ_c)(σ+3σ_c)` (Chai, Babcock & Knauss 1981; Hutchinson & Suo 1992) set equal to a Benzeggagh-Kenane mixed-mode toughness (`G_Ic`/`G_IIc` at a 52.1° phase angle, η = 2; 0.48 N/mm for IM7/8552). The weakest interface governs, and the result is capped by the Soutis empirical strength. Whitney-Nuismer is retained for TAI. The governing interface is `AnalysisResults.critical_sublaminate` and its buckling load (N/mm) is in `AnalysisResults.buckling_eigenvalues`.

### 3D FE tier

A structured hexahedral mesh is built for the damaged laminate. Delaminated interfaces are approximated by a **component-wise stiffness-reduction model** (true cohesive surfaces deferred to a future release): elements at a delaminated interface carry an out-of-plane factor (`DAMAGE_OOP_FACTOR ≈ 0.05`) that scales the through-thickness and transverse-shear stiffness. The damage zone as a whole (every ply under the projected delamination footprint, where matrix cracks and fibre damage accompany the delaminations) is a soft inclusion with in-plane stiffness ×`DAMAGE_ZONE_INPLANE_FACTOR` (0.15, fitted to three measured CAI coupons; see [Limitations](#limitations); Soutis & Curtis 1996). The delaminations, and so the zone, are kept inside the impacted panel (`impact_panel`, else `panel`; for an ASTM D7136 coupon the 125 × 75 mm support window): each template ellipse is shrunk, aspect ratio kept, until it fits, and the union is re-solved for the same area. The optional fiber-break core under the impact site (`fiber_break_eta`, off in every preset) reduces in-plane stiffness the same way. Elements straddling the zone's edge take the factor weighted by their area fraction inside it. First-ply failure is evaluated on stress rotated into each ply's material frame, using LaRC05 (CAI) and Tsai-Wu (TAI), with the failure index averaged over a 2.5 mm disk within each element layer, each element weighted by the fraction of its area inside the disk (`FPF_AVERAGING_RADIUS_MM`; average-stress criterion, Whitney & Nuismer 1974), so the result converges as the mesh is refined. A damage zone only a few elements across is under-resolved at the default 5 mm mesh (see [Limitations](#limitations)). Damage-zone elements fail at the strain intact material would (strain equivalence; Lemaitre 1992): the criterion sees their effective stress, the in-plane stress the pristine ply would carry at the same strain, so their strength drops with their stiffness. The failure stress is the reaction force on the loaded edge at failure over the gross section, so a softened panel is not credited with load it cannot carry. For CAI, the lower of first-ply failure and the semi-analytical delamination growth stress governs. Buckling onset (whole panel or sublaminate, from the Navier sine-basis Rayleigh-Ritz closed form, issue #129 — the 3D K_g eigensolve previously used here was retired because 3D Hex on thin laminates locks too aggressively at affordable mesh sizes) is reported in `AnalysisResults.buckling_eigenvalues` but does not set the strength. Whole-panel buckling is left out of the knockdown because CAI fixtures carry anti-buckling guides; its stress is recorded in `AnalysisResults.notes`.

### Knockdown definition and cross-tier comparability

`AnalysisResults.knockdown` is computed in exactly one place — `BvidAnalysis.run()` (`src/bvidfe/analysis/bvid.py`):

<!-- docs-test: skip (fragment) -->
```python
knockdown = residual_strength_MPa / pristine_strength_MPa
```

**The denominator is identical across all three tiers.** `_pristine_strength()` (`src/bvidfe/analysis/bvid.py`) is the unnotched laminate strength from classical lamination theory, by ply-discount last-ply failure (`src/bvidfe/failure/laminate_strength.py`): the laminate is loaded uniaxially; a ply whose matrix fails (σ₂ against Yt/Yc, |τ₁₂| against S12) loses its E22 and G12 and sheds its load; the laminate fails when any ply's fibres fail (σ₁ against Xt/Xc) or when every ply has cracked. It matches the NCAMP unnotched tension and compression strengths of AS4/8552 and IM7/8552 laminates (25/50/25, 10/80/10, 50/40/10) within 0.83–1.14 (`tests/validation/test_unnotched_strength_ncamp.py`).

**The numerator (residual strength) is what differs between tiers** (`σ₀` is the pristine strength, `σ_growth` the delamination growth stress of the weakest sublaminate, `σ_FPF` the fe3d first-ply-failure stress):

| Tier | CAI residual stress | TAI residual stress |
|---|---|---|
| `empirical` | Soutis: `σ₀ / (1 + k_s·(DPA/A_ref)^m)`, `A_ref` = 15,000 mm² | Whitney-Nuismer point stress on a hole `tai_notch_fraction` × footprint width across, finite-width corrected |
| `semi_analytical` | `min(Soutis, σ_growth)` — adds the delamination growth stress of the weakest sublaminate | Delegates to the empirical TAI model (mathematically identical to `empirical`) |
| `fe3d` | `σ₀ × min(1, S_damaged / S_undamaged)` with `S = min(σ_growth, σ_FPF,LaRC05)` | `σ₀ × min(1, S_damaged / S_undamaged)` with `S = σ_FPF,Tsai-Wu` |

fe3d's own failure stresses (first-ply failure on the 3D mesh, delamination growth) are not the same kind of quantity as σ₀, so fe3d is normalised by an undamaged fe3d run of the same panel: an undamaged panel gives a knockdown of exactly 1.0 in every tier, and fe3d runs its solve twice. The raw damaged and undamaged fe3d stresses are recorded in `AnalysisResults.notes`.

**What this means for users:**

- All three tiers report knockdown on the **same scale** (ratio relative to the same pristine baseline), so values are *qualitatively comparable*.
- They are **not** numerically interchangeable: each tier captures different failure mechanisms.
  - For **TAI**, `empirical` and `semi_analytical` are mathematically identical; `fe3d` differs and is far more conservative (about one third of the closed-form TAI at 15–30 J; not recommended, see [Limitations](#limitations)).
  - For **CAI**, `semi_analytical ≤ empirical` always (the growth stress only lowers the residual). `fe3d` divides by its own undamaged strength, which is below σ₀ (8-ply QI IM7/8552, 150×100 mm, 5 mm mesh: 599.3 vs 601.3 MPa in compression, 494.5 vs 753.2 MPa in tension; fe3d is first-ply failure, which CLT puts at 601.3 and 502.5 MPa), so the same delamination growth stress gives a slightly higher knockdown in `fe3d` than in `semi_analytical`.
- For **energy-scaling studies**, prefer `empirical` (Soutis scales with DPA). `fe3d` also falls with damage area, through its softened damage zone (16-ply QI tension: 0.282/0.253 at 15/30 J, against empirical 0.81/0.72; fe3d tension is far more conservative than the closed-form TAI, see [Limitations](#limitations)), but its damage-zone stiffness is fitted to only three CAI coupons. The delamination growth stress changes little with delamination size once the sublaminate buckles well below it, so where it governs the knockdown is flatter in energy. It governs `semi_analytical` mainly at low energy and on cross-ply laminates: on a 16-ply [0/45/-45/90]2s IM7/8552 200×150 mm panel at 5 J (0.597 against empirical 0.726; at 10/15/20 J the Soutis bound governs, 0.503/0.427/0.381), on 8-ply QI 150×100 mm panels at 2–5 J (0.500 against 0.546 at 5 J) and on 8-ply cross-ply 150×100 mm panels at 2–15 J (0.424/0.389/0.388 against 0.544/0.447/0.394 at 5/10/15 J). On 16- and 24-ply QI IM7/8552 150×100 mm panels it governs just above the onset energy (16-ply at 5 J: 0.556 against 0.651; 24-ply at 10 J: 0.500 against 0.540) and the Soutis bound from 10 J (16-ply) and 15 J (24-ply).
- A few silent fallbacks affect interpretation:
  - `fe3d` CAI ignores whole-panel buckling (CAI fixtures carry anti-buckling guides). A slender panel without guides can buckle well below the reported residual; the panel buckling stress is recorded in `AnalysisResults.notes`.
  - DPA is globally capped at 80% of panel area (`src/bvidfe/impact/mapping.py`); above this damage threshold all three tiers saturate. The delaminations are also kept inside the impacted panel; on a 150×100 mm IM7/8552 panel that holds the area at about 11,500 mm² from about 110 J on symmetric quasi-isotropic laminates (with a `UserWarning`), before the cap.

### Damage severity heatmap

The "Damage Severity" tab in the Streamlit app is **not** a simple count of damaged interfaces, nor is it a continuous physical damage variable (e.g. a Kachanov-style scalar). It is a through-thickness accumulation of the per-element stiffness-reduction metric used by the `fe3d` mesh:

1. Each hex element in the damaged mesh carries an **out-of-plane** stiffness factor `damage_factor`: `1.0` if pristine, or `DAMAGE_OOP_FACTOR ≈ 0.05` if it is intersected by a delamination interface inside an ellipse footprint, or sits inside the fiber-break core. (Every element in a column under a delamination footprint also carries a reduced **in-plane** factor `in_plane_damage_factor = DAMAGE_ZONE_INPLANE_FACTOR = 0.15`, weighted by the fraction of its area inside the footprint where it straddles the edge, and fiber-break-core elements `DAMAGE_FIBER_BREAK_INPLANE_FACTOR = 0.30`.) The out-of-plane and fiber-break factors are **categorical per element** (geometric overlap tests in `bvidfe.analysis.fe_mesh.build_fe_mesh`), not a continuum damage variable.
2. The heatmap metric is `1 − damage_factor` (the OOP factor) — `0.0` for pristine, `0.95` for fully delaminated.
3. For each in-plane column `(x, y)` the metric is **summed over the through-thickness elements** to produce the heatmap value (see the "Damage Severity" tab in `app.py`).

The colorbar is therefore in units of "stacked OOP-stiffness loss contributions": `0` means no delamination at any interface in that column. A delaminated interface marks the two elements that touch it, so with one element per ply a column under a single delamination reads `1.9` (2 × 0.95), and a column under delaminations at all `n_plies − 1` interfaces reads `0.95 × n_plies` (7.6 for 8 plies). It is a visualization aid analogous to a C-scan operator's depth-projected damage map, not a quantitative continuum-damage field. The in-plane reduction is a separate channel (`mesh.in_plane_damage_factors`) and is not currently overlaid on the heatmap.

## Validation

`validation/validate_bvid_public.py` runs a tier on every case in
`validation/datasets/` and reports the mean absolute error (MAE) between
predicted and measured strength. Three datasets are measured CAI tests and
one a measured TAI test; details and sources are in [`validation/README.md`](https://github.com/elhajjar1/BVID-FE/blob/main/validation/README.md).

| Dataset | Material and layup | Coupons | Measured strength (MPa) | Empirical MAE | Semi-analytical MAE |
|---|---|---|---|---|---|
| `ncamp_as4_8552_cai` (NCAMP CAM-RP-2010-002) | AS4/8552 [45/0/-45/90]3s, about 30 J | 7 | 171–182 | 1.5% (max 3.0%) | 1.5% (max 3.0%) |
| `lovejoy_scotti_im7_8552_cai` (NASA LaRC 2019) | IM7/8552 [45/0/-45/90]4s at 54 J; hard [45/0/90/0/-45/0/45/0/-45/0]2s at 68 J | 8 | 133–147; 236–249 | 8.8% (max 17.4%) | 8.8% (max 17.4%) |
| `hasebe_t800s_3900_2b_cai` (Hasebe et al. 2025), not gated | T800S/3900-2B, 8 layups (8–24 plies), 12.7–31.8 mm hemispheres, a 15.9 mm flat-ended impactor and a 60° cone at 2–21 J; 80×50 mm JIS K 7089 coupons | 343 | 81–592 | 21.5% (max 143.1%) | 23.3% (max 143.1%) |
| `korbelin_m21_t800s_tai` (Körbelin et al. 2020), TAI, not gated | M21/T800S [45/0/-45/90]2s at 25 J, 250×60 mm tension coupons (T800S/3900-2B card) | 5 | 884–1014 (tension) | 3.5% (max 9.8%), in-sample | 3.5% (max 9.8%) |

Read these numbers with care:

- They are in-sample. `soutis_k_s` = 4.8 was fitted to these same two
  datasets, and the default `olsson_alpha` (0.0855) to the median of Lovejoy
  & Scotti's C-scan areas. NCAMP reports no damage area, so its damage area
  is an independent prediction. Fitted on one dataset alone, `k_s` predicts
  the other with an MAE of 8.8% (NCAMP to Lovejoy & Scotti) and 5.0%
  (Lovejoy & Scotti to NCAMP). Expect errors of about 5–9% on a new
  material or layup tested to ASTM D7136/D7137.
- The Soutis bound sets the semi-analytical result on every NCAMP and
  Lovejoy & Scotti coupon (on the hard laminate the delamination growth
  stress is 252 MPa against a 224 MPa Soutis bound), so the two tiers give
  the same strengths there.
- The fe3d tier gives 141.0 MPa on the first NCAMP coupon (measured 177.7
  MPa, 20.6% low), 148.6 MPa on the first Lovejoy & Scotti quasi-isotropic
  coupon (146.8 MPa, 1.2% high) and 236.4 MPa on the first hard-laminate
  coupon (235.7 MPa, 0.3% high). These are in-sample: its damage-zone
  in-plane factor (0.15) is the factor with the lowest mean error on these
  three coupons (7.4%). No single factor fits both laminate families
  within 5% (NCAMP wants about 0.22, Lovejoy & Scotti about 0.15), so
  fe3d strengths are qualitative. On the measured TAI dataset
  (`korbelin_m21_t800s_tai`) fe3d gives a knockdown of 0.335 against the
  measured 0.889 (MAE 62.2%, against 3.5% in-sample for the closed-form
  tiers); fe3d TAI is not recommended (see [Limitations](#limitations)).
  It takes about 40 s per 24-ply coupon (about 80 s for the 40-ply hard
  laminate) and is not run on the measured data in CI.
- **The Hasebe set is a known miss.** The predicted damage area matches
  the C-scans (mean 410 against 364 mm² for the 60° cone, 52 against 52 mm²
  for the flat-ended impactor, 282 against 281 mm² for the 12.7 mm
  hemisphere), but a constant knockdown (the set's mean, 0.73) scores
  30.3%, against 21.5% for the empirical tier. The DPA-based knockdown does
  not capture how much the knockdown depends on layup (cross-ply c16 keeps
  0.64 of its strength, quasi-isotropic q16 0.84, at similar damage area),
  and it is too low on average (0.63 predicted, 0.73 measured). It is
  reported, not gated. Its undamaged coupons fail at 0.4–0.6× the
  unnotched strength (thin coupons in the JIS fixture), so it is scored on
  knockdown times the measured undamaged strength
  (`measured_pristine_MPa`).
- **TAI: one measured condition, `korbelin_m21_t800s_tai` (not gated).**
  Five M21/T800S [45/0/-45/90]2s coupons (Körbelin, Dreiner & Fiedler
  2020, CC BY 4.0) impacted at 25 J keep 0.889 of their tension strength.
  `tai_notch_fraction` is fitted to them, so the 3.5% MAE is in-sample;
  the previous whole-footprint hole scored 57%. There is no M21/T800S
  preset, so the T800S/3900-2B card stands in.
- A fifth dataset, `synthetic_selfcheck`, is generated by the empirical
  tier itself (MAE 0% by construction). It only checks that the harness and
  the empirical pipeline have not changed.

CI fails if the empirical MAE on any gated dataset, or the semi-analytical
MAE on a gated measured dataset, exceeds 1.25× the dataset's 15% target.
The pristine strength is checked separately against NCAMP unnotched tension
and compression strengths (within 0.83–1.14×, see above).

```bash
python validation/validate_bvid_public.py --tier empirical
python validation/validate_bvid_public.py --tier semi_analytical \
    --dataset ncamp_as4_8552_cai lovejoy_scotti_im7_8552_cai hasebe_t800s_3900_2b_cai korbelin_m21_t800s_tai
```

## Limitations

- `soutis_k_s` and the default `olsson_alpha` are calibrated on two CAI datasets (AS4/8552 and IM7/8552, 24- to 40-ply laminates impacted at 30–68 J): `soutis_k_s` on their strengths, `olsson_alpha` on Lovejoy & Scotti's C-scan areas; see [Validation](https://github.com/elhajjar1/BVID-FE#validation). `olsson_alpha` for T800S/3900-2B (0.0605) and the impactor-shape and diameter factors are fitted to Hasebe et al.'s C-scan areas, the only open set that varies the impactor. `dent_beta` and `dent_gamma` are fitted to the measured dents of Lovejoy & Scotti and Hasebe et al. (median predicted / measured 1.01 and 0.94, scatter a factor of about 1.7). The other impact-mapping constants (the impactor-mass factor and related parameters) are uncalibrated defaults. The predicted damage area is within 14% of every Lovejoy & Scotti C-scan (in-sample) and at a median 1.03× of Hasebe et al.'s, but on Hasebe et al.'s coupons the knockdown still misses the effect of layup. Calibrate against coupon data for the material before using results in certification.
- The onset energy uses a single-term large-deflection plate law, `P = k_b w + k_m w³`: one assumed mode per edge condition, the membrane coefficient for a central impact used at every impact point, and edges free to slide in-plane on a simply supported panel (held on a clamped one). It is 0.8 J on an 8-ply, 3.3 J on a 16-ply and 5.7 J on a 24-ply 150 × 100 mm IM7/8552 panel, and 1.2 J on an 8-ply 300 × 200 mm one (29 J with the linear plate). The threshold load is checked against measured threshold loads only on 4–6 mm quasi-isotropic laminates, and no measured onset energy is checked (`docs/physics_models.md`).
- LaRC05 is implemented as a minimal Hashin-3D reduction. Full plane-search fiber-kinking is deferred to a future release.
- The `fe3d` tier uses component-wise stiffness reduction (out-of-plane at delaminated interfaces, in-plane over the damage zone) instead of true cohesive surfaces with bilinear traction-separation laws and a continuum damage model. Cohesive surfaces are deferred to a future release. The damage-zone in-plane factor is uniform over the footprint. Its value, 0.15, is fitted in-sample to three measured CAI coupons, with the damage kept inside the ASTM D7136 impact window: it puts the first Lovejoy & Scotti quasi-isotropic coupon and their first 40-ply hard laminate 1.2% and 0.3% high and the first NCAMP coupon 20.6% low, the lowest mean error (7.4%) of the factors tried (0.13 / 0.17 / 0.20 give 13.8% / 10.5% / 12.4%). NCAMP wants about 0.22 and Lovejoy & Scotti about 0.15, so no single factor fits both families within 5% (with the damage template before it took ply angles modulo 180°, 0.17 did; see `validation/README.md`). Before the footprint was confined the hard laminate's damage band spanned the 102 mm coupon and no single factor fitted both families (0.20 left it 47% low). A damage area at the 80%-of-panel cap would put the soft zone across the loaded edges.
- **`fe3d` TAI is not recommended and not validated.** On 150×100 mm IM7/8552 panels the fe3d tension knockdown is about one third of the closed-form TAI at 15–30 J on quasi-isotropic laminates (16-ply QI 0.282 / 0.253 at 15 / 30 J against 0.813 / 0.718; 24-ply QI 0.246 against 0.731 at 30 J, damage area 2,649 mm²), and on the one measured TAI dataset (`korbelin_m21_t800s_tai`, 25 J) it gives 0.335 against the measured 0.889 (MAE 62.2%). The damage zone's strain-equivalent strength reduction applies in every in-plane direction, so in tension the soft zone fails as if its fibres were cut, while delaminations cut none; the damage-zone factor is fitted on CAI coupons only. Use the closed-form tiers for TAI.
- **The `fe3d` first-ply failure length scale is fixed, not fitted.** The failure index is averaged over a 2.5 mm disk (`FPF_AVERAGING_RADIUS_MM`, within the 1-4 mm characteristic distances of CFRP notched strength), each element weighted by the fraction of its area inside the disk; apart from that length, the stress concentration at the edge of the soft damage zone depends on the zone's shape and on the panel width, not its absolute size, so fe3d's notch-size effect is weaker than that of the empirical tier's Soutis and Whitney-Nuismer models. At the default 5 mm mesh the disk holds only the element itself. Cohesive surfaces and proper load-introduction BCs are left to a future release.
- **At the default 5 mm mesh a small `fe3d` damage zone is under-resolved.** On 8-ply IM7/8552 150×100 mm panels at 15 J (damage area about 1,533 mm² QI, 1,538 mm² cross-ply) the knockdown at 5 / 2.5 / 2.0 mm elements is 0.256 / 0.274 / 0.268 for QI compression, 0.242 / 0.273 / 0.269 for QI tension, 0.167 / 0.168 / 0.163 for cross-ply compression and 0.156 / 0.168 / 0.162 for cross-ply tension: 2.5 to 2.0 mm within about 2–4%, and the 5 mm value within about 3–10% of the 2.0 mm one (mostly below it). Smaller zones, only a few elements across at 5 mm, are less well resolved: with earlier code, at damage areas of about 710–855 mm², the 5 mm cross-ply compression knockdown was about 25% above its converged value. The calibration coupons' damage areas (2,800–6,500 mm²) are well resolved at 5 mm. Use 2.5 mm or finer when the damage area is below about 1,000 mm² (about 35 s for an 8-ply 150×100 mm panel); `FE3D_MAX_DOF` (150,000) admits such a panel down to about 1.67 mm.
- **The delamination growth stress is a 1D thin-film estimate.** It uses the straight-sided blister energy release rate over the delamination's enclosing rectangle, a fixed mixed-mode phase angle (52.1°), and treats the rest of the laminate as rigid. It is least accurate when the weakest sublaminate is near the mid-plane, which is common: the minimum over interfaces often lands on a 3-5 ply sublaminate. Once a sublaminate buckles at less than about 0.3× its steady-state growth stress `√(2E_f·G_c/h)`, the growth stress stays between 0.87× and 1× that value, so larger delaminations at the same interface lower it little and can raise it slightly. Knockdowns from this channel therefore follow impact energy mainly when a different interface becomes the weakest. In `fe3d` the first-ply failure of the softened damage zone falls with energy as well: on IM7/8552 150×100 mm panels (5 mm mesh) the fe3d CAI knockdown at 5 / 15 / 30 J is 0.430 / 0.302 / 0.263 on 16-ply QI (empirical 0.651 / 0.417 / 0.322), 0.333 / 0.256 / 0.225 on 8-ply QI and 0.307 / 0.167 / 0.134 on an 8-ply cross-ply panel, every energy above the onset.
- Validation covers CAI: two carbon/8552 systems (gated) and T800S/3900-2B (reported, not gated). **TAI is not validated**: `tai_notch_fraction` is fitted in-sample to one measured condition (M21/T800S quasi-isotropic coupons at 25 J, run with the T800S/3900-2B card), so its energy, layup and thickness trends are untested, and a fixed fraction of the footprint cannot follow extra fibre failure (from a sharper impactor, impact at elevated temperature or energies near penetration). The AS4/3501-6, T700/2510 and T800/epoxy presets have no measured residual-strength comparison. Candidate datasets are listed in [`validation/README.md`](https://github.com/elhajjar1/BVID-FE/blob/main/validation/README.md).

## Citation

If you use BVID-FE in your research, please cite:

```bibtex
@software{elhajjar2026bvidfe,
  author    = {Elhajjar, Rani},
  title     = {{BVID-FE}: Barely Visible Impact Damage residual-strength analysis
               for composite laminates},
  year      = {2026},
  version   = {0.2.0},
  publisher = {GitHub},
  url       = {https://github.com/elhajjar1/BVID-FE},
  note      = {University of Wisconsin-Milwaukee}
}
```

## References

- Olsson, R. (2001). Analytical prediction of large mass impact damage in composite laminates. *Composites Part A*, 32(9), 1207-1215.
- Soutis, C. (1996). Compressive strength of unidirectional composites: measurement and prediction. *ASTM STP*, 1242, 168-176.
- Whitney, J.M. & Nuismer, R.J. (1974). Stress fracture criteria for laminated composites containing stress concentrations. *Journal of Composite Materials*, 8(3), 253-265.
- Tan, S.C. (1988). Finite-width correction factors for anisotropic plate containing a central opening. *Journal of Composite Materials*, 22(11), 1080-1097.
- Körbelin, J., Dreiner, C. & Fiedler, B. (2020). Impact of temperature on LVI-damage and tensile and compressive residual strength of CFRP. *Composites Part C: Open Access*, 3, 100074.
- Tsai, S.W. & Wu, E.M. (1971). A general theory of strength for anisotropic materials. *Journal of Composite Materials*, 5(1), 58-80.
- Davila, C.G., Camanho, P.P., & Rose, C.A. (2005). Failure criteria for FRP laminates. NASA/TM-2005-213530 (LaRC05).
- Soutis, C. & Curtis, P.T. (1996). Prediction of the post-impact compressive strength of CFRP laminated composites. *Composites Science and Technology*, 56(6), 677-684.
- Chai, H., Babcock, C.D., & Knauss, W.G. (1981). One dimensional modelling of failure in laminated plates by delamination buckling. *International Journal of Solids and Structures*, 17(11), 1069-1083.
- Hutchinson, J.W. & Suo, Z. (1992). Mixed mode cracking in layered materials. *Advances in Applied Mechanics*, 29, 63-191.
- Benzeggagh, M.L. & Kenane, M. (1996). Measurement of mixed-mode delamination fracture toughness of unidirectional glass/epoxy composites with mixed-mode bending apparatus. *Composites Science and Technology*, 56(4), 439-449.
- Lemaitre, J. (1992). *A Course on Damage Mechanics*. Springer.
- Daniel & Ishai (2006). *Engineering Mechanics of Composite Materials*, sec. 7.6 (ply-discount laminate strength).
- NCAMP CAM-RP-2010-002 Rev A (2011). Hexcel 8552 AS4 unitape qualification, sections 2.3.31 and 4.31.
- Lovejoy & Scotti (2019). Potential Weight Benefits of IM7/8552 Hybrid Thin-Ply Composites for Aircraft Structures. NASA Langley Research Center, NTRS 20200002432.
- Hasebe, Higuchi, Yokozeki & Takeda (2025). Dataset of compression after impact testing on carbon fiber reinforced plastic laminates. *Data in Brief*, 60, 111509. Per-specimen data: Mendeley Data, doi:10.17632/8scdmfdcfb.3 (CC BY 4.0).

## Contributing

See [CONTRIBUTING.md](https://github.com/elhajjar1/BVID-FE/blob/main/CONTRIBUTING.md) for guidelines on reporting bugs, suggesting features, and submitting pull requests.

## License

MIT License. See [LICENSE](https://github.com/elhajjar1/BVID-FE/blob/main/LICENSE) for details.
