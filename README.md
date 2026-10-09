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

Low-velocity impacts on composite structures can create internal delaminations that are invisible to the naked eye yet significantly degrade compression and tension strength. Engineers certifying composite airframes and pressure vessels need fast, reliable estimates of how much strength is lost and how the knockdown depends on layup, panel size, and impact energy. BVID-FE provides three modeling tiers — from a 30-millisecond empirical lookup to a full 3D finite element solution — so the right level of fidelity can be chosen for each stage of the design process.

BVID-FE is the third in a family of defect-specific composite tools, joining **PorosityFE** (porosity defects) and **WrinkleFE** (fiber waviness). The three tools share material models, laminate theory, failure criteria, and documentation conventions.

## Features

- **Two workflow paths** converging on a shared `DamageState`:
  - *Impact-driven*: Olsson quasi-static threshold + peanut-template DPA distribution + empirical dent model
  - *Inspection-driven*: C-scan JSON import per the documented schema (`docs/cscan_schema.md`)
- **Three modeling tiers** for residual strength after BVID:
  - *Empirical*: Soutis CAI knockdown + Whitney-Nuismer TAI (seconds)
  - *Semi-analytical*: buckling-driven delamination growth (Rayleigh-Ritz sublaminate buckling + thin-film energy release rate), capped by Soutis; Whitney-Nuismer for TAI (seconds)
  - *3D FE*: First-ply-failure on a damaged hexahedral mesh; LaRC05 for CAI, Tsai-Wu for TAI (minutes)
- **CAI and TAI loading modes** (Compression-After-Impact and Tension-After-Impact)
- **Per-interface ellipse damage model** using `DelaminationEllipse` with shapely-union projected damage area
- **Five material presets**: AS4/3501-6, AS4/8552, IM7/8552, T700/2510, T800/epoxy
- **CLI** for single-run and batch use
- **Streamlit web app** with sidebar-driven `AnalysisConfig` and result tabs — Summary, Damage Map, Knockdown Curve (live empirical sweep around the current energy), 3D Damage (Plotly hex mesh coloured by damage factor), Buckling (tier-specific buckling indicator: `N_cr` in N/mm for `semi_analytical`, `σ_crit` in MPa for `fe3d`), Damage Severity (through-thickness sum of per-element `1 − damage_factor`; see [Physics Models](#damage-severity-heatmap)), and Sweep (CSV-exportable parametric energy sweep). Cached runs, fe3d mesh-size guards, JSON/CSV downloads
- **Parametric sweeps** over impact energy, layup, or ply thickness with CSV output
- **2D plots**: damage map, knockdown curves, tier comparison charts
- **Plotly 3D viz** for the damaged hex mesh (`bvidfe.viz.plotly_3d`)

## Installation

```bash
git clone https://github.com/elhajjar1/BVID-FE.git
cd bvid-fe
pip install -e ".[dev]"
pytest -v
```

## Quick Start

### Try it in Jupyter

The fastest tour of the physics is [`examples/quickstart.ipynb`](examples/quickstart.ipynb) — a runnable notebook that walks through the Olsson damage threshold, peanut-template DPA distribution, Soutis CAI, and Whitney-Nuismer TAI with LaTeX derivations and inline Plotly damage map / knockdown sweep / 3D damage mesh visualisations. Open it with `jupyter notebook examples/quickstart.ipynb` (after `pip install -e . jupyter`).

### Streamlit web app

```bash
streamlit run app.py
```

Opens the BVID-FE web UI in your browser: configure material + layup, panel size, impact event or C-scan inspection, pick a tier, click **Run analysis**. Results appear across the Summary / Damage Map / 3D Damage / Knockdown / Buckling / Damage Severity / Sweep tabs. Use the JSON / CSV download buttons to export results. See [DEPLOYMENT_STREAMLIT.md](DEPLOYMENT_STREAMLIT.md) for deploying a public URL on Streamlit Community Cloud.

### Command-line interface

```bash
bvidfe --material IM7/8552 \
       --layup "0,45,-45,90,90,-45,45,0" \
       --thickness 0.152 \
       --panel 150x100 \
       --loading compression \
       --energy 30
```

All lengths are in **millimetres**. `--panel` is `<Lx>x<Ly>` (in-plane
dimensions, x then y, literal `x` separator, no spaces — e.g. `150x100`
is Lx = 150 mm, Ly = 100 mm). `--thickness` is the per-ply thickness in
mm, `--layup` is a comma-separated list of ply angles in degrees, and
`--energy` is the impact energy in joules.

### Python API — impact-driven path

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
        energy_J=30,
        impactor=ImpactorGeometry(diameter_mm=16),
        mass_kg=5.5,
    ),
    loading="compression",
    tier="semi_analytical",
)
result = BvidAnalysis(cfg).run()
print(result.summary())
# -> residual CAI strength, knockdown, DPA, dent depth, per-interface delaminations
```

### Inspection-driven path (C-scan import)

```python
from bvidfe.damage.io import load_cscan_json
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry

damage_state = load_cscan_json("scan.json")

cfg = AnalysisConfig(
    material="AS4/3501-6",
    layup_deg=[0, 45, -45, 90] * 4,
    ply_thickness_mm=0.125,
    panel=PanelGeometry(Lx_mm=200, Ly_mm=150),
    damage=damage_state,
    loading="compression",
    tier="empirical",
)
result = BvidAnalysis(cfg).run()
print(f"Knockdown: {result.knockdown:.3f}")
```

### Parametric sweep

```python
from bvidfe.sweep.parametric_sweep import sweep_energies

results_df = sweep_energies(
    cfg,
    energies_J=[10, 20, 30, 40, 50],
    csv_path="knockdown_vs_energy.csv",
)
```

## Physics Models

### Empirical tier

Soutis open-hole-equivalent CAI model relates the delamination-affected stress concentration around the projected damage area to a net-section failure. Whitney-Nuismer point-stress and average-stress criteria are used for TAI. Both models are closed-form and run in milliseconds.

### Semi-analytical tier

At each delaminated interface the thinner sublaminate over the largest delamination is treated as a plate over the delamination's enclosing rectangle, and a Rayleigh-Ritz closed form gives its buckling load. Buckling alone does not fail the laminate (a thin sublaminate over a BVID-sized delamination buckles at a few MPa); the CAI strength is the far-field stress at which the buckled sublaminate grows its delamination, from the thin-film energy release rate `G = h/(2E_f)·(σ−σ_c)(σ+3σ_c)` (Chai, Babcock & Knauss 1981; Hutchinson & Suo 1992) set equal to a Benzeggagh-Kenane mixed-mode toughness (`G_Ic`/`G_IIc` at a 52.1° phase angle, η = 2; 0.48 N/mm for IM7/8552). The weakest interface governs, and the result is capped by the Soutis empirical strength. Whitney-Nuismer is retained for TAI. The governing interface is `AnalysisResults.critical_sublaminate` and its buckling load (N/mm) is in `AnalysisResults.buckling_eigenvalues`.

### 3D FE tier

A structured hexahedral mesh is built for the damaged laminate. Delaminated interfaces are approximated by a **component-wise stiffness-reduction model** (true cohesive surfaces deferred to a future release): elements at a delaminated interface carry an out-of-plane factor (`DAMAGE_OOP_FACTOR ≈ 0.05`) that scales the through-thickness and transverse-shear stiffness. The damage zone as a whole (every ply under the projected delamination footprint, where matrix cracks and fibre damage accompany the delaminations) is a soft inclusion with in-plane stiffness ×`DAMAGE_ZONE_INPLANE_FACTOR` (0.30, a calibration value; Soutis & Curtis 1996). The optional fiber-break core under the impact site (`fiber_break_eta`, off in every preset) reduces in-plane stiffness the same way. First-ply-failure is evaluated at all Gauss points, on stress rotated into each ply's material frame, using LaRC05 (CAI) and Tsai-Wu (TAI). Damage-zone elements fail at the strain intact material would (strain equivalence; Lemaitre 1992): the criterion sees their effective stress, the in-plane stress the pristine ply would carry at the same strain, so their strength drops with their stiffness. The failure stress is the reaction force on the loaded edge at failure over the gross section, so a softened panel is not credited with load it cannot carry. For CAI, the lower of first-ply failure and the semi-analytical delamination growth stress governs. Buckling onset (whole panel or sublaminate, from the Navier sine-basis Rayleigh-Ritz closed form, issue #129 — the 3D K_g eigensolve previously used here was retired because 3D Hex on thin laminates locks too aggressively at affordable mesh sizes) is reported in `AnalysisResults.buckling_eigenvalues` but does not set the strength. Whole-panel buckling is left out of the knockdown because CAI fixtures carry anti-buckling guides; its stress is recorded in `AnalysisResults.notes`.

### Knockdown definition and cross-tier comparability

`AnalysisResults.knockdown` is computed in exactly one place — `BvidAnalysis.run()` (`src/bvidfe/analysis/bvid.py`):

```python
knockdown = residual_strength_MPa / pristine_strength_MPa
```

**The denominator is identical across all three tiers.** `_pristine_strength()` (`src/bvidfe/analysis/bvid.py`) is a thickness-weighted ply-average of the lamina-level strengths from the material card:

- CAI: `Σ tᵢ (Xc·cos²θᵢ + Yc·sin²θᵢ) / Σ tᵢ`
- TAI: `Σ tᵢ (Xt·cos²θᵢ + Yt·sin²θᵢ) / Σ tᵢ`

**The numerator (residual strength) is what differs between tiers:**

| Tier | CAI residual stress | TAI residual stress |
|---|---|---|
| `empirical` | Soutis: `σ₀ / (1 + k_s·(DPA/A_panel)^m)` | Whitney-Nuismer point-stress on equivalent hole |
| `semi_analytical` | `min(Soutis, σ_growth)` — adds the delamination growth stress of the weakest sublaminate | Delegates to Whitney-Nuismer (mathematically identical to `empirical`) |
| `fe3d` | `σ₀ × min(σ_growth, FPF_LaRC05)` damaged ÷ the same quantity undamaged | `σ₀ × FPF_Tsai-Wu` damaged ÷ undamaged |

fe3d's own failure stresses (first-ply failure on the 3D mesh, delamination growth) are not the same kind of quantity as σ₀, so fe3d is normalised by an undamaged fe3d run of the same panel: an undamaged panel gives a knockdown of exactly 1.0 in every tier, and fe3d runs its solve twice. The raw damaged and undamaged fe3d stresses are recorded in `AnalysisResults.notes`.

**What this means for users:**

- All three tiers report knockdown on the **same scale** (ratio relative to the same pristine baseline), so values are *qualitatively comparable*.
- They are **not** numerically interchangeable: each tier captures different failure mechanisms.
  - For **TAI**, `empirical` and `semi_analytical` are mathematically identical; `fe3d` differs.
  - For **CAI**, `semi_analytical ≤ empirical` always (the growth stress only lowers the residual). `fe3d` divides by its own undamaged strength (first-ply failure, typically below σ₀), so for the same growth stress its knockdown is higher than `semi_analytical`'s — see the flat-vs-energy caveat in [Limitations](#limitations).
- For **energy-scaling studies**, prefer `empirical` (Soutis scales with DPA). `fe3d` also falls with damage area, through its softened damage zone (16-ply QI tension: 0.46/0.33/0.20 at 5/15/30 J), but its damage-zone stiffness is uncalibrated. The delamination growth stress in `semi_analytical` changes little with delamination size once the sublaminate buckles well below it, so that tier is flatter in energy.
- A few silent fallbacks affect interpretation:
  - `fe3d` CAI ignores whole-panel buckling (CAI fixtures carry anti-buckling guides). A slender panel without guides can buckle well below the reported residual; the panel buckling stress is recorded in `AnalysisResults.notes`.
  - DPA is globally capped at 80% of panel area (`src/bvidfe/impact/mapping.py`); above this damage threshold all three tiers saturate.

### Damage severity heatmap

The "Damage Severity" tab in the Streamlit app is **not** a simple count of damaged interfaces, nor is it a continuous physical damage variable (e.g. a Kachanov-style scalar). It is a through-thickness accumulation of the per-element stiffness-reduction metric used by the `fe3d` mesh:

1. Each hex element in the damaged mesh carries an **out-of-plane** stiffness factor `damage_factor`: `1.0` if pristine, or `DAMAGE_OOP_FACTOR ≈ 0.05` if it is intersected by a delamination interface inside an ellipse footprint, or sits inside the fiber-break core. (Every element in a column under a delamination footprint also carries a reduced **in-plane** factor `in_plane_damage_factor = DAMAGE_ZONE_INPLANE_FACTOR ≈ 0.30`, as do fiber-break-core elements.) The factors are **categorical per element** (geometric overlap tests in `bvidfe.analysis.fe_mesh.build_fe_mesh`), not a continuum damage variable.
2. The heatmap metric is `1 − damage_factor` (the OOP factor) — `0.0` for pristine, `0.95` for fully delaminated.
3. For each in-plane column `(x, y)` the metric is **summed over the through-thickness elements** to produce the heatmap value (see the "Damage Severity" tab in `app.py`).

The colorbar is therefore in units of "stacked OOP-stiffness loss contributions": `0` means no delamination at any interface in that column, and the maximum (`≈ n_plies − 1` for a single fully-delaminated interface; higher if multiple interfaces are delaminated and `elements_per_ply > 1`) means every element through the thickness sits inside a damaged region. It is a visualization aid analogous to a C-scan operator's depth-projected damage map, not a quantitative continuum-damage field. The in-plane reduction is a separate channel (`mesh.in_plane_damage_factors`) and is not currently overlaid on the heatmap.

## Limitations

- Material calibration constants (`olsson_alpha`, `soutis_k_s`, `dent_beta`, and related parameters) are reasonable defaults for typical CFRP systems. Precise values need to be calibrated against material-specific coupon test data before use in certification.
- LaRC05 is implemented as a minimal Hashin-3D reduction. Full plane-search fiber-kinking is deferred to a future release.
- The `fe3d` tier uses component-wise stiffness reduction (out-of-plane at delaminated interfaces, in-plane over the damage zone) instead of true cohesive surfaces with bilinear traction-separation laws and a continuum damage model. Cohesive surfaces are deferred to a future release. The damage-zone in-plane factor (0.30) is uncalibrated against test data and uniform over the footprint. Once the damage area hits its 80%-of-panel cap the zone reaches the loaded edges, and fe3d knockdowns fall below the empirical tier's (16- and 24-ply QI tension at 30 J: 0.20 and 0.16 vs 0.35).
- **The `fe3d` first-ply failure has no characteristic length.** The stress concentration at the edge of the soft damage zone depends on the zone's shape and on the panel width, not on the zone's absolute size, so fe3d lacks the notch-size effect of the empirical tier's Soutis and Whitney-Nuismer models. It also depends on the mesh near the zone edge: 1-8% between 5 mm and 2.5 mm elements at 15 J on 8-ply panels. Cohesive surfaces and proper load-introduction BCs are v0.3.0 scope.
- **The delamination growth stress is a 1D thin-film estimate.** It uses the straight-sided blister energy release rate over the delamination's enclosing rectangle, a fixed mixed-mode phase angle (52.1°), and treats the rest of the laminate as rigid. It is least accurate when the weakest sublaminate is near the mid-plane, which is common: the minimum over interfaces often lands on a 3-5 ply sublaminate. Once a sublaminate buckles at less than about 0.3× its steady-state growth stress `√(2E_f·G_c/h)`, the growth stress stays between 0.87× and 1× that value, so larger delaminations at the same interface lower it little and can raise it slightly. Knockdowns from this channel therefore follow impact energy mainly when a different interface becomes the weakest: on 16-ply QI IM7/8552 (150×100 mm) the fe3d knockdown falls from 0.53 at 5 J to 0.36 at 30 J, while on 8-ply QI and cross-ply panels it is flat or rises slightly (0.46 → 0.48, 0.39 → 0.42).
- No validated datasets included; comparison against published Soutis, Caprino, and NASA datasets is on the roadmap.

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
- Tsai, S.W. & Wu, E.M. (1971). A general theory of strength for anisotropic materials. *Journal of Composite Materials*, 5(1), 58-80.
- Davila, C.G., Camanho, P.P., & Rose, C.A. (2005). Failure criteria for FRP laminates. NASA/TM-2005-213530 (LaRC05).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on reporting bugs, suggesting features, and submitting pull requests.

## License

MIT License. See [LICENSE](LICENSE) for details.
