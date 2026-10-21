# BVID-FE Python API quick reference

For users who want to script BVID-FE workflows beyond what the Streamlit
app and CLI provide. Covers the public types, typical idioms, and where to
find more.

## Top-level entry points

```python
from bvidfe.analysis import (
    AnalysisConfig,      # config dataclass — material, layup, panel, tier, ...
    BvidAnalysis,        # orchestrator — .run() executes a single analysis
    AnalysisResults,     # result dataclass — knockdown, damage, eigs, notes, ...
    MeshParams,          # fe3d-only mesh resolution (elements_per_ply, in_plane_size_mm)
    FieldResults,        # reserved for fe3d field output; not populated yet (always None)
)
```

The most used types are also re-exported from the package root:
`from bvidfe import AnalysisConfig, AnalysisResults, BvidAnalysis,
MeshParams, ImpactEvent, ImpactorGeometry, PanelGeometry,
MATERIAL_LIBRARY, OrthotropicMaterial, DamageState, DelaminationEllipse,
load_cscan_json`.

The core idiom is:

```python
from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, PanelGeometry

config = AnalysisConfig(material="IM7/8552", layup_deg=[0,45,-45,90,90,-45,45,0],
                        ply_thickness_mm=0.152, panel=PanelGeometry(150, 100),
                        loading="compression", tier="empirical",
                        impact=ImpactEvent(energy_J=10))
result = BvidAnalysis(config).run()
print(result.knockdown)   # 0.4609...
```

## Workflow paths

### Impact-driven (the most common)

```python
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry

config = AnalysisConfig(
    material="IM7/8552",
    layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ply_thickness_mm=0.152,
    panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
    impact=ImpactEvent(
        energy_J=10.0,
        impactor=ImpactorGeometry(diameter_mm=16),
        mass_kg=5.5,
    ),
    loading="compression",
    tier="empirical",
)
result = BvidAnalysis(config).run()
print(result)
# AnalysisResults(tier='empirical', kd=0.461, sigma=277.2/601.3MPa, dpa=928.3mm2, notes=0, warnings=0)
```

Internally `BvidAnalysis.run()` dispatches `impact_to_damage()` (Olsson
threshold + energy-balance damage area + peanut-template distribution +
dent model) to build a `DamageState`, then runs the selected tier's
residual-strength engine. The predicted damage area is capped at 80% of
the panel area, with a `DPACapClipWarning`, and the delaminations are kept
inside the impacted panel. On this panel the area grows by about 99 mm²
per joule above onset; from about 114 J the footprint fills what the
panel can hold (about 11,200 mm², with a `UserWarning`), before the
12,000 mm² cap.

`ImpactorGeometry` defaults to the ASTM D7136 impactor (`diameter_mm=16.0`,
`shape="hemispherical"`; also `"flat"` and `"conical"`), and `ImpactEvent`
to `mass_kg=5.5` at the panel centre. `PanelGeometry` takes
`boundary="simply_supported"` (default), `"clamped"` or `"free"`. State the
impactor explicitly when you know it: diameter, mass and shape all change
the predicted damage area.

### Inspection-driven (you have a C-scan)

```python
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry
from bvidfe.damage.io import load_cscan_json

damage = load_cscan_json("examples/sample_cscan.json")
cscan_config = AnalysisConfig(
    material="IM7/8552",
    layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ply_thickness_mm=0.152,
    panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
    damage=damage,            # note: provide damage OR impact, not both
    loading="compression",
    tier="semi_analytical",
)
cscan_result = BvidAnalysis(cscan_config).run()
print(cscan_result)
# AnalysisResults(tier='semi_analytical', kd=0.447, sigma=268.7/601.3MPa, dpa=1040.8mm2, notes=0, warnings=0)
print(cscan_result.critical_sublaminate)   # 4: the interface whose delamination grows first
```

`DamageState` and `DelaminationEllipse` are the universal handoff types
— see the [C-scan schema](cscan_schema.md) for the JSON format.
`AnalysisConfig` raises `ValueError` unless exactly one of `impact` and
`damage` is set, and for an unknown `tier` or `loading`.
`ply_thickness_mm` may also be a list with one thickness per ply.
`impact_panel` (a `PanelGeometry`, with `impact` only) is the panel the
impact was applied on when it is not the one tested, e.g. a 60 × 60 mm
clamped window and an 80 × 50 mm compression coupon: the damage is mapped
on `impact_panel` and centred on `panel`. For an ASTM D7136 coupon it is
the 125 × 75 mm support window (`PanelGeometry(125, 75)`). It also bounds
the damage footprint: every delamination ellipse is kept inside it,
measured from the impact point. Without `impact_panel` the tested panel
bounds it.

## Results

`AnalysisResults` fields:

| Field | Meaning |
|---|---|
| `residual_strength_MPa` | Damaged strength from the selected tier |
| `pristine_strength_MPa` | Unnotched laminate strength (same for every tier) |
| `knockdown` | `residual_strength_MPa / pristine_strength_MPa` |
| `damage`, `dpa_mm2` | The `DamageState` used and its projected damage area (mm²) |
| `tier_used`, `config_snapshot` | Provenance |
| `buckling_eigenvalues` | Buckling *onset*, not strength: sublaminate `N_cr` (N/mm) for `semi_analytical` CAI, lower of panel and sublaminate buckling stress (MPa) for `fe3d` CAI, else `None` |
| `critical_sublaminate` | CAI only: interface whose delamination grows at the lowest stress |
| `notes` | Human-readable diagnostics; fe3d records its raw damaged / undamaged stresses here |
| `warnings` | Machine-readable tags; no backend emits one at present |
| `field_results` | Always `None` for now |

`result.summary()` returns a printable block and `result.to_dict()` a
JSON-serialisable dict.

## Tier choice

| tier | Runtime | Validation status | Good for |
|---|---|---|---|
| `empirical` | < 1 s | CAI calibrated (`soutis_k_s` and the default `olsson_alpha`) on two measured datasets: MAE 1.5% (NCAMP AS4/8552) and 9.0% (Lovejoy & Scotti IM7/8552); expect about 6–9% on a new material, and 21% on Hasebe et al.'s T800S/3900-2B coupons (not gated; a constant knockdown scores 30%). TAI not validated. | Design allowables, energy sweeps, quick screening. Soutis scales with DPA. |
| `semi_analytical` | < 1 s | Same datasets: MAE 1.5% and 9.0% (23% on Hasebe et al.'s coupons). Identical to `empirical` for TAI. | Checking whether delamination growth governs before the Soutis bound; never less conservative than `empirical` for CAI. |
| `fe3d` | Tens of seconds (two FE solves; about 30 s for the 8-ply 150×100 mm panel at the default 5 mm mesh) | **Qualitative** — the damage-zone stiffness (×0.17) is fitted to three measured CAI coupons, quasi-isotropic and hard (within 5%); not run on the measured datasets in CI | Stress-field context, damage-through-thickness. Knockdown falls with damage area but has no notch-size effect. |

## Knockdown semantics

`AnalysisResults.knockdown == residual_strength_MPa / pristine_strength_MPa`.
The pristine strength is the unnotched laminate strength (ply-discount
last-ply failure on classical lamination theory) and is **identical for all
three tiers** — only the residual-strength numerator differs:

- `empirical`: Soutis (CAI) / Whitney-Nuismer (TAI) closed-form.
- `semi_analytical`: `min(Soutis, delamination growth)` for CAI, where
  delamination growth is the far-field stress at which the weakest buckled
  sublaminate grows its delamination; delegates to Whitney-Nuismer for TAI
  (so its TAI knockdown is mathematically identical to `empirical`).
- `fe3d`: `min(delamination growth, first-ply-failure)` for CAI and
  first-ply failure for TAI, each divided by the same quantity on the
  undamaged panel (capped at 1) and applied to the pristine reference.
  First-ply failure runs on a mesh whose damage zone (the delamination
  footprints, every ply) is softened in-plane ×0.17 with strain-equivalent
  strengths, and its stress is the reaction force at failure over the gross
  section. Buckling onset (Rayleigh-Ritz closed form, issue #129) is
  reported in `result.buckling_eigenvalues` but does not set the strength;
  whole-panel buckling is left out because CAI fixtures carry anti-buckling
  guides, and its stress is recorded in `result.notes`.

Cross-tier expectations:

- For **CAI**, `semi_analytical.knockdown ≤ empirical.knockdown` always
  (the growth stress only lowers the residual). `fe3d` is independent
  and can sit above or below `empirical`.
- For **TAI**, `empirical.knockdown == semi_analytical.knockdown` exactly;
  `fe3d` differs.
- All three are on the same scale and qualitatively comparable, but
  numerical agreement is not guaranteed.

See [Physics Models → Knockdown definition and cross-tier
comparability](physics_models.md#knockdown-definition-and-cross-tier-comparability)
for the full breakdown and the limitations.

## Parametric sweeps

```python
from bvidfe.sweep.parametric_sweep import (
    sweep_energies, sweep_layups, sweep_thicknesses,
)

# Each returns a pandas DataFrame with columns
# [energy_J, knockdown, residual_MPa, pristine_MPa, dpa_mm2, dent_mm,
#  n_delaminations, tier_used, error]
df = sweep_energies(config, energies_J=[5, 10, 15, 20, 25, 30])  # the impact-driven config above
df.to_csv("sweep.csv")
```

`sweep_energies` needs a config with `impact` set. `sweep_layups(config,
layups=[...])` and `sweep_thicknesses(config, ply_thicknesses_mm=[...])`
put `layup` or `ply_thickness_mm` in the first column instead. All three
take `csv_path=` to write the CSV directly, `on_error=` (`"raise"`
(default), `"skip"` or `"warn"`; a failed point becomes a row of NaNs with
the message in `error`), and `progress_callback=` (called as
`cb(i, n)` or `cb(i, n, inputs)`).

## Visualization

```python
from bvidfe.viz.plots_2d import (
    plot_damage_map,           # ellipse footprints + panel outline
    plot_knockdown_curve,      # single tier
    plot_tier_comparison,      # overlaid multi-tier
)

fig = plot_damage_map(result.damage, config.panel)
fig.savefig("damage_map.png", dpi=300)
```

`plot_knockdown_curve(energies_J, knockdowns, tier_label="")` and
`plot_tier_comparison(energies_J, {"empirical": [...], ...})` take plain
sequences, for example two columns of a sweep DataFrame. All three return
a matplotlib Figure.

For an interactive 3D view of the damaged hex mesh (the same figure the
Streamlit app's "3D Damage" tab shows), build the mesh and render it with
Plotly:

```python
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.viz.plotly_3d import mesh_damage_figure  # FeMesh -> plotly Figure

mesh = build_fe_mesh(config, result.damage)
mesh_damage_figure(mesh).write_html("damage_3d.html")
```

The colour is `1 − damage_factor`, where `damage_factor` is each element's
out-of-plane stiffness factor (0.05 at delaminated interfaces). The
in-plane damage-zone factor is in `mesh.in_plane_damage_factors`;
`stress_field_figure(mesh, element_scalar, label=...)` renders any
per-element array.

## Failure criteria and building blocks

The ply failure indices take a material and a material-frame Voigt stress
`[σ1, σ2, σ3, τ23, τ13, τ12]` (MPa); an index ≥ 1 means failure:

```python
from bvidfe import MATERIAL_LIBRARY
from bvidfe.failure.tsai_wu import tsai_wu_index
from bvidfe.failure.larc05 import larc05_index
from bvidfe.failure.puck import puck_index

m = MATERIAL_LIBRARY["IM7/8552"]
s = [-800.0, 20.0, 0.0, 0.0, 0.0, 30.0]
print(tsai_wu_index(m, s), larc05_index(m, s), puck_index(m, s))
# 0.7227... 0.2531... 0.5031...
```

`FailureEvaluator(material, criterion)` (`bvidfe.failure.evaluator`;
`criterion` is `"tsai_wu"`, `"larc05"` or `"puck"`) applies one to an
`(n_elem, n_gp, 6)` stress array and reports the maximum index and where
it occurs.

The closed-form pieces behind the tiers can be called directly:

```python
from bvidfe import MATERIAL_LIBRARY, load_cscan_json
from bvidfe.core.laminate import Laminate
from bvidfe.failure.laminate_strength import unnotched_strength
from bvidfe.analysis.semi_analytical import weakest_sublaminate_growth

lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
print(unnotched_strength(lam, "compression"))    # 601.28... MPa (pristine strength)
print(lam.effective_engineering_constants())     # (Ex, Ey, Gxy, nu_xy): 62434 MPa, 62434 MPa, 23885 MPa, 0.307

g = weakest_sublaminate_growth(lam, load_cscan_json("examples/sample_cscan.json"))
print(g.interface_index, g.growth_stress_MPa)    # 4 275.88... MPa
```

`weakest_sublaminate_growth` returns a `SublaminateGrowth` with the
interface, its buckling load (N/mm), the film buckling and growth stresses
and the far-field growth stress (MPa); see
[Physics Models](physics_models.md#semi-analytical-tier).

## CLI

Available flags:
```
bvidfe --help
bvidfe --version
bvidfe --list-materials
bvidfe --material NAME --layup 0,45,-45,90 --thickness 0.152 \
       --panel 150x100 --loading compression --tier empirical \
       --energy 30                                # impact-driven
bvidfe ... --energy 30 --impactor-diameter 16 --mass 5.5   # impactor (defaults shown)
bvidfe ... --cscan path.json                      # inspection-driven
bvidfe ... --tier empirical,semi_analytical,fe3d  # several tiers in one call
bvidfe ... --quick                                # knockdown scalar only
bvidfe ... --quick-json                           # knockdown, strengths and tier as one-line JSON
```

`--thickness` also takes a comma-separated list with one thickness per ply.
Without `--quick` / `--quick-json` the full result is printed as JSON.

## Adding a custom material

Presets: `AS4/3501-6`, `AS4/8552`, `IM7/8552`, `T700/2510`, `T800S/3900-2B`, `T800/epoxy`
(`MATERIAL_LIBRARY`, or `bvidfe --list-materials`).

```python
from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, PanelGeometry
from bvidfe.core.material import OrthotropicMaterial

my_mat = OrthotropicMaterial(
    name="my_cfrp",
    E11=150e3, E22=9e3, nu12=0.32,
    G12=4.5e3, G13=4.5e3, G23=3.2e3,
    Xt=2200, Xc=1500, Yt=65, Yc=180, S12=85, S23=55,
    G_Ic=0.27, G_IIc=0.85, rho=1.58e-6,  # rho in kg/mm^3 (1.58 g/cm^3)
    # Optional impact/calibration constants (defaults shown):
    olsson_alpha=0.078,   # fraction of the energy above onset that delaminates (DPA)
    soutis_k_s=4.7,       # Soutis CAI constant (calibrated on measured CAI data)
    soutis_m=0.5,         # Soutis exponent
    wn_d0_mm=1.0,         # Whitney-Nuismer characteristic distance
    dent_beta=0.05,       # Dent-depth prefactor
    dent_gamma=0.5,       # Dent-depth exponent
)

config = AnalysisConfig(
    material=my_mat,  # pass the object directly instead of a preset name
    layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ply_thickness_mm=0.152,
    panel=PanelGeometry(150, 100),
    impact=ImpactEvent(energy_J=10),
)
print(BvidAnalysis(config).run().knockdown)
```

Moduli and strengths are in MPa, toughnesses `G_Ic` / `G_IIc` in N/mm.
Other optional fields: `Zt`, `Zc`, `S13` (through-thickness strengths;
default to `Yt`, `Yc`, `S12`), `puck_p_nt_plus` / `puck_p_nt_minus` (0.30 /
0.25), and `fiber_break_eta` / `fiber_break_E_threshold` (fibre-break core,
off by default). The defaults of `soutis_k_s` and `olsson_alpha` are
calibrated against measured CAI data on the 8552-matrix presets (AS4/8552
and IM7/8552), and `olsson_alpha` for T800S/3900-2B (0.047) against Hasebe
et al.'s C-scan areas; treat the others as defaults.

## Logging

The fe3d tier streams stage-by-stage timing lines to stderr under the
`bvidfe.fe3d` logger: mesh build, element build, the first-ply-failure
solve and its result, and the closed-form buckling stress. The empirical
and semi-analytical tiers do not log.

```bash
BVIDFE_LOG_LEVEL=WARNING bvidfe ...    # DEBUG, INFO (default), WARNING, ERROR
```

## Environment variables

- `BVIDFE_FE3D_MAX_DOF` — hard cap on fe3d problem size (default 500,000).
  Above it `BvidAnalysis.run()` raises `FE3DSizeError`
  (`bvidfe.analysis.fe_tier`). Raise at your own risk; past this size
  scipy's sparse solvers can exhaust memory and SIGSEGV.
- `BVIDFE_LOG_LEVEL` — level of the `bvidfe.fe3d` logger; see above. An
  invalid value falls back to `INFO` with a warning on stderr.

## See also

- `examples/01_empirical_quick.py` — minimal end-to-end
- `examples/02_tier_comparison.py` — all three tiers on one damage
- `examples/03_energy_sweep.py` — sweep + plot
- `examples/04_inspection_driven.py` — C-scan JSON loading
- `examples/05_tier_comparison_sweep.py` — empirical vs. semi_analytical over energies
- [C-scan schema](cscan_schema.md) — C-scan JSON format spec
- [Physics Models](physics_models.md) — equations, defaults and limitations
- `validation/README.md` — measured datasets and the `soutis_k_s` and `olsson_alpha` calibration
- `ARCHITECTURE.md` — module dependency graph and data flow
- `CHANGELOG.md` — what changed in each release, with measured effects
