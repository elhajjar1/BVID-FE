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
print(f"{result.knockdown:.3f}")   # 0.593
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
# AnalysisResults(tier='empirical', kd=0.593, sigma=356.4/601.3MPa, dpa=307.5mm2, notes=0, warnings=0)
```

Internally `BvidAnalysis.run()` dispatches `impact_to_damage()` (Olsson
threshold + energy-balance damage area + peanut-template distribution +
dent model) to build a `DamageState`, then runs the selected tier's
residual-strength engine. The predicted damage area is capped at 80% of
the panel area, with a `DPACapClipWarning`, and the delaminations are kept
inside the impacted panel. On this panel no damage is predicted below
the onset energy (7.2 J), the area then grows by about 110 mm² per joule;
from about 109 J the footprint fills what the panel can hold (about
11,200 mm², with a `UserWarning`), before the 12,000 mm² cap.

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
# AnalysisResults(tier='semi_analytical', kd=0.442, sigma=265.5/601.3MPa, dpa=1040.8mm2, notes=0, warnings=0)
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

## Configuration reference

`AnalysisConfig` fields (lengths in mm, angles in degrees):

| Field | Default | Meaning |
|---|---|---|
| `material` | required | Preset name from `MATERIAL_LIBRARY` or an `OrthotropicMaterial` |
| `layup_deg` | required | Ply angles in stacking order; `0` is the load direction (x) |
| `ply_thickness_mm` | required | One thickness for every ply, or a list with one per ply (`len(layup_deg)` entries) |
| `panel` | required | `PanelGeometry(Lx_mm, Ly_mm, boundary="simply_supported")`, the tested panel; x is the load direction |
| `loading` | `"compression"` | `"compression"` (CAI) or `"tension"` (TAI) |
| `tier` | `"empirical"` | `"empirical"`, `"semi_analytical"` or `"fe3d"` |
| `impact` | `None` | `ImpactEvent(energy_J, impactor=ImpactorGeometry(), mass_kg=5.5, location_xy_mm=(0.0, 0.0))`; `(0, 0)` means the panel centre |
| `damage` | `None` | A `DamageState` in the frame of `panel` (exactly one of `impact` and `damage`) |
| `mesh` | `None` | `MeshParams` for fe3d; filled with `MeshParams()` when `tier="fe3d"`, ignored by the other tiers |
| `impact_panel` | `None` | The panel the impact was applied on, when it is not `panel` (see above); its `boundary` is the one the impact model uses |

`MeshParams` fields (fe3d only):

| Field | Default | Meaning |
|---|---|---|
| `elements_per_ply` | `1` | Hex elements through each ply (a positive `int`) |
| `in_plane_size_mm` | `5.0` | Target element edge; the mesh has `ceil(Lx_mm / size)` × `ceil(Ly_mm / size)` equal elements |
| `cohesive_zone_factor` | `1.0` | Reserved for cohesive surfaces; validated (> 0) but not used |

The first-ply failure index is averaged over a 2.5 mm disk
(`FPF_AVERAGING_RADIUS_MM` in `bvidfe.analysis.fe_tier`), so at the
default 5 mm mesh it is checked element by element; finer meshes average
over neighbouring elements. Check the size of a mesh before running it:

```python
from dataclasses import replace
from bvidfe import MeshParams
from bvidfe.analysis.fe_mesh import estimate_fe_mesh_size

fine = replace(config, tier="fe3d", mesh=MeshParams(in_plane_size_mm=2.5))
print(estimate_fe_mesh_size(fine)["n_dof"])   # 67527 (fe3d refuses more than 150,000)
```

## Results

`AnalysisResults` fields:

| Field | Meaning |
|---|---|
| `residual_strength_MPa` | Damaged strength from the selected tier |
| `pristine_strength_MPa` | Unnotched laminate strength (same for every tier) |
| `knockdown` | `residual_strength_MPa / pristine_strength_MPa` |
| `damage`, `dpa_mm2` | The `DamageState` used and its projected damage area (mm²) |
| `tier_used`, `config_snapshot` | Provenance |
| `buckling_eigenvalues` | Buckling *onset*, not strength, as a one-element list: sublaminate `N_cr` (N/mm) for `semi_analytical` CAI, lower of panel and sublaminate buckling stress (MPa) for `fe3d` CAI; `None` for `empirical`, for TAI and for `semi_analytical` without delaminations |
| `critical_sublaminate` | `semi_analytical` and `fe3d` CAI: interface whose delamination grows at the lowest stress; else `None` |
| `notes` | Human-readable diagnostics; fe3d records its raw damaged / undamaged stresses here (and, for CAI, the delamination-growth and whole-panel buckling stresses). Impact-mapping regimes are Python warnings, not notes |
| `warnings` | Every `UserWarning` raised during the run, as `"<WarningClass>: <message>"` (also emitted as Python warnings) |
| `field_results` | Always `None` for now |

`result.summary()` returns a printable block and `result.to_dict()` a
JSON-serialisable dict (every field but `field_results`, the damage in the
C-scan layout; the CLI prints it).

## Tier choice

| tier | Runtime | Validation status | Good for |
|---|---|---|---|
| `empirical` | Tens of milliseconds | CAI calibrated (`soutis_k_s` and the default `olsson_alpha`) on two measured datasets: MAE 1.6% (NCAMP AS4/8552) and 8.7% (Lovejoy & Scotti IM7/8552); expect about 5–9% on a new material, and 21% on Hasebe et al.'s T800S/3900-2B coupons (not gated; a constant knockdown scores 30%). TAI not validated. | Design allowables, energy sweeps, quick screening. Soutis scales with DPA. |
| `semi_analytical` | Tens of milliseconds | Same datasets: MAE 1.6% and 8.7% (23% on Hasebe et al.'s coupons). Identical to `empirical` for TAI. | Checking whether delamination growth governs before the Soutis bound; never less conservative than `empirical` for CAI. |
| `fe3d` | Seconds to a minute (two FE solves; about 5 s for the 8-ply 150×100 mm panel at the default 5 mm mesh, about 35 s at 2.5 mm, longer for thick laminates) | **Qualitative** — the damage-zone stiffness (×0.15) is fitted in-sample to three measured CAI coupons, quasi-isotropic and hard (Lovejoy & Scotti within 3%, NCAMP 19% low; no single factor fits both datasets within 5%); not run on the measured datasets in CI | Stress-field context, damage-through-thickness. Knockdown falls with damage area but has no notch-size effect. For damage areas below about 1,500 mm² use `MeshParams(in_plane_size_mm=2.5)` or finer: at 5 mm such a zone is under-resolved (knockdown up to about 25% high). |

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
  footprints, every ply) is softened in-plane ×0.15 with strain-equivalent
  strengths, an element partly inside the zone softened by its area
  fraction. The failure index is averaged over a 2.5 mm disk in each layer,
  each element weighted by its overlap with the disk, before it is checked, and the stress is the reaction force at failure
  over the gross section. Buckling onset (Rayleigh-Ritz closed form, issue #129) is
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
take `csv_path=` to write the CSV directly (no index column), `on_error=`
(`"raise"` (default), `"skip"` or `"warn"`; a failed point keeps its input,
NaN numbers, `n_delaminations=0`, `tier_used=""` and the message in
`error`, which is NaN on success; a point whose config fails validation,
such as a layup that does not match a per-ply thickness list, fails the
same way), and `progress_callback=` (called as
`cb(i, n)` or `cb(i, n, inputs)`, where `inputs` is e.g.
`{"energy_J": 10.0}`). `sweep_thicknesses` sets a uniform thickness, so it
replaces a per-ply list. Each impact warning category is shown once per
sweep.

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
print(f"{tsai_wu_index(m, s):.3f} {larc05_index(m, s):.3f} {puck_index(m, s):.3f}")
# 0.705 0.253 0.503
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

The tiers' residual-strength functions take the laminate, the damage and
the pristine strength (MPa):

```python
from bvidfe.failure.soutis_openhole import SOUTIS_REFERENCE_AREA_MM2, soutis_cai
from bvidfe.analysis.semi_analytical import semi_analytical_cai

damage = load_cscan_json("examples/sample_cscan.json")
sigma_0 = unnotched_strength(lam, "compression")
dpa = damage.projected_damage_area_mm2
print(SOUTIS_REFERENCE_AREA_MM2)                  # 15000.0 (the 150 x 100 mm coupon)
print(f"{soutis_cai(lam.material, dpa, SOUTIS_REFERENCE_AREA_MM2, sigma_0):.0f}")   # 266 MPa
r = semi_analytical_cai(lam, damage, sigma_0, boundary="clamped")
print(f"{r.residual_strength_MPa:.0f}", r.critical_interface_index)   # 266 4
```

`soutis_cai(m, dpa_mm2, A_ref_mm2, sigma_pristine_MPa)` is the empirical
CAI strength `sigma_pristine_MPa / (1 + m.soutis_k_s * (dpa_mm2 /
A_ref_mm2) ** m.soutis_m)`; `BvidAnalysis` passes
`SOUTIS_REFERENCE_AREA_MM2` (15,000 mm², fixed, not the panel area).
`semi_analytical_cai(lam, damage, sigma_pristine_MPa,
A_ref_mm2=SOUTIS_REFERENCE_AREA_MM2, boundary="simply_supported")` returns
a `SemiAnalyticalResult` (`residual_strength_MPa`,
`critical_interface_index`, `critical_buckling_load_N` in N/mm), the
lower of that Soutis bound and the delamination growth stress; `boundary`
scales the sublaminate buckling load (clamped ×1.9, free ×0.5).
`whitney_nuismer_tai(m, dpa_mm2, sigma_pristine_MPa, Kt_inf)` (same
module) is the TAI counterpart; pass `Kt_inf=lekhnitskii_kt_infinity(lam)`.

`impact_to_damage(event, lam, panel)` (`bvidfe.impact.mapping`) is the
impact model on its own, and `save_cscan_json(damage, path)`
(`bvidfe.damage.io`) writes any `DamageState`, for example
`result.damage`, as C-scan JSON for an inspection-driven rerun.

## CLI

`bvidfe` (or `python -m bvidfe.cli`) runs one analysis per tier and prints
the result:

```
bvidfe --material IM7/8552 --layup 0,45,-45,90,90,-45,45,0 --thickness 0.152 \
       --panel 150x100 --loading compression --energy 10 --quick
# 0.593 (printed to six decimals)
```

| Flag | Default | Meaning |
|---|---|---|
| `--material NAME` | required | Preset: `AS4/3501-6`, `AS4/8552`, `IM7/8552`, `T700/2510`, `T800/epoxy`, `T800S/3900-2B` |
| `--layup A,B,...` | required | Comma-separated ply angles in degrees |
| `--thickness T` | required | Ply thickness in mm, or a comma-separated list with one per ply |
| `--panel 150x100` | required | Panel size in mm, e.g. `150x100` (simply supported; the CLI has no boundary option) |
| `--loading` | required | `compression` or `tension` |
| `--tier` | `empirical` | One tier, or a comma-separated list such as `empirical,semi_analytical,fe3d` |
| `--energy J` | — | Impact energy in J (impact-driven); exactly one of `--energy` and `--cscan` |
| `--cscan PATH` | — | C-scan JSON file (inspection-driven) |
| `--impactor-diameter D` | `16.0` | Impactor diameter in mm (hemispherical; the CLI has no shape option) |
| `--mass M` | `5.5` | Impactor mass in kg |
| `--quick` | off | Print only the knockdown (6 decimals) |
| `--quick-json` | off | Print `{knockdown, residual_strength_MPa, pristine_strength_MPa, tier_used}` as one-line JSON |
| `--list-materials` | — | List the presets (E11, E22, Xt, Xc, Yt, Yc in MPa) and exit |
| `--version`, `-h/--help` | — | Print the version or the help and exit |

Without `--quick` / `--quick-json` the full result is printed as JSON
(`AnalysisResults.to_dict()`). With several tiers the output is a JSON
array, one JSON line per tier (`--quick-json`), or `tier<TAB>knockdown`
lines (`--quick`). `--quick` and `--quick-json` cannot be combined.

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
    olsson_alpha=0.0865,  # fraction of the energy above onset that delaminates (DPA)
    soutis_k_s=4.8,       # Soutis CAI constant (calibrated on measured CAI data)
    soutis_m=0.5,         # Soutis exponent
    wn_d0_mm=1.0,         # Whitney-Nuismer characteristic distance
    dent_beta=0.044,      # Dent depth d = dent_beta (E - E_onset)^dent_gamma, mm (E in J)
    dent_gamma=0.57,      # Dent-depth exponent (both fitted to measured dents)
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
default to `Yt`, `Yc`, `S12`), `nu23` (transverse Poisson ratio of the fe3d
3D stiffness; defaults to `E22 / (2 G23) - 1`, 0.41 above),
`puck_p_nt_plus` / `puck_p_nt_minus` (0.30 /
0.25), and `fiber_break_eta` / `fiber_break_E_threshold` (fibre-break core,
off by default). The defaults of `soutis_k_s` and `olsson_alpha` are
calibrated against measured CAI data on the 8552-matrix presets (AS4/8552
and IM7/8552), and `olsson_alpha` for T800S/3900-2B (0.0648) against Hasebe
et al.'s C-scan areas; treat the others as defaults.

## Logging

The fe3d tier logs stage-by-stage timing lines under the `bvidfe.fe3d`
logger: mesh build, element build, the first-ply-failure solve and its
result, and the closed-form buckling stress. The empirical and
semi-analytical tiers do not log. The `bvidfe` CLI and the Streamlit app
print these lines to stderr; as a library, bvidfe attaches only a
`NullHandler`, so from Python call
`bvidfe.analysis.fe_tier.enable_stderr_logging()` (or configure `logging`
yourself) to see them.

```bash
BVIDFE_LOG_LEVEL=WARNING bvidfe ...    # DEBUG, INFO (default), WARNING, ERROR
```

## Environment variables

- `BVIDFE_FE3D_MAX_DOF` — hard cap on fe3d problem size (default 150,000).
  Above it `BvidAnalysis.run()` raises `FE3DSizeError`
  (`bvidfe.analysis.fe_tier`). Raise at your own risk: a run needs about
  0.5 GB at 67k DOFs and ~1 GB at 150k, but the direct-solver fallback
  (taken only if the iterative solve fails) needs several times more, and
  an out-of-memory error inside scipy's native code kills the process. A
  value that is not a positive integer falls back to the default with a
  warning on stderr. It is read once, at import.
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
