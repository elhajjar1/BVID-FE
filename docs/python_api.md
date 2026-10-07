# BVID-FE Python API quick reference

For users who want to script BVID-FE workflows beyond what the GUI and CLI
provide. Covers the public types, typical idioms, and where to find more.

## Top-level entry points

```python
from bvidfe.analysis import (
    AnalysisConfig,      # config dataclass — material, layup, panel, tier, ...
    BvidAnalysis,        # orchestrator — .run() executes a single analysis
    AnalysisResults,     # result dataclass — knockdown, damage, eigs, ...
    MeshParams,          # fe3d-only mesh resolution (elements_per_ply, in_plane_size_mm)
    FieldResults,        # fe3d-only per-element stress/strain/failure fields
)
```

The core idiom is:

```python
config = AnalysisConfig(material="IM7/8552", layup_deg=[0,45,-45,90,90,-45,45,0],
                        ply_thickness_mm=0.152, panel=PanelGeometry(150, 100),
                        loading="compression", tier="empirical",
                        impact=ImpactEvent(energy_J=20, ...))
result = BvidAnalysis(config).run()
print(result.knockdown)
```

## Workflow paths

### Impact-driven (the most common)

```python
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry

config = AnalysisConfig(
    material="IM7/8552",
    layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ply_thickness_mm=0.152,
    panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
    impact=ImpactEvent(
        energy_J=30.0,
        impactor=ImpactorGeometry(diameter_mm=16),
        mass_kg=5.5,
    ),
    loading="compression",
    tier="empirical",
)
result = BvidAnalysis(config).run()
```

Internally `BvidAnalysis.run()` dispatches `impact_to_damage()` (Olsson
threshold + peanut-template DPA distribution + dent model) to build a
`DamageState`, then runs the selected tier's residual-strength engine.

### Inspection-driven (you have a C-scan)

```python
from bvidfe.damage.io import load_cscan_json

damage = load_cscan_json("path/to/measurement.json")
config = AnalysisConfig(
    material="IM7/8552",
    layup_deg=[...],
    ply_thickness_mm=0.152,
    panel=PanelGeometry(Lx_mm=150, Ly_mm=100),
    damage=damage,            # note: provide damage OR impact, not both
    loading="compression",
    tier="semi_analytical",
)
result = BvidAnalysis(config).run()
```

`DamageState` and `DelaminationEllipse` are the universal handoff types
— see `docs/cscan_schema.md` for the JSON format.

## Tier choice

| tier | Runtime | Knockdown trustworthy? | Good for |
|---|---|---|---|
| `empirical` | < 1 s | Yes — closed-form Soutis (CAI) / Whitney-Nuismer (TAI) | Design allowables, energy sweeps, quick screening. Soutis scales with DPA. |
| `semi_analytical` | ~ 1 s | Yes — lower of Soutis and delamination growth for CAI; identical to `empirical` for TAI | More conservative where delamination growth governs; flatter in energy than `empirical`. |
| `fe3d` | ~ 10 s | **Qualitative** — the damage-zone stiffness (×0.30) is uncalibrated; see note below | Stress-field context, damage-through-thickness. Knockdown falls with damage area but has no notch-size effect. |

## Knockdown semantics

`AnalysisResults.knockdown == residual_strength_MPa / pristine_strength_MPa`.
The pristine strength is a thickness-weighted ply-average of the lamina-level
material strengths and is **identical for all three tiers** — only the
residual-strength numerator differs:

- `empirical`: Soutis (CAI) / Whitney-Nuismer (TAI) closed-form.
- `semi_analytical`: `min(Soutis, delamination growth)` for CAI, where
  delamination growth is the far-field stress at which the weakest buckled
  sublaminate grows its delamination; delegates to Whitney-Nuismer for TAI
  (so its TAI knockdown is mathematically identical to `empirical`).
- `fe3d`: `min(delamination growth, first-ply-failure)` for CAI and
  first-ply failure for TAI, each divided by the same quantity on the
  undamaged panel and applied to the pristine reference. First-ply failure
  runs on a mesh whose damage zone (the delamination footprints, every ply)
  is softened in-plane ×0.3 with strain-equivalent strengths, and its
  stress is the reaction force at failure over the gross section. Buckling onset
  (Rayleigh-Ritz closed form, issue #129) is reported in
  `result.buckling_eigenvalues` but does not set the strength; whole-panel
  buckling is left out because CAI fixtures carry anti-buckling guides,
  and its stress is recorded in `result.notes`.

Cross-tier expectations:

- For **CAI**, `semi_analytical.knockdown ≤ empirical.knockdown` always
  (the growth stress only lowers the residual). `fe3d` is independent
  and can sit above or below `empirical`.
- For **TAI**, `empirical.knockdown == semi_analytical.knockdown` exactly;
  `fe3d` differs.
- All three are on the same scale and qualitatively comparable, but
  numerical agreement is not guaranteed.

See README → "Knockdown definition and cross-tier comparability" for the
full breakdown and limitation cross-links.

## Parametric sweeps

```python
from bvidfe.sweep.parametric_sweep import (
    sweep_energies, sweep_layups, sweep_thicknesses,
)

# Each returns a pandas DataFrame with columns
# [energy_J, knockdown, residual_MPa, pristine_MPa, dpa_mm2, dent_mm, ...]
df = sweep_energies(config, energies_J=[5, 10, 15, 20, 25, 30])
df.to_csv("sweep.csv")
```

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

For an interactive 3D view of the damaged hex mesh (the same figure the
Streamlit app's "3D Damage" tab shows), build the mesh and render it with
Plotly:

```python
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.viz.plotly_3d import mesh_damage_figure  # FeMesh -> plotly Figure

mesh = build_fe_mesh(config, result.damage)
mesh_damage_figure(mesh).write_html("damage_3d.html")
```

## CLI

Available flags:
```
bvidfe --help
bvidfe --version
bvidfe --list-materials
bvidfe --material NAME --layup 0,45,-45,90 --thickness 0.152 \
       --panel 150x100 --loading compression --tier empirical \
       --energy 30                                # impact-driven
bvidfe ... --cscan path.json                      # inspection-driven
bvidfe ... --quick                                # knockdown scalar only
```

## Adding a custom material

```python
from bvidfe.core.material import OrthotropicMaterial

my_mat = OrthotropicMaterial(
    name="my_cfrp",
    E11=150e3, E22=9e3, nu12=0.32,
    G12=4.5e3, G13=4.5e3, G23=3.2e3,
    Xt=2200, Xc=1500, Yt=65, Yc=180, S12=85, S23=55,
    G_Ic=0.27, G_IIc=0.85, rho=1.58e-6,  # rho in kg/mm^3 (1.58 g/cm^3)
    # Optional impact/calibration constants (defaults for CFRP):
    olsson_alpha=0.8,     # DPA scaling
    soutis_k_s=2.5,       # Soutis CAI constant
    soutis_m=0.5,         # Soutis exponent
    wn_d0_mm=1.0,         # Whitney-Nuismer characteristic distance
    dent_beta=0.05,       # Dent-depth prefactor
    dent_gamma=0.5,       # Dent-depth exponent
)

config = AnalysisConfig(material=my_mat, ...)  # pass the object directly
```

## Logging

Both CLI and GUI stream stage-by-stage timing lines to stderr under:
- `bvidfe.fe3d` — mesh build, K assembly, FPF analytic solve; buckling
  channel is a closed-form delegation (no log line of its own)
- `bvidfe.gui` — worker heartbeats, tab update timings

```bash
BVIDFE_LOG_LEVEL=DEBUG bvidfe ...    # DEBUG, INFO (default), WARNING, ERROR
```

## Environment variables

- `BVIDFE_FE3D_MAX_DOF` — hard cap on fe3d problem size (default 500,000).
  Raise at your own risk; past this size scipy's sparse solvers can
  exhaust memory and SIGSEGV.
- `BVIDFE_LOG_LEVEL` — see above.
- `QT_QPA_PLATFORM=offscreen` — headless testing of the GUI.

## See also

- `examples/01_empirical_quick.py` — minimal end-to-end
- `examples/02_tier_comparison.py` — all three tiers on one damage
- `examples/03_energy_sweep.py` — sweep + plot
- `examples/04_inspection_driven.py` — C-scan JSON loading
- `examples/05_tier_comparison_sweep.py` — empirical vs. semi_analytical over energies
- `docs/cscan_schema.md` — C-scan JSON format spec
- `ARCHITECTURE.md` — module dependency graph and data flow
- `CHANGELOG.md` — tier limitations and v0.2.0-dev known issues
