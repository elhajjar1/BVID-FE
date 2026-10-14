---
name: add-material
description: Add a new composite material preset (e.g. AS4/3501-6, IM7/8552, T700/2510, T800/epoxy) to the BVID-FE material library. Use when the user wants to add a material, register a new ply, support a new prepreg, or extend MATERIAL_LIBRARY.
---

# Add a new material to MATERIAL_LIBRARY

Materials live in `src/bvidfe/core/material.py`. The library is a dict
named `MATERIAL_LIBRARY` (NOT `MATERIAL_PRESETS` — older docs use that
name, but the actual symbol is `MATERIAL_LIBRARY`; `_resolve_material`
in `analysis/bvid.py` and `analysis/fe_tier.py` looks up by string key).
There are six presets today: AS4/3501-6, AS4/8552, IM7/8552, T700/2510,
T800S/3900-2B, T800/epoxy.

## Pre-flight

Before adding the entry, the user must provide a **literature source**
for every constant. BVID-FE refuses to ship a material without one. If
they don't have a source, stop and ask — vendor datasheets, published
papers, or NASA/CMH-17 handbooks all qualify.

Confirm with the user:

1. Material name (the dict key, e.g. `"AS4/3501-6"`)
2. The source(s) — citations go in a comment block directly above the
   entry, one per constant or group of constants, as the `AS4/8552` entry
   does (test report and table, conditioning, and how any derived value
   was obtained, e.g. `G23 = E22 / (2 (1 + nu23))` with the cited `nu23`).
   A derived or borrowed constant (such as a toughness taken from another
   material with the same matrix) is allowed only if the comment says so
   and cites the source of the assumption.

## Required constants and units

All constants are required and must be `> 0`. Units:

| Field | Unit | Notes |
|---|---|---|
| `E11`, `E22`, `G12`, `G13`, `G23` | MPa | Moduli |
| `nu12` | dimensionless | Must satisfy `-1.0 < nu12 < 0.5` |
| `Xt`, `Xc`, `Yt`, `Yc`, `S12`, `S23` | MPa | Strengths (`S23` is the transverse 2-3 shear strength) |
| `G_Ic`, `G_IIc` | N/mm | Interlaminar fracture toughness |
| `rho` | **kg/mm³** | e.g. `1.58e-6` for CFRP. NOT g/cm³ — `__post_init__` rejects values outside `[1e-7, 1e-5]` because the old wrong docstring caused silent 1000x errors |

Optional (leave as `None` for transverse-isotropic fallback):

- `Zt`, `Zc` — through-thickness tensile/compressive strengths (fall
  back to `Yt`/`Yc`)
- `S13` — 1-3 plane shear (falls back to `S12`)

Optional calibration fields (`olsson_alpha`, `dent_beta`, `dent_gamma`,
`fiber_break_eta`, `fiber_break_E_threshold`, `soutis_k_s`, `soutis_m`,
`wn_d0_mm`, `puck_p_nt_plus`, `puck_p_nt_minus`) keep their defaults
unless the user has material-specific calibration data. The
`soutis_k_s` default (4.7) is calibrated on the measured CAI datasets in
`validation/datasets/`, and the `olsson_alpha` default (0.078) on Lovejoy &
Scotti's IM7/8552 C-scan areas; do not override either per material
without measured data for that material (CAI strengths for `soutis_k_s`,
C-scan damage areas for `olsson_alpha`, as for T800S/3900-2B), and check
the fitted value on data that was not used to fit it (see
`validation/README.md`).

## Steps

### 1. Add the entry to `MATERIAL_LIBRARY`

Append the new entry to the dict in `src/bvidfe/core/material.py`,
matching the formatting of the existing entries (one field per line,
alphabetical-ish but really just match the existing pattern). The
`name` field on the dataclass must equal the dict key.

### 2. Add tests

`test_preset_name_matches_key_and_compliance_inverts_stiffness` in
`tests/core/test_material.py` is parametrized over every key in
`MATERIAL_LIBRARY`, so the new entry automatically gets the
`name == key`, matrix-shape and `S @ C ≈ I` checks, and importing the
module runs `__post_init__` validation. In addition:

- Add the name to the preset list in
  `test_material_library_has_six_presets` (`tests/core/test_material.py`)
  and in `test_cli_list_materials_prints_all_presets`
  (`tests/test_cli_version.py`).
- Add a test that pins the cited values, mirroring
  `test_as4_8552_lamina_properties_are_ncamp_rtd_values` and
  `test_as4_8552_derived_constants`, so a later edit cannot silently drift
  from the source.

Mirror the structure of the existing tests in the file — don't invent
a new pattern.

### 3. Validate and format

```bash
ruff check src tests
black src tests
pytest tests/core/test_material.py tests/test_cli_version.py -v
```

### 4. Mention downstream surfaces (do NOT edit unless the user asks)

The new key automatically flows through:

- `AnalysisConfig.material` (`Union[str, OrthotropicMaterial]`)
- The Streamlit material picker in `app.py` (`MATERIAL_NAMES`)
- The CLI `--material` choices and `bvidfe --list-materials`

The JSON Schema does not list material names (`material` is a free
string), so adding a preset needs no schema regeneration. Adding or
changing an `OrthotropicMaterial` *field* does: use the
`regenerate-schemas` skill.

If `AnalysisConfig` ever gains a typed `Literal` of material names (grep
`MATERIAL_LIBRARY` to be sure), extend it too. As of this writing it has
none; materials are looked up by string at runtime.

Prose that counts or lists the presets needs a manual update: `README.md`
("Six material presets"), `docs/index.md`, and the `core/material.py`
row in `ARCHITECTURE.md`. Add an entry under `## [Unreleased]` in
`CHANGELOG.md`.

Nothing requires a validation dataset for a new preset. If the source
also publishes CAI/TAI coupon data or unnotched laminate strengths for the
material, offer to add them (a dataset under `validation/datasets/`, see
`validation/README.md`, or a check like
`tests/validation/test_unnotched_strength_ncamp.py`).

## Things this skill must never do

- Never invent material constants. If the user can't cite a source,
  stop and ask.
- Never use g/cm³ for `rho` — kg/mm³ only (`1.58e-6`, not `1.58`).
- Never skip the test addition. CI runs `pytest` and an unguarded
  entry will silently break the next time someone tweaks the dataclass.
