# Animations

BVID-FE renders 1920 × 1080 videos of the damage it predicts, for talks,
reports and teaching. Three scenes are available:

- **Energy sweep**: the delaminations of every interface growing as the
  impact energy rises, beside the damage area, the dent depth and the
  knockdown against energy.
- **Damage stack**: the delaminations through the thickness, as an exploded
  3D stack of plies with the camera orbiting it.
- **fe3d load ramp**: the damaged panel loaded to its residual strength in
  the fe3d model, each element coloured by how close it is to failing,
  until the first ply fails.

Each is rendered by the `bvidfe-animate` command, from Python, or in the
Streamlit app's **Animate** tab. The
default output is an MP4 (H.264, `yuv420p`, 30 fps, CRF 18, `+faststart`),
which plays in PowerPoint, browsers, Teams and GitHub. A PNG sequence and a
reduced-size GIF preview are also available.

Every frame carries the citation of BVID-FE and its version in the footer.
If you show or publish an animation, please also cite the software (see
[Citation](#citation)).

## Install

MP4 output needs an ffmpeg executable. The `anim` extra installs one
(through `imageio-ffmpeg`, a bundled ffmpeg binary):

```bash
pip install "bvidfe[anim] @ git+https://github.com/elhajjar1/BVID-FE.git"
```

BVID-FE looks for ffmpeg in this order:

1. the `BVIDFE_FFMPEG` environment variable (an error if it names no
   executable);
2. matplotlib's `rcParams["animation.ffmpeg_path"]`, when it is set to a
   path that exists;
3. `ffmpeg` on `PATH`;
4. the binary of `imageio-ffmpeg` (the `anim` extra).

Without ffmpeg, `--format mp4` stops with an error that names the extra,
and `--format auto` writes a PNG sequence instead, with a warning. GIF and
PNG output need no ffmpeg.

## Energy sweep

![Energy sweep preview: the delaminations grow with the impact energy, beside the damage area, dent depth and knockdown curves](animations/energy_sweep.gif)

*Preview (640 × 360, 12 fps). The 1080p MP4 is about 1.8 MB for 10 s.*

```bash
bvidfe-animate energy-sweep --material IM7/8552 \
    --layup 0,45,-45,90,0,45,-45,90,0,45,-45,90,0,45,-45,90 \
    --thickness 0.152 --panel 150x100 --energy-max 40 --out energy_sweep.mp4
```

The preview above is the same command with
`--format gif --duration 8 --out docs/animations/energy_sweep.gif`.

What the frame shows:

- **Left**: the plan view of the panel from the impact face. Each
  interface's delamination ellipse is drawn in its own colour, with a line
  style as well beyond the third interface. An interface that has not yet
  delaminated is muted in the key. When the impact was applied on a window
  (`AnalysisConfig.impact_panel`), its outline is drawn. When even the
  largest damage is small on the panel, the view zooms in on it, and the
  title says so.
- **Right**: the projected damage area, the dent depth and the knockdown
  (CAI or TAI) against energy. The curve up to the current energy is solid,
  the rest faint, and a marker follows the energy. The dent panel has the
  0.30 mm detectability threshold (`BVID_DENT_THRESHOLD_MM`) and says
  "now visible" once the dent reaches it.
- **Top**: the material, layup, impactor and panel, the current energy and
  the model tier.

The sweep runs one `BvidAnalysis` per keyframe energy: by default 60
(`empirical`) or 40 (`semi_analytical`) energies from 0.8 × the onset
energy to `--energy-max` (or `--energy`), plus the onset energy itself.
`--energy-min`, `--energy-max` and `--n-energies` change the grid. The clip
holds at the first energy, sweeps to the last and holds there. The `fe3d`
tier is refused: it would solve two finite-element models per energy.

## Damage stack

<video src="../animations/damage_stack.webm" autoplay loop muted playsinline width="640">
  Damage stack preview (WebM): docs/animations/damage_stack.webm
</video>

*Preview (640 × 360 WebM). The 1080p MP4 is about 4.5 MB for 10 s.*

```bash
bvidfe-animate damage-stack --material IM7/8552 \
    --layup 0,45,-45,90,0,45,-45,90,0,45,-45,90,0,45,-45,90 \
    --thickness 0.152 --panel 150x100 --energy 30 --out damage_stack.mp4
```

The preview above was made from that MP4 with
`ffmpeg -i damage_stack.mp4 -vf scale=640:-2:flags=lanczos -c:v libvpx-vp9 -crf 40 -b:v 0 -an docs/animations/damage_stack.webm`.

What the frame shows:

- Every ply as a sheet, the impact face (ply 0) on top. The delamination at
  each interface lies between its two plies, filled in its interface colour
  (the colours of the energy sweep). It elongates along the fibres of the
  ply below it, so the stack shows the helix that the plan view hides.
- Short lines in each ply's fibre direction, on the corner of the ply
  farthest from the camera, where they do not cover the delaminations.
  Stacked, these corners read as the stacking sequence.
- The impactor above the dent on the impact face.
- On the right, the damage area, dent depth, onset energy and each
  interface's ellipse size.

The intro drops the plies into place one by one, from the back face up;
then the camera orbits 360°. `--no-intro` leaves the intro out, for an orbit
that loops. With `--cscan FILE` the stack shows a measured damage state
instead of the impact model's (a C-scan has no impactor, so none is drawn).

## fe3d load ramp

<video src="../animations/cai_ramp.webm" autoplay loop muted playsinline width="640">
  fe3d load ramp preview (WebM): docs/animations/cai_ramp.webm
</video>

*Preview (640 × 360 WebM). The 1080p MP4 is about 1.7 MB for 10 s.*

```bash
bvidfe-animate cai-ramp --material IM7/8552 \
    --layup 0,45,-45,90,0,45,-45,90,0,45,-45,90,0,45,-45,90 \
    --thickness 0.152 --panel 150x100 --energy 20 --in-plane-size 5 \
    --out cai_ramp.mp4
```

The preview above was made from that MP4 with
`ffmpeg -i cai_ramp.mp4 -vf scale=640:-2:flags=lanczos -c:v libvpx-vp9 -crf 40 -b:v 0 -an docs/animations/cai_ramp.webm`.
The render takes about 1.5 min on one core: about 12 s for the two fe3d
solves, then 0.25 s per frame.

The load on the damaged panel rises from zero to its residual strength,
and each element of the fe3d mesh is coloured by its **load ratio to local
failure**, λ / c_e: the applied load over the load at which that element
fails. The first element to reach 1 is where fe3d puts first-ply failure.
`--tier` is not an option: the ramp is always fe3d's, on the mesh of
`--in-plane-size` and `--elements-per-ply` (default 5 mm, one element per
ply). It works for CAI and TAI (`--loading tension`), for an impact or a
C-scan, and for a biaxial load.

What the frame shows:

- **Left, in 3D**: the damaged panel, impact face up, cut open along the
  load through the first element to fail. Only the faces the camera sees
  are drawn (the impact face, the cut face and one end of the block that
  is kept), so the drawing cannot put a face in front of one it is
  behind; the cut face shows every ply. The displacement and the
  thickness are exaggerated, both factors printed under the view. The
  outline on the impact face is the delamination footprint (the union of
  every interface's ellipse: fe3d's damage zone). The first element to
  fail is outlined, in red once it fails, with its ply and angle.
- **Left, below**: top views of the highest ratio through the thickness,
  of the damaged panel and of the undamaged one at the same applied
  stress (it stays far from 1), the footprint and the cut on the first.
  A top view has no depth to sort, so it is the check on the 3D view.
- **Right**: the stress against the end shortening (or extension in TAI)
  of both panels, with the first-ply-failure point; the highest ratio of
  each panel against 1; the residual strength, knockdown, pristine
  strength and governing channel of the fe3d result.
- **Top**: the applied stress, rising frame by frame.

The clip holds unloaded briefly, ramps the stress (eased) and holds 1.5 s
at the end. When **delamination growth governs** (its closed-form stress is
below first-ply failure, as `bvidfe` decides), the ramp ends at the growth
stress and says so: the growth has no FE field, so no element reaches 1,
and the first-ply-failure point is marked "not reached".

`--pointwise` colours by each element's lowest Gauss point instead of the
disk-averaged check that sets first-ply failure; that ratio can pass 1
before the end of the ramp.

The 3D view draws at most 12,000 faces (`cai_ramp.MAX_FACES`; a
150 × 100 mm, 16-ply panel draws about 1,100 at 5 mm and 2,800 at 2.5 mm);
a finer mesh is refused before solving, with the fe3d size cap
(`BVIDFE_FE3D_MAX_DOF`) as for `bvidfe`. Frame time hardly depends on the
mesh: 0.26 s at 5 and at 2.5 mm, 0.30 s at 5,200 faces.

## In the app

The Streamlit app's **Animate** tab renders the same three scenes from the
sidebar's configuration, at 1080p:

1. Pick the scene: *Damage vs impact energy* (Impact event mode; the tab
   picks `empirical` or `semi_analytical` for the sweep), *Exploded
   delamination stack* (impact or C-scan) or *fe3d load ramp to first-ply
   failure* (the sidebar config solved with fe3d, on the sidebar's fe3d
   mesh when the tier is fe3d and the 5 mm default otherwise; the same
   mesh cap as an fe3d run).
2. **Preview stills** runs the model and shows the frames at t = 0, 0.5
   and 1, with an estimate of the render time timed on the server.
3. **Render 1080p MP4** encodes the 10 s clip with a progress bar, then
   plays it, with a **Download MP4** button. A session renders one video
   at a time.

Rendering needs ffmpeg on the server: `requirements.txt` installs
`imageio-ffmpeg`, and on Streamlit Community Cloud `packages.txt` installs
the system ffmpeg (see [Deployment](deployment.md)). Without it the tab
says how to install it; the stills work regardless.

## Output options

| Flag | Default | Meaning |
|---|---|---|
| `--out PATH` | `energy_sweep.mp4`, `damage_stack.mp4`, `cai_ramp.mp4` | Output file. A PNG sequence goes into `<stem>_frames/` (or `PATH` itself when it has no suffix) |
| `--format` | `mp4` | `mp4`, `gif` (a preview), `png` (one file per frame) or `auto` (MP4 if ffmpeg is found, else PNG) |
| `--fps N` | 30 (GIF 12) | Frames per second |
| `--duration S` | 10 | Clip length in seconds |
| `--width`, `--height` | 1920 × 1080 (GIF 640 × 360) | Frame size: 16:9, even. Give one, and the other follows |
| `--preview` | off | Write three stills (t = 0, 0.5 and 1) as `<stem>_t0.00.png` ... instead of the clip |
| `--energy-min`, `--energy-max`, `--n-energies` | see above | The keyframe energies of `energy-sweep` |
| `--no-intro` | off | `damage-stack` without the assembly intro |
| `--pointwise` | off | `cai-ramp` coloured by the lowest Gauss point of each element |

The analysis options are those of `bvidfe` (`--material`, `--layup`,
`--thickness`, `--panel`, `--boundary`, `--loading`,
`--transverse-stress-ratio`, `--energy` or `--cscan`, `--impactor-*`,
`--mass`, `--location`; `--tier` for `energy-sweep`, the fe3d mesh
options `--in-plane-size` and `--elements-per-ply` for `cai-ramp`), with
the same validation. Rendering takes about 0.2 s per 1080p frame, so about a minute
for 10 s at 30 fps, plus a few seconds for the analyses (about 0.4 s per
energy for `semi_analytical`). Smaller frames render faster.

## From Python

```python
from bvidfe.viz.anim import (
    AnimationSpec,
    render_cai_ramp,
    render_damage_stack,
    render_energy_sweep,
)

render_energy_sweep(config, path="energy_sweep.mp4")
render_damage_stack(config, path="damage_stack.mp4", spec=AnimationSpec(duration_s=8))
render_cai_ramp(config, path="cai_ramp.mp4")
```

See [Python API → Animations](python_api.md#animations) for the timeline
and still-preview functions.

## What the frames say, and what they do not

- **The numbers are model outputs.** The model runs once per keyframe
  energy. Between two keyframes only the drawing is interpolated (the
  ellipses' size, the marker's position, an interface fading in as it
  first delaminates). Every number printed on a frame, the energy
  included, is the nearest keyframe's model output. The footer gives the
  number of model runs.
- **The model tier** is on every frame of the energy sweep. The damage
  stack shows the impact model's damage, which is the same for every tier.
  The load ramp is fe3d's, and its footer gives the mesh.
- **The load ramp is linear elastic to first-ply failure**, with no
  progressive damage and nothing after the first failure; fe3d is
  qualitative. Every frame says so. Its frames are exact scalings of one
  solve per panel (the model is linear), not interpolations; the stresses
  are on the scale `bvidfe` reports (fe3d's applied stress times the
  pristine strength over the undamaged panel's first-ply failure, the
  factor printed), so the residual strength and knockdown on the frame are
  `BvidAnalysis(config).run()`'s with `tier="fe3d"`. The deformation and the
  thickness are exaggerated by the printed factors.
- **Below the onset energy** the plan view says "below damage onset
  (x.x J)", and so does the stack when there is no delamination.
- **Biaxial loads** (`--transverse-stress-ratio`, β ≠ 0) carry
  "experimental, not validated" on every frame, and TAI "not validated".
- **The stack is not to scale through the thickness.** The plies are spread
  apart (the frame gives the factor) and the dent is drawn on the impact
  face only, its depth shortened. In plan, ellipses, panel and impactor are
  to scale.
- **The delamination shapes are the model's**: one ellipse per interface
  from the impact mapping's template, not a measured or simulated shape.
  A C-scan input is drawn as supplied.
- **The same inputs give the same frames** on a given platform.
  Videos are made on demand and are not kept in the repository; the
  previews on this page are reduced-size copies.

## Citation

```text
Elhajjar, R. (2026). BVID-FE: Barely Visible Impact Damage residual-strength
analysis for composite laminates. https://github.com/elhajjar1/BVID-FE
```

The BibTeX entry is in the [README](https://github.com/elhajjar1/BVID-FE#citation),
and `CITATION.cff` is in the repository. The footer text comes from
`bvidfe.citation`, which a test keeps in step with `CITATION.cff`.
