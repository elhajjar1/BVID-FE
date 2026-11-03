# Deployment

BVID-FE ships with a Streamlit web app (`app.py`) that wraps
`bvidfe.analysis.BvidAnalysis`, and a `requirements.txt` tuned for the
cloud environment. This page shows how to publish it at a public URL on
Streamlit Community Cloud.

## Local smoke test

```bash
pip install -r requirements.txt
pip install -e .
streamlit run app.py
```

A browser tab should open at <http://localhost:8501>. Tweak parameters in
the sidebar, click **Run analysis**, and verify the result tabs
(Summary, Damage Map, Knockdown Curve, 3D Damage, Buckling, Damage
Severity, Sweep, Export NCR) populate.

## Streamlit Community Cloud

1. Push your branch to GitHub (`elhajjar1/BVID-FE`). You can deploy a
   feature branch directly or merge it to `main` first.
2. Sign in at <https://share.streamlit.io> with **Sign in with GitHub**. On
   first sign-in, grant access to the GitHub account that owns the repo
   (`elhajjar1`).
3. **Create app** → **Deploy a public app from GitHub**:
    - **Repository:** `elhajjar1/BVID-FE`
    - **Branch:** `main` (or your feature branch)
    - **Main file path:** `app.py`
    - **App URL:** pick a sub-domain (e.g. `bvidfe.streamlit.app`)
4. In **Advanced settings** set Python to **3.11**. It satisfies
   `pyproject.toml`'s `requires-python = ">=3.10"`, is in the CI test
   matrix (3.10–3.12), and is the version the CI validation job uses.
5. Click **Deploy**.

Streamlit Cloud then clones the repo, creates a virtualenv, runs
`pip install -r requirements.txt` and starts `streamlit run app.py`. It
does not necessarily install the `bvidfe` package itself, so `app.py` adds
the repository's `src/` directory to `sys.path` before importing it. The
first build takes 2–4 minutes (numpy/scipy/shapely wheels); later
redeploys reuse the cache and finish in about 30 s.

Once the build finishes the URL is live (e.g.
`https://bvidfe.streamlit.app`); no GitHub or Streamlit account is needed
to use a public app.

Every push to the configured branch triggers an automatic redeploy. The
**Manage app** panel (bottom-right of the live app) shows the live logs and
lets you reboot the app, edit secrets and other settings, or pause or
delete it.

## Resource limits

Streamlit Cloud (Community tier) provisions per app:

- 1 GB RAM
- 1 vCPU
- shared, ephemeral filesystem
- no GPU
- public repos only by default

Implications for BVID-FE:

- The **fe3d** tier solves a 3D FE problem and can easily blow past 1 GB
  on dense meshes. It also solves twice per run (damaged and undamaged
  panel, to normalise the knockdown): the default 8-ply 150 × 100 mm
  panel at the default 5 mm mesh (4,800 elements, about 17,600 DOF)
  takes about 5 s (about 35 s at 2.5 mm), and thick laminates or finer
  meshes take longer. The Streamlit
  app defaults to the **empirical** tier and refuses an fe3d run whose
  estimated mesh exceeds `BVIDFE_FE3D_MAX_DOF` (150,000 degrees of freedom
  by default). An 8-ply 150 × 100 mm panel peaks at about 0.2 GB at 5 mm
  and 0.5 GB at 2.5 mm, but the direct-solver fallback needs ~1.3 GB at
  2.5 mm. If you allow fe3d in the cloud, keep `elements_per_ply = 1` and
  `in_plane_size_mm >= 5` (the app's defaults).
- Don't write user-supplied data anywhere except `tempfile`-scoped paths;
  the filesystem resets on every reboot.

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| Build fails on `PyQt6` / `pyvista` | They should not be in `requirements.txt`. Streamlit Cloud has no display server. |
| App boots but `import bvidfe` fails | `app.py` adds `src/` to `sys.path`, so check that `src/bvidfe/` is on the deployed branch, and that `pip install -e .` works locally. |
| Plots don't render | Make sure `matplotlib.use("Agg")` is called before any `pyplot` import. `app.py` already does this. |
| Slow first run | Expected — numpy/scipy/shapely wheels are ~100 MB combined. Subsequent runs are cached. |
| Out of memory on large fe3d run | Increase `in_plane_size_mm`, keep `elements_per_ply = 1`, or switch to the empirical / semi_analytical tier. |
| `FE3DSizeError` or "exceeding the safe cap" | The mesh exceeds `BVIDFE_FE3D_MAX_DOF` (default 150,000). Coarsen the mesh — Streamlit Cloud doesn't have the headroom to raise this cap. |

## Optional next steps

- The Python version is the one picked under **Advanced settings** at
  deploy time; see Streamlit's documentation for changing it later.
- Configure a custom domain in Streamlit Cloud's settings.
- Add Streamlit secrets via the Manage panel if the app later needs API keys.
- CI already imports `app.py` on every push (`tests/test_streamlit_app.py`);
  a fuller smoke test could click **Run analysis** through
  `streamlit.testing.v1.AppTest`.
