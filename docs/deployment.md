# Deployment

BVID-FE ships with a Streamlit web app (`app.py`) that wraps
`bvidfe.analysis.BvidAnalysis`. This page summarises how to publish a public
URL on Streamlit Community Cloud — the full guide lives in
[`DEPLOYMENT_STREAMLIT.md`](https://github.com/elhajjar1/BVID-FE/blob/main/DEPLOYMENT_STREAMLIT.md).

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

1. Push your branch (typically `main`) to GitHub.
2. Sign in at <https://share.streamlit.io> with the GitHub account that owns
   the repo (`elhajjar1`).
3. **Create app** → **Deploy a public app from GitHub**:
    - **Repository:** `elhajjar1/BVID-FE`
    - **Branch:** `main`
    - **Main file path:** `app.py`
    - **App URL:** pick a sub-domain (e.g. `bvidfe.streamlit.app`)
4. In **Advanced settings** set Python to **3.11**.
5. Click **Deploy**. First build takes 2–4 minutes; subsequent redeploys
   finish in ~30 s.

Every push to the configured branch triggers an automatic redeploy.

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
  takes about 15-30 s, and thick laminates or finer meshes take
  minutes. The Streamlit
  app defaults to the **empirical** tier and refuses an fe3d run whose
  estimated mesh exceeds `BVIDFE_FE3D_MAX_DOF` (500,000 degrees of freedom
  by default). If you allow fe3d in the cloud, keep `elements_per_ply = 1`
  and `in_plane_size_mm >= 5` (the app's defaults).
- Don't write user-supplied data anywhere except `tempfile`-scoped paths;
  the filesystem resets on every reboot.

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| Build fails on `PyQt6` / `pyvista` | They should not be in `requirements.txt`. Streamlit Cloud has no display server. |
| App boots but `import bvidfe` fails | `app.py` adds `src/` to `sys.path`, so check that `src/bvidfe/` is on the deployed branch, and that `pip install -e .` works locally. |
| Plots don't render | Make sure `matplotlib.use("Agg")` is called before any `pyplot` import. `app.py` already does this. |
| Slow first run | Expected — numpy/scipy/shapely wheels are ~100 MB combined. |
| Out of memory on large fe3d run | Increase `in_plane_size_mm`, keep `elements_per_ply = 1`, or switch to the empirical / semi_analytical tier. |
| `FE3DSizeError` or "exceeding the safe cap" | The mesh exceeds `BVIDFE_FE3D_MAX_DOF` (default 500,000). Coarsen the mesh — Streamlit Cloud doesn't have the headroom to raise this cap. |

## Optional next steps

- The Python version is the one picked under **Advanced settings** at
  deploy time; see Streamlit's documentation for changing it later.
- Configure a custom domain in Streamlit Cloud's settings.
- Add Streamlit secrets via the Manage panel if the app later needs API keys.
- CI already imports `app.py` on every push (`tests/test_streamlit_app.py`);
  a fuller smoke test could click **Run analysis** through
  `streamlit.testing.v1.AppTest`.
