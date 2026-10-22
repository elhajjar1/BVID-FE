# Deploying BVID-FE to Streamlit Community Cloud

This repo includes `app.py` (a Streamlit wrapper around
`bvidfe.analysis.BvidAnalysis`) and a `requirements.txt` tuned for the
cloud environment. Follow the steps below to publish a public URL that
your engineers can hit from a browser.

## 1. Pre-flight checks

Before deploying, confirm everything works locally:

```bash
pip install -r requirements.txt
pip install -e .
streamlit run app.py
```

A browser tab should open at <http://localhost:8501>. Tweak parameters in
the sidebar, click **Run analysis**, and verify the result tabs
(Summary, Damage Map, Knockdown Curve, 3D Damage, Buckling, Damage
Severity, Sweep, Export NCR) populate.

## 2. Push to GitHub

Streamlit Cloud reads from a GitHub repository. Make sure your branch is
pushed to GitHub and you know the branch name. For this repo:

```
elhajjar1/BVID-FE
```

If you've been working on a feature branch, either deploy directly from
that branch or merge to `main` first.

## 3. Sign in to Streamlit Community Cloud

1. Go to <https://share.streamlit.io>.
2. Click **Sign in with GitHub** and authorize the Streamlit app.
3. On first sign-in, grant access to the GitHub org/account that owns
   the repo (`elhajjar1`).

## 4. Create the app

1. Click **Create app** (top right) → **Deploy a public app from GitHub**.
2. Fill in:
   - **Repository:** `elhajjar1/BVID-FE`
   - **Branch:** `main` (or your feature branch)
   - **Main file path:** `app.py`
   - **App URL:** pick a sub-domain like `bvidfe.streamlit.app`
3. Click **Advanced settings** and set Python version to **3.11**. It
   satisfies `pyproject.toml`'s `requires-python = ">=3.10"`, is in the CI
   test matrix (3.10–3.12), and is the version the CI validation job uses.
4. Click **Deploy**.

Streamlit Cloud will:
- clone the repo
- create a virtualenv
- run `pip install -r requirements.txt`
- start `streamlit run app.py`

Streamlit Cloud does not necessarily install the `bvidfe` package itself,
so `app.py` adds the repository's `src/` directory to `sys.path` before
importing it.

The first build typically takes 2–4 minutes (numpy/scipy/shapely wheels).
Subsequent redeploys reuse the cache and finish in ~30 s.

## 5. Test the live URL

Once the build finishes the URL becomes live (e.g.
`https://bvidfe.streamlit.app`). Send it to your engineers — no GitHub
or Streamlit account is needed to use a public app.

## 6. Iterate

Every push to the configured branch triggers an automatic redeploy. Use
the **Manage app** panel (bottom-right gear in the live app) to:
- view live logs
- reboot the app if it gets stuck
- edit secrets and other app settings
- pause / delete the app

## 7. Resource limits to know about

Streamlit Cloud (Community tier) gives every app:
- **1 GB RAM**
- **1 vCPU**
- **shared filesystem**, ephemeral
- no GPU
- public-only repo by default (private requires a paid plan)

Implications for BVID-FE:
- The **fe3d** tier solves a 3D FE problem and can easily blow past 1 GB
  on dense meshes. It also solves twice per run (damaged and undamaged
  panel, to normalise the knockdown): the default 8-ply 150 × 100 mm
  panel at the default 5 mm mesh (4,800 elements, about 17,600 DOF)
  takes about 5 s (about 30 s at 2.5 mm), and thick laminates or finer
  meshes take longer. The Streamlit
  app defaults to the **empirical** tier and refuses an fe3d run whose
  estimated mesh exceeds `BVIDFE_FE3D_MAX_DOF` (150,000 degrees of freedom
  by default). An 8-ply 150 × 100 mm panel peaks at about 0.2 GB at 5 mm
  and 0.5 GB at 2.5 mm, but the direct-solver fallback needs ~1.3 GB at
  2.5 mm. If you want to allow fe3d in the cloud, keep
  `elements_per_ply = 1` and `in_plane_size_mm >= 5` (the app's defaults).
- Don't write user-supplied data anywhere except `tempfile`-scoped paths;
  the filesystem resets on every reboot.

## 8. Troubleshooting

| Symptom | Fix |
|---|---|
| Build fails on `PyQt6` / `pyvista` | They should no longer be in `requirements.txt`. Streamlit Cloud has no display server. |
| App boots but `import bvidfe` fails | `app.py` adds `src/` to `sys.path`, so check that `src/bvidfe/` is on the deployed branch, and that `pip install -e .` works locally. |
| Plots don't render | Make sure `matplotlib.use("Agg")` is called before any `pyplot` import. `app.py` already does this. |
| Slow first run | Expected — numpy/scipy/shapely wheels are ~100 MB combined. Subsequent runs are cached. |
| Out of memory on large fe3d run | Increase `in_plane_size_mm`, keep `elements_per_ply = 1`, or switch to the empirical / semi_analytical tier. |
| `FE3DSizeError` or "exceeding the safe cap" | The mesh exceeds `BVIDFE_FE3D_MAX_DOF` (default 150,000). Coarsen the mesh — Streamlit Cloud doesn't have the headroom to raise this cap. |

## 9. Optional next steps

- The Python version is the one picked under **Advanced settings** at
  deploy time; see Streamlit's documentation for changing it on a
  deployed app.
- Add Streamlit secrets via the Manage panel if the app needs API keys
  later.
- Configure a custom domain in Streamlit Cloud's settings.
- CI already imports `app.py` on every push (`tests/test_streamlit_app.py`);
  a fuller smoke test could click **Run analysis** through
  `streamlit.testing.v1.AppTest`.
