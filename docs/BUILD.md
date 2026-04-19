# Building BVID-FE executables (Windows + macOS)

This document covers how to produce a standalone desktop executable for
Windows (`BVID-FE.exe`) or macOS (`BVID-FE.app`). The source-install
Python flow (`pip install bvidfe` → `bvidfe` CLI / `bvidfe-gui`) still
works on all three major OSes; this doc is only about the bundled
executable.

**PyInstaller can't cross-compile.** You cannot build a Windows `.exe`
on macOS or vice versa. To get a Windows binary while working on a Mac
(or a Mac `.app` while working on a PC), use GitHub Actions (Options 1 &
2). If you have the target OS physically available, use Option 3.

## Option 1 — GitHub Actions, artifact only (recommended for iteration)

Use this when you want a Windows `.exe` or Mac `.app` to test, without
creating a permanent GitHub Release.

1. Push your working branch to GitHub.
2. Navigate to `Actions` → **build-artifacts** in the repo UI.
3. Click **Run workflow**, pick your branch, and hit the green **Run workflow** button.
4. Wait ~10 minutes. Two matrix jobs (macOS + Windows) build in parallel.
5. Download the artifacts from the run page:
   - `BVID-FE-macOS` — unzips to `BVID-FE.app`
   - `BVID-FE-Windows` — unzips to a folder containing `BVID-FE.exe`

Artifacts are retained for 90 days on the Actions run page. They are NOT
attached to a permanent GitHub Release.

Alternative trigger: push to any branch whose name starts with `build-`
(e.g. `build-test1`) and the workflow runs automatically.

## Option 2 — GitHub Actions, tagged release (permanent download)

Use this for public releases that you want to ship with release notes.

1. Bump `src/bvidfe/__init__.py` → `__version__` and update `CHANGELOG.md`.
2. Commit and push.
3. Tag and push the tag:
       git tag v0.2.0
       git push origin v0.2.0
4. The `release` workflow builds for macOS + Windows, creates a GitHub
   Release at that tag, and attaches:
   - `BVID-FE-v0.2.0-macOS.zip`
   - `BVID-FE-v0.2.0-Windows.zip`

If you want to test the release flow WITHOUT cutting a real tag:

1. Go to `Actions` → **release** → **Run workflow**.
2. Enter a fake tag like `v0.2.0-dev-manual` and hit Run.
3. The workflow produces a **pre-release** at that tag name (safely
   marked "pre-release" so it doesn't look like a production release).

## Option 3 — Local build (you have the target OS)

### Windows

Requires Python 3.9+ on PATH (`python --version`). From `cmd.exe` or
PowerShell in the repo root:

    scripts\build_windows.bat

The script:

1. Creates a `.venv\` if missing.
2. Runs `pip install -e ".[all]"` + `pyinstaller`.
3. Runs the test suite headlessly (offscreen Qt) as a smoke check.
4. Invokes `pyinstaller BvidFE.spec --noconfirm --clean`.

Output: `dist\BVID-FE\BVID-FE.exe` (and supporting DLLs in the same
folder). Zip the entire `dist\BVID-FE\` folder to distribute.

### macOS

Requires Python 3.9+ (`brew install python@3.11` or the official
installer). From the repo root:

    bash scripts/build_macos.sh

Output: `dist/BVID-FE.app`. Distribute with:

    cd dist && zip -r BVID-FE-macOS.zip BVID-FE.app

Recipients who download the zip from a browser may see the macOS
quarantine-flag gatekeeper warning. They can clear it once with:

    xattr -rd com.apple.quarantine BVID-FE.app

## Troubleshooting

**"SmartScreen prevented an unrecognized app from starting"** (Windows)
The build is unsigned. Click **More info** → **Run anyway**. Signing
the `.exe` requires a code-signing certificate and is deferred past
v0.2.0.

**"BVID-FE.app can't be opened because Apple cannot check it for
malicious software"** (macOS) — same story. Clear the quarantine flag
as above, or right-click → Open → Open (allows one-time open).

**The app launches but crashes on startup** — the most common cause is
a missing hidden import. Add it to `BvidFE.spec` under `hidden_imports`
and rebuild. Check `dist/BVID-FE/warn-BVID-FE.txt` (Windows) or
`build/BVID-FE/warn-BVID-FE.txt` (macOS) for missing-module warnings.

**The Windows app is very large (>500 MB)** — PyInstaller bundles scipy,
pandas, matplotlib, pyvista, VTK, etc. This is expected. UPX
compression (`upx=True` in the spec) is enabled; toggle to `upx=False`
if you have UPX issues.

## Why no cross-compile?

PyInstaller embeds the Python interpreter for the runtime OS + architecture.
There is no supported way to produce a Windows `.exe` on macOS/Linux or
vice versa. This is not a BVID-FE limitation — it applies to every
PyInstaller-packaged Python app. GitHub Actions' free matrix runners
solve this cleanly.
