# -*- mode: python ; coding: utf-8 -*-
# BvidFE.spec — PyInstaller spec for the BVID-FE desktop app
# Build: pyinstaller BvidFE.spec --noconfirm --clean
#
# Output (macOS):  dist/BVID-FE.app
# Output (Windows): dist/BVID-FE/

import sys
import os

from PyInstaller.utils.hooks import collect_submodules, collect_data_files

block_cipher = None

# Project root (directory containing this spec file)
PROJECT_ROOT = os.path.dirname(os.path.abspath(SPEC))

# Collect all scipy submodules to avoid dynamic-import misses at runtime
scipy_hiddenimports = collect_submodules("scipy")

hidden_imports = (
    # ---- bvidfe package — explicit listing avoids missed sub-packages ----
    [
        "bvidfe",
        "bvidfe.core",
        "bvidfe.core.material",
        "bvidfe.core.laminate",
        "bvidfe.core.geometry",
        "bvidfe.impact",
        "bvidfe.impact.olsson",
        "bvidfe.impact.shape_templates",
        "bvidfe.impact.dent_model",
        "bvidfe.impact.mapping",
        "bvidfe.damage",
        "bvidfe.damage.state",
        "bvidfe.damage.io",
        "bvidfe.elements",
        "bvidfe.elements.hex8",
        "bvidfe.elements.hex8i",
        "bvidfe.elements.gauss",
        "bvidfe.elements.cohesive",
        "bvidfe.solver",
        "bvidfe.solver.static",
        "bvidfe.solver.assembler",
        "bvidfe.solver.boundary",
        "bvidfe.solver.buckling",
        "bvidfe.failure",
        "bvidfe.failure.larc05",
        "bvidfe.failure.tsai_wu",
        "bvidfe.failure.soutis_openhole",
        "bvidfe.failure.evaluator",
        "bvidfe.analysis",
        "bvidfe.analysis.config",
        "bvidfe.analysis.bvid",
        "bvidfe.analysis.results",
        "bvidfe.analysis.semi_analytical",
        "bvidfe.analysis.fe_tier",
        "bvidfe.analysis.fe_mesh",
        "bvidfe.viz",
        "bvidfe.viz.plots_2d",
        "bvidfe.viz.plots_3d",
        "bvidfe.viz.style",
        "bvidfe.sweep",
        "bvidfe.sweep.parametric_sweep",
        "bvidfe.cli",
        "bvidfe.gui",
        "bvidfe.gui.app",
        "bvidfe.gui.main_window",
        "bvidfe.gui.config_io",
        "bvidfe.gui.workers",
        "bvidfe.gui.panels",
        "bvidfe.gui.panels.analysis_panel",
        "bvidfe.gui.panels.damage_panel",
        "bvidfe.gui.panels.impact_panel",
        "bvidfe.gui.panels.input_mode_panel",
        "bvidfe.gui.panels.material_panel",
        "bvidfe.gui.panels.panel_panel",
        "bvidfe.gui.panels.sweep_panel",
        "bvidfe.gui.tabs",
        "bvidfe.gui.tabs.damage_map_tab",
        "bvidfe.gui.tabs.knockdown_tab",
        "bvidfe.gui.tabs.placeholder_tab",
        "bvidfe.gui.tabs.summary_tab",
        # v0.2.0-dev additions — lazy-imported inside BvidMainWindow.__init__
        # so PyInstaller's static analysis may miss them without this hint:
        "bvidfe.gui.tabs.buckling_tab",
        "bvidfe.gui.tabs.stress_field_tab",
        "bvidfe.gui.tabs.mesh_3d_tab",
    ]
    # ---- matplotlib backends ----
    + [
        "matplotlib",
        "matplotlib.backends",
        "matplotlib.backends.backend_qtagg",
        "matplotlib.backends.backend_qt5agg",
        "matplotlib.backends.backend_agg",
        "matplotlib.backends.backend_svg",
        "matplotlib.backends.backend_pdf",
        "mpl_toolkits",
        "mpl_toolkits.mplot3d",
    ]
    # ---- PyQt6 ----
    + [
        "PyQt6",
        "PyQt6.QtCore",
        "PyQt6.QtGui",
        "PyQt6.QtWidgets",
        "PyQt6.sip",
    ]
    # ---- pyvista / vtk (collect_submodules catches dynamic vtk plugin loading) ----
    + collect_submodules("pyvista")
    + collect_submodules("vtkmodules")
    + [
        "pyvista",
        "pyvistaqt",
    ]
    # ---- shapely ----
    + collect_submodules("shapely")
    # ---- scipy (full submodule tree) ----
    + scipy_hiddenimports
    # ---- array API compat shim (scipy 1.11+) ----
    + [
        "scipy._lib.array_api_compat.numpy.fft",
    ]
    # ---- pandas ----
    + collect_submodules("pandas")
    # ---- packaging / pkg_resources ----
    + [
        "appdirs",
        "packaging",
        "packaging.version",
        "packaging.specifiers",
        "packaging.requirements",
        "pkg_resources",
    ]
)

datas = (
    collect_data_files("pyvista", include_py_files=False)
    + [
        (os.path.join(PROJECT_ROOT, "README.md"), "."),
    ]
)

a = Analysis(
    [os.path.join(PROJECT_ROOT, "src", "bvidfe", "gui", "app.py")],
    pathex=[os.path.join(PROJECT_ROOT, "src")],
    binaries=[],
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "_tkinter",
        "PyQt5",
        "PySide6",
        "PySide2",
        "IPython",
        "jupyter",
        "notebook",
        "pytest",
        "sphinx",
        # NOTE: do NOT exclude pydoc / xmlrpc / docutils —
        # scipy._lib._docscrape imports pydoc, and rich-rst pulls in
        # docutils. Excluding them breaks the windowed bundle at
        # "import scipy" time (ModuleNotFoundError: No module named 'pydoc').
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="BVID-FE",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    windowed=True,
    disable_windowed_traceback=False,
    argv_emulation=True,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="BVID-FE",
)

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="BVID-FE.app",
        icon=None,
        bundle_identifier="edu.uwm.elhajjar.bvidfe",
        info_plist={
            "CFBundleName": "BVID-FE",
            "CFBundleDisplayName": "BVID-FE",
            # Read version from installed bvidfe package so it stays in sync
            "CFBundleShortVersionString": __import__("bvidfe").__version__,
            "CFBundleVersion": __import__("bvidfe").__version__,
            "NSHighResolutionCapable": True,
            "NSPrincipalClass": "NSApplication",
            "LSMinimumSystemVersion": "12.0",
            "CFBundleDocumentTypes": [
                {
                    "CFBundleTypeName": "BVID-FE Project",
                    "CFBundleTypeExtensions": ["json"],
                    "CFBundleTypeRole": "Viewer",
                }
            ],
        },
    )
