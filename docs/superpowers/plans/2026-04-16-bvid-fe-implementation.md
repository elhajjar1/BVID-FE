# BVID-FE Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `bvidfe` — a Python package + PyQt6 desktop app that predicts residual compression-after-impact (CAI) and tension-after-impact (TAI) strength of composite laminates containing Barely Visible Impact Damage, per `docs/superpowers/specs/2026-04-16-bvid-fe-design.md`.

**Architecture:** Modular Python package with strict layer ordering `core → impact → damage → elements → solver → failure → analysis → (viz, sweep, gui)`. Three modeling tiers (empirical / semi-analytical / 3D FE) dispatched from a single `BvidAnalysis(AnalysisConfig).run()` entry point. Two workflow paths (impact-driven and inspection-driven) converge on a shared `DamageState` data model.

**Tech Stack:** Python 3.9+, NumPy, SciPy (sparse linalg + eigsh), Matplotlib, PyQt6, pyvista/pyvistaqt, pytest + pytest-qt, PyInstaller.

**Source conventions:** `src/bvidfe/` layout (matches WrinkleFE at `/Users/elhajjar/Library/CloudStorage/OneDrive-UWM/AI/Double_Wrinkle/wrinklefe/src/wrinklefe`). Port CLT, hex8, hex8i, gauss quadrature, and orthotropic material from WrinkleFE wherever cleanly reusable; otherwise write from scratch against this plan.

**Reference repos (read-only):**
- WrinkleFE: `/Users/elhajjar/Library/CloudStorage/OneDrive-UWM/AI/Double_Wrinkle/wrinklefe`
- PorosityFE: `/Users/elhajjar/Library/CloudStorage/OneDrive-UWM/AI/Porosity_FE_App`

**Conventions for every task:**
- Follow @superpowers:test-driven-development — failing test first, minimal implementation, test passes, commit.
- Follow @superpowers:verification-before-completion — run the tests and confirm output before checking a step off.
- Use Python type hints everywhere, dataclasses for value objects, `Literal[...]` for string enums.
- Line length 100, `black`-formatted, `ruff` clean. Configured in phase 1.
- Import order: stdlib → third-party → first-party (`bvidfe`).

**Phase map:**
| Phase | Scope | Deliverable |
|---|---|---|
| 1 | Project scaffold | pyproject, CI, dev deps, empty package importable |
| 2 | core/material, core/laminate, core/geometry | Orthotropic materials, CLT ABD, panel geometry |
| 3 | damage/state, damage/io | `DamageState` + `DelaminationEllipse` + JSON round-trip |
| 4 | impact/olsson, impact/shape_templates, impact/dent_model, impact/mapping | Full `ImpactEvent → DamageState` forward mapping |
| 5 | failure/soutis_openhole, failure/tsai_wu, failure/larc05, failure/evaluator | Failure criteria library |
| 6 | analysis/config + empirical tier + analysis/bvid | First end-to-end path: `BvidAnalysis(config).run()` returns empirical CAI/TAI |
| 7 | elements (hex8, hex8i, gauss, cohesive) | Element library for the FE tier |
| 8 | solver (static, assembler, boundary, buckling) | Static FE + linear buckling eigensolve |
| 9 | Semi-analytical tier (Rayleigh-Ritz sublaminate buckling + Soutis notch) | Middle-tier CAI/TAI wired into `BvidAnalysis` |
| 10 | 3D FE tier | Full-fidelity CAI/TAI with cohesive surfaces wired into `BvidAnalysis` |
| 11 | viz (plots_2d, plots_3d, style) | Damage maps, knockdown curves, mesh + mode shape + stress |
| 12 | sweep | Parametric sweeps + CSV export |
| 13 | GUI (panels, workers, main window) | PyQt6 desktop app |
| 14 | Validation (datasets + script) | Soutis / Caprino / Sanchez-Saez / NASA against spec MAE targets |
| 15 | Packaging (PyInstaller spec + CLI entry points) | Standalone app + `bvidfe` and `bvidfe-gui` commands |
| 16 | Docs (README, ARCHITECTURE, CLAUDE.md, screenshots) | Ship-ready repository |

---

## Phase 1 — Project Scaffold

### Task 1.1: Package skeleton

**Files:**
- Create: `pyproject.toml`
- Create: `requirements.txt`
- Create: `.gitignore`
- Create: `src/bvidfe/__init__.py`
- Create: `tests/__init__.py`
- Create: `tests/test_package.py`

- [ ] **Step 1: Write the failing test**

`tests/test_package.py`:
```python
def test_bvidfe_importable():
    import bvidfe
    assert bvidfe.__version__ == "0.1.0.dev0"
```

- [ ] **Step 2: Run test — expect failure (module not found)**

```
pytest tests/test_package.py -v
```

- [ ] **Step 3: Write pyproject.toml, requirements.txt, .gitignore**

`pyproject.toml`:
```toml
[build-system]
requires = ["setuptools>=64", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "bvidfe"
version = "0.1.0.dev0"
description = "Barely Visible Impact Damage residual-strength analysis for composite laminates"
readme = "README.md"
license = {text = "MIT"}
requires-python = ">=3.9"
authors = [{name = "Rani Elhajjar"}]
dependencies = [
    "numpy>=1.24",
    "scipy>=1.10",
    "matplotlib>=3.7",
]

[project.optional-dependencies]
gui = ["PyQt6>=6.5", "pyvista>=0.42", "pyvistaqt>=0.11"]
dev = ["pytest>=7", "pytest-qt>=4.3", "ruff>=0.5", "black>=24"]
all = ["bvidfe[gui,dev]"]

[project.scripts]
bvidfe = "bvidfe.cli:main"
bvidfe-gui = "bvidfe.gui.app:main"

[tool.setuptools.packages.find]
where = ["src"]

[tool.black]
line-length = 100
target-version = ["py39"]

[tool.ruff]
line-length = 100
target-version = "py39"

[tool.pytest.ini_options]
testpaths = ["tests"]
```

`requirements.txt`:
```
numpy>=1.24
scipy>=1.10
matplotlib>=3.7
PyQt6>=6.5
pyvista>=0.42
pyvistaqt>=0.11
```

`.gitignore` (port from WrinkleFE): `__pycache__/`, `*.pyc`, `.pytest_cache/`, `build/`, `dist/`, `*.egg-info/`, `.venv/`, `.DS_Store`, `validation/figures/*.png`, `screenshots/*.png`.

`src/bvidfe/__init__.py`:
```python
"""BVID-FE — Barely Visible Impact Damage analysis for composite laminates."""
__version__ = "0.1.0.dev0"
```

- [ ] **Step 4: Install editable and run tests**

```
pip install -e ".[all]"
pytest tests/test_package.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml requirements.txt .gitignore src/bvidfe/__init__.py tests/
git commit -m "Scaffold bvidfe package with pyproject.toml and smoke test"
```

### Task 1.2: CI workflow + ruff/black config

**Files:**
- Create: `.github/workflows/tests.yml`
- Create: `.pre-commit-config.yaml` (optional but recommended; add only if repo uses pre-commit)

- [ ] **Step 1: Write GitHub Actions CI**

`.github/workflows/tests.yml`:
```yaml
name: tests
on: [push, pull_request]
jobs:
  test:
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        python-version: ["3.9", "3.10", "3.11", "3.12"]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install -e ".[dev]"
      - run: ruff check src tests
      - run: black --check src tests
      - run: pytest -v
```

- [ ] **Step 2: Verify locally**

```
ruff check src tests
black --check src tests
pytest -v
```
All should pass.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/tests.yml
git commit -m "Add CI workflow: ruff, black, pytest on macOS/Ubuntu/Windows"
```

---

## Phase 2 — Core Modules

### Task 2.1: `core/material.py` — OrthotropicMaterial + MaterialLibrary

**Files:**
- Create: `src/bvidfe/core/__init__.py`
- Create: `src/bvidfe/core/material.py`
- Create: `tests/core/__init__.py`
- Create: `tests/core/test_material.py`

- [ ] **Step 1: Write failing tests**

`tests/core/test_material.py`:
```python
import pytest
from bvidfe.core.material import OrthotropicMaterial, MATERIAL_LIBRARY

def test_material_library_has_four_presets():
    for name in ["AS4/3501-6", "IM7/8552", "T700/2510", "T800/epoxy"]:
        assert name in MATERIAL_LIBRARY
        m = MATERIAL_LIBRARY[name]
        assert isinstance(m, OrthotropicMaterial)
        assert m.E11 > m.E22 > 0

def test_orthotropic_material_stiffness_matrix_is_symmetric():
    m = MATERIAL_LIBRARY["IM7/8552"]
    C = m.get_stiffness_matrix()   # 6x6
    assert C.shape == (6, 6)
    # Symmetric within tolerance
    import numpy as np
    assert np.allclose(C, C.T, atol=1e-6)

def test_material_rejects_negative_modulus():
    with pytest.raises(ValueError):
        OrthotropicMaterial(
            name="bad", E11=-1.0, E22=10.0, nu12=0.3,
            G12=5.0, G13=5.0, G23=3.0,
            Xt=1.0, Xc=1.0, Yt=1.0, Yc=1.0, S12=1.0, S23=1.0,
            G_Ic=0.1, G_IIc=0.5, rho=1.6e-6,
        )
```

- [ ] **Step 2: Run — expect import failure**

```
pytest tests/core/test_material.py -v
```

- [ ] **Step 3: Implement OrthotropicMaterial**

`src/bvidfe/core/material.py` (approx):
```python
from dataclasses import dataclass
import numpy as np

@dataclass
class OrthotropicMaterial:
    name: str
    E11: float; E22: float; nu12: float
    G12: float; G13: float; G23: float
    Xt: float; Xc: float; Yt: float; Yc: float
    S12: float; S23: float
    G_Ic: float; G_IIc: float
    rho: float
    # Impact-mapping calibration (defaults chosen for typical CFRP; per-material override):
    olsson_alpha: float = 0.8
    dent_beta: float = 0.05
    dent_gamma: float = 0.5
    fiber_break_eta: float = 0.0
    fiber_break_E_threshold: float = 1e9   # effectively disabled unless overridden
    soutis_k_s: float = 2.5
    soutis_m: float = 0.5
    wn_d0_mm: float = 1.0                  # Whitney-Nuismer characteristic distance

    def __post_init__(self):
        for k in ("E11","E22","G12","G13","G23","Xt","Xc","Yt","Yc","S12","S23","G_Ic","G_IIc","rho"):
            if getattr(self, k) <= 0:
                raise ValueError(f"{k} must be > 0 (got {getattr(self, k)})")
        if not (-1 < self.nu12 < 0.5):
            raise ValueError(f"nu12 out of physical range (got {self.nu12})")

    @property
    def nu21(self) -> float:
        return self.nu12 * self.E22 / self.E11

    def get_compliance_matrix(self) -> np.ndarray:
        S = np.zeros((6, 6))
        S[0,0] = 1 / self.E11
        S[1,1] = 1 / self.E22
        S[2,2] = 1 / self.E22
        S[0,1] = S[1,0] = -self.nu12 / self.E11
        S[0,2] = S[2,0] = -self.nu12 / self.E11
        S[1,2] = S[2,1] = -self.nu21 / self.E22
        S[3,3] = 1 / self.G23
        S[4,4] = 1 / self.G13
        S[5,5] = 1 / self.G12
        return S

    def get_stiffness_matrix(self) -> np.ndarray:
        return np.linalg.inv(self.get_compliance_matrix())


MATERIAL_LIBRARY: dict[str, OrthotropicMaterial] = {
    "AS4/3501-6": OrthotropicMaterial(
        name="AS4/3501-6",
        E11=138000, E22=9000, nu12=0.30, G12=6900, G13=6900, G23=3450,
        Xt=2280, Xc=1440, Yt=57, Yc=228, S12=71, S23=50,
        G_Ic=0.26, G_IIc=1.0, rho=1.58e-6,
    ),
    "IM7/8552": OrthotropicMaterial(
        name="IM7/8552",
        E11=165000, E22=8400, nu12=0.34, G12=5600, G13=5600, G23=2800,
        Xt=2560, Xc=1590, Yt=73, Yc=185, S12=90, S23=55,
        G_Ic=0.28, G_IIc=0.79, rho=1.57e-6,
    ),
    "T700/2510": OrthotropicMaterial(
        name="T700/2510",
        E11=127000, E22=8400, nu12=0.31, G12=4200, G13=4200, G23=2500,
        Xt=2100, Xc=1200, Yt=58, Yc=175, S12=66, S23=45,
        G_Ic=0.22, G_IIc=0.60, rho=1.55e-6,
    ),
    "T800/epoxy": OrthotropicMaterial(
        name="T800/epoxy",
        E11=155000, E22=8500, nu12=0.33, G12=5000, G13=5000, G23=2600,
        Xt=2700, Xc=1680, Yt=68, Yc=200, S12=85, S23=52,
        G_Ic=0.25, G_IIc=0.70, rho=1.58e-6,
    ),
}
```

- [ ] **Step 4: Run tests — expect PASS**

```
pytest tests/core/test_material.py -v
```

- [ ] **Step 5: Commit**

```bash
git add src/bvidfe/core/ tests/core/
git commit -m "Add OrthotropicMaterial dataclass and 4-preset material library"
```

### Task 2.2: `core/laminate.py` — CLT ABD

**Files:**
- Create: `src/bvidfe/core/laminate.py`
- Create: `tests/core/test_laminate.py`

- [ ] **Step 1: Failing tests for CLT**

`tests/core/test_laminate.py`:
```python
import numpy as np
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.core.laminate import Laminate

def test_symmetric_laminate_has_zero_B():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0, 45, -45, 90, 90, -45, 45, 0], ply_thickness_mm=0.152)
    A, B, D = lam.abd_matrices()
    assert np.allclose(B, 0.0, atol=1e-6)

def test_quasi_isotropic_A_is_isotropic():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0, 45, -45, 90]*2, ply_thickness_mm=0.125)
    A, _, _ = lam.abd_matrices()
    # A11 == A22, A16 == A26 == 0 for quasi-iso
    assert abs(A[0,0] - A[1,1]) / A[0,0] < 0.02
    assert abs(A[0,2]) / A[0,0] < 0.02
    assert abs(A[1,2]) / A[0,0] < 0.02

def test_effective_Ex_matches_ply_when_all_zero():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0]*8, ply_thickness_mm=0.125)
    Ex, Ey, Gxy, nuxy = lam.effective_engineering_constants()
    assert abs(Ex - m.E11) / m.E11 < 0.01
```

- [ ] **Step 2: Run — expect failure**

- [ ] **Step 3: Implement Laminate with CLT ABD**

Port from WrinkleFE `src/wrinklefe/core/laminate.py`. Key methods:
- `_Qbar(angle_deg) -> (3,3)` — in-plane reduced stiffness rotated to laminate frame
- `abd_matrices() -> (A, B, D)` — assembled over plies via z-offsets
- `effective_engineering_constants() -> (Ex, Ey, Gxy, nuxy)` — from `A^-1`
- `flexural_rigidity_Deff() -> float` — geometric mean `sqrt(D11 * D22)` (used by Olsson)

Store `layup_deg: list[float]`, `ply_thickness_mm: float`, `material: OrthotropicMaterial`. Expose `thickness_mm` property.

- [ ] **Step 4: Run — expect PASS**

- [ ] **Step 5: Commit**

```bash
git add src/bvidfe/core/laminate.py tests/core/test_laminate.py
git commit -m "Add Laminate with CLT ABD matrices and effective engineering constants"
```

### Task 2.3: `core/geometry.py` — PanelGeometry, ImpactorGeometry

**Files:**
- Create: `src/bvidfe/core/geometry.py`
- Create: `tests/core/test_geometry.py`

- [ ] **Step 1: Tests**

```python
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry

def test_panel_default_boundary_simply_supported():
    p = PanelGeometry(Lx_mm=150, Ly_mm=100)
    assert p.boundary == "simply_supported"

def test_panel_rejects_nonpositive_dims():
    import pytest
    with pytest.raises(ValueError):
        PanelGeometry(Lx_mm=0, Ly_mm=100)

def test_impactor_default_is_hemispherical_16mm():
    i = ImpactorGeometry()
    assert i.diameter_mm == 16.0
    assert i.shape == "hemispherical"
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement two frozen dataclasses with `__post_init__` validation**

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add PanelGeometry and ImpactorGeometry dataclasses"
```

---

## Phase 3 — Damage Model

### Task 3.1: `damage/state.py`

**Files:**
- Create: `src/bvidfe/damage/__init__.py`
- Create: `src/bvidfe/damage/state.py`
- Create: `tests/damage/__init__.py`
- Create: `tests/damage/test_state.py`

- [ ] **Step 1: Tests covering ellipse area, DPA union, per-interface map**

```python
import math
from bvidfe.damage.state import DelaminationEllipse, DamageState

def test_ellipse_area():
    e = DelaminationEllipse(interface_index=0, centroid_mm=(0,0),
                            major_mm=10, minor_mm=5, orientation_deg=0)
    assert math.isclose(e.area_mm2, math.pi * 10 * 5)

def test_damage_state_projected_area_two_identical_ellipses_union():
    e = DelaminationEllipse(0, (0,0), 10, 5, 0)
    e2 = DelaminationEllipse(1, (0,0), 10, 5, 0)   # same footprint, different interface
    ds = DamageState([e, e2], dent_depth_mm=0.3)
    # Union should equal single-ellipse area (they overlap 100% in plan view)
    assert math.isclose(ds.projected_damage_area_mm2, math.pi * 10 * 5, rel_tol=1e-2)

def test_damage_state_projected_area_two_non_overlapping():
    e1 = DelaminationEllipse(0, (0, 0), 10, 5, 0)
    e2 = DelaminationEllipse(1, (100, 100), 10, 5, 0)   # far apart
    ds = DamageState([e1, e2], dent_depth_mm=0.3)
    assert math.isclose(ds.projected_damage_area_mm2, 2 * math.pi * 10 * 5, rel_tol=1e-2)

def test_damage_state_per_interface_area():
    e1 = DelaminationEllipse(0, (0,0), 10, 5, 0)
    e2 = DelaminationEllipse(2, (0,0), 8, 4, 0)
    ds = DamageState([e1, e2], dent_depth_mm=0.3)
    d = ds.per_interface_area
    assert math.isclose(d[0], math.pi*10*5)
    assert math.isclose(d[2], math.pi*8*4)
    assert 1 not in d
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement**

```python
from dataclasses import dataclass, field
import math
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
# NOTE: adds shapely to dependencies; add to pyproject.toml in this task.

def _ellipse_polygon(e: "DelaminationEllipse", n_pts: int = 72) -> Polygon:
    theta = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    c, s = np.cos(np.deg2rad(e.orientation_deg)), np.sin(np.deg2rad(e.orientation_deg))
    x = e.major_mm * np.cos(theta)
    y = e.minor_mm * np.sin(theta)
    xr = c*x - s*y + e.centroid_mm[0]
    yr = s*x + c*y + e.centroid_mm[1]
    return Polygon(zip(xr, yr))

@dataclass
class DelaminationEllipse:
    interface_index: int
    centroid_mm: tuple[float, float]
    major_mm: float
    minor_mm: float
    orientation_deg: float
    def __post_init__(self):
        if self.major_mm <= 0 or self.minor_mm <= 0:
            raise ValueError("ellipse semi-axes must be positive")
        if self.interface_index < 0:
            raise ValueError("interface_index must be >= 0")
    @property
    def area_mm2(self) -> float:
        return math.pi * self.major_mm * self.minor_mm

@dataclass
class DamageState:
    delaminations: list[DelaminationEllipse] = field(default_factory=list)
    dent_depth_mm: float = 0.0
    fiber_break_radius_mm: float = 0.0

    @property
    def projected_damage_area_mm2(self) -> float:
        if not self.delaminations:
            return 0.0
        polys = [_ellipse_polygon(e) for e in self.delaminations]
        return unary_union(polys).area

    @property
    def per_interface_area(self) -> dict[int, float]:
        out: dict[int, float] = {}
        for e in self.delaminations:
            out[e.interface_index] = out.get(e.interface_index, 0.0) + e.area_mm2
        return out
```

Add `shapely>=2.0` to `pyproject.toml` dependencies and `requirements.txt`. Reinstall.

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml requirements.txt src/bvidfe/damage/ tests/damage/
git commit -m "Add DamageState and DelaminationEllipse with union-based DPA"
```

### Task 3.2: `damage/io.py` — C-scan JSON import/export

**Files:**
- Create: `src/bvidfe/damage/io.py`
- Create: `tests/damage/test_io.py`
- Create: `docs/cscan_schema.md`

- [ ] **Step 1: Tests — round-trip + validation errors**

```python
import json, pytest, tempfile, pathlib
from bvidfe.damage.io import (
    damage_state_to_dict, damage_state_from_dict,
    load_cscan_json, save_cscan_json, CScanSchemaError,
)
from bvidfe.damage.state import DamageState, DelaminationEllipse

def _make_state():
    return DamageState([
        DelaminationEllipse(3, (75, 50), 28, 18, 45),
        DelaminationEllipse(4, (78, 52), 32, 20, 50),
    ], dent_depth_mm=0.45, fiber_break_radius_mm=3.0)

def test_round_trip_to_dict_and_back():
    ds = _make_state()
    ds2 = damage_state_from_dict(damage_state_to_dict(ds))
    assert ds2.dent_depth_mm == ds.dent_depth_mm
    assert len(ds2.delaminations) == 2
    assert ds2.delaminations[0].interface_index == 3

def test_round_trip_file(tmp_path):
    ds = _make_state()
    fp = tmp_path / "cscan.json"
    save_cscan_json(ds, fp)
    ds2 = load_cscan_json(fp)
    assert ds2.delaminations[1].major_mm == 32

def test_rejects_bad_schema_version(tmp_path):
    fp = tmp_path / "bad.json"
    fp.write_text(json.dumps({"schema_version": "99.0", "delaminations": [], "dent_depth_mm": 0}))
    with pytest.raises(CScanSchemaError):
        load_cscan_json(fp)

def test_rejects_negative_ellipse(tmp_path):
    fp = tmp_path / "bad.json"
    fp.write_text(json.dumps({
        "schema_version": "1.0", "dent_depth_mm": 0.0,
        "delaminations": [{"interface_index": 0, "centroid_mm":[0,0],
                           "major_mm": -1, "minor_mm": 5, "orientation_deg": 0}]
    }))
    with pytest.raises(CScanSchemaError):
        load_cscan_json(fp)
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement** — functions `damage_state_to_dict`, `damage_state_from_dict`, `save_cscan_json`, `load_cscan_json`, custom `CScanSchemaError`. Schema version `"1.0"`. Validation wraps exceptions into `CScanSchemaError`.

- [ ] **Step 4: Write `docs/cscan_schema.md`** documenting the format with an example (same JSON as the spec).

- [ ] **Step 5: Run tests — pass**

- [ ] **Step 6: Commit**

```bash
git add src/bvidfe/damage/io.py tests/damage/test_io.py docs/cscan_schema.md
git commit -m "Add C-scan JSON I/O with schema validation"
```

---

## Phase 4 — Impact Mapping

### Task 4.1: `impact/olsson.py` — threshold load + onset energy

**Files:**
- Create: `src/bvidfe/impact/__init__.py`
- Create: `src/bvidfe/impact/olsson.py`
- Create: `tests/impact/__init__.py`
- Create: `tests/impact/test_olsson.py`

- [ ] **Step 1: Tests**

```python
import math, pytest
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.core.laminate import Laminate
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry
from bvidfe.impact.olsson import threshold_load, onset_energy, NAVIER_N

def test_threshold_load_scales_with_sqrt_G_IIc():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0,45,-45,90]*4, 0.152)
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    Pc1 = threshold_load(lam, pan, imp)
    m2 = MATERIAL_LIBRARY["IM7/8552"]
    # Clone with 4x G_IIc via dataclasses.replace
    import dataclasses
    m4 = dataclasses.replace(m2, G_IIc=m2.G_IIc * 4)
    lam2 = Laminate(m4, [0,45,-45,90]*4, 0.152)
    Pc2 = threshold_load(lam2, pan, imp)
    assert math.isclose(Pc2 / Pc1, 2.0, rel_tol=0.05)

def test_onset_energy_positive_and_monotonic_in_Pc():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam_thin = Laminate(m, [0,90]*4, 0.125)
    lam_thick = Laminate(m, [0,90]*12, 0.125)
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    E_thin = onset_energy(lam_thin, pan, imp)
    E_thick = onset_energy(lam_thick, pan, imp)
    assert E_thin > 0 and E_thick > 0
    assert E_thick > E_thin      # thicker plate harder to damage

def test_navier_n_is_11_by_default():
    assert NAVIER_N == 11
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement**

`src/bvidfe/impact/olsson.py`:
```python
"""Olsson quasi-static impact threshold load and onset energy.

Pc   = pi * sqrt(8 * G_IIc * D_eff / 9)
Eons = Pc^2 / (2 * k_cb)
k_cb = 1 / (1/k_bending + 1/k_contact)
"""
import math
import numpy as np
from bvidfe.core.laminate import Laminate
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry

NAVIER_N: int = 11   # Navier series truncation (N x N modes)

def _k_bending_ssss(lam: Laminate, pan: PanelGeometry, x0: float, y0: float,
                   n_modes: int = NAVIER_N) -> float:
    """Navier series point-load stiffness of a simply-supported rectangular
    orthotropic plate at (x0, y0)."""
    _, _, D = lam.abd_matrices()
    D11, D22, D12, D66 = D[0,0], D[1,1], D[0,1], D[2,2]
    a, b = pan.Lx_mm, pan.Ly_mm
    w_over_P = 0.0
    for m in range(1, n_modes+1):
        for n in range(1, n_modes+1):
            sin_mx = math.sin(m*math.pi*x0/a)
            sin_ny = math.sin(n*math.pi*y0/b)
            Dmn = (D11*(m*math.pi/a)**4 + 2*(D12+2*D66)*(m*math.pi/a)**2*(n*math.pi/b)**2
                   + D22*(n*math.pi/b)**4)
            w_over_P += (sin_mx*sin_ny)**2 / Dmn
    w_over_P *= 4 / (a*b)
    return 1.0 / w_over_P

def _k_contact_hertz(lam: Laminate, imp: ImpactorGeometry) -> float:
    """Hertzian contact stiffness: k = (4/3) * sqrt(R) * E_eff.
    E_eff combines impactor (assumed steel, E=200 GPa) and transverse laminate modulus."""
    R = imp.diameter_mm / 2.0     # mm
    E_steel = 200e3               # MPa = N/mm^2
    nu_steel = 0.3
    E_plate = lam.material.E22
    nu_plate = 0.3
    inv_E = (1 - nu_steel**2)/E_steel + (1 - nu_plate**2)/E_plate
    E_eff = 1.0 / inv_E
    return (4.0/3.0) * math.sqrt(R) * E_eff     # N/mm^(3/2), linearized below

def threshold_load(lam: Laminate, pan: PanelGeometry, imp: ImpactorGeometry) -> float:
    """Olsson threshold load Pc (N). Uses geometric-mean flexural rigidity D_eff."""
    _, _, D = lam.abd_matrices()
    D_eff = math.sqrt(D[0,0] * D[1,1])   # N*mm
    G_IIc = lam.material.G_IIc           # N/mm
    return math.pi * math.sqrt(8 * G_IIc * D_eff / 9.0)

def onset_energy(lam: Laminate, pan: PanelGeometry, imp: ImpactorGeometry,
                location_xy_mm: tuple[float, float] | None = None) -> float:
    """Impact energy (J) at which damage onsets. Uses Pc and the series stiffness
    k_cb combining plate bending + Hertzian contact at the impact location."""
    if location_xy_mm is None:
        location_xy_mm = (pan.Lx_mm/2, pan.Ly_mm/2)
    x0, y0 = location_xy_mm
    k_b = _k_bending_ssss(lam, pan, x0, y0, n_modes=NAVIER_N)
    k_c = _k_contact_hertz(lam, imp)
    k_cb = 1.0 / (1.0/k_b + 1.0/k_c)
    Pc = threshold_load(lam, pan, imp)
    # Energy in mJ (N*mm), convert to J
    E_mJ = Pc**2 / (2.0 * k_cb)
    return E_mJ * 1e-3
```

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add Olsson quasi-static threshold load and onset energy"
```

### Task 4.2: `impact/shape_templates.py` — peanut template + DPA conservation

**Files:**
- Create: `src/bvidfe/impact/shape_templates.py`
- Create: `tests/impact/test_shape_templates.py`

- [ ] **Step 1: Tests**

```python
import math, pytest
from bvidfe.impact.shape_templates import distribute_damage
from bvidfe.damage.state import DamageState

def test_one_ellipse_per_interface():
    layup = [0, 45, -45, 90, 0, 90, -45, 45, 0]   # 9 plies => 8 interfaces
    ellipses = distribute_damage(layup_deg=layup, target_dpa_mm2=1000.0,
                                 dent_depth_mm=0.4, fiber_break_radius_mm=0.0)
    assert len({e.interface_index for e in ellipses}) == 8

def test_dpa_conservation_within_tolerance():
    layup = [0, 45, -45, 90, 90, -45, 45, 0]      # 8 plies => 7 interfaces
    target = 800.0
    ellipses = distribute_damage(layup_deg=layup, target_dpa_mm2=target,
                                 dent_depth_mm=0.3, fiber_break_radius_mm=0.0)
    ds = DamageState(ellipses, 0.3, 0.0)
    assert abs(ds.projected_damage_area_mm2 - target) / target < 0.01

def test_aspect_ratio_grows_with_ply_angle_mismatch():
    # Single 0/0 interface vs 0/90 interface — 0/90 should yield higher AR
    ellipses_aligned = distribute_damage([0, 0, 0], 400.0, 0.3, 0.0)
    ellipses_cross   = distribute_damage([0, 90, 0], 400.0, 0.3, 0.0)
    def _ar(es): return max(e.major_mm/e.minor_mm for e in es)
    assert _ar(ellipses_cross) > _ar(ellipses_aligned)
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement**

`src/bvidfe/impact/shape_templates.py`:
```python
"""Layup-dependent 'peanut' template distributing DPA into per-interface ellipses."""
import math
from typing import List
from scipy.optimize import brentq
from bvidfe.damage.state import DelaminationEllipse, DamageState

def _aspect_ratio(delta_theta_deg: float) -> float:
    """Ellipse AR (major/minor) as a function of neighbor-ply angle mismatch.
    Linear fit: AR = 1 + 0.025 * |delta_theta|, clipped to [1, 4]."""
    return min(4.0, max(1.0, 1.0 + 0.025 * abs(delta_theta_deg)))

def _orientation_deg(lower_ply_deg: float, upper_ply_deg: float) -> float:
    """Ellipse major-axis orientation = bisector of the two neighbor-ply angles."""
    return 0.5 * (lower_ply_deg + upper_ply_deg)

def _relative_size(interface_index: int, n_interfaces: int) -> float:
    """Through-thickness weighting: largest near back face (i close to n_interfaces-1)."""
    z = (interface_index + 1) / n_interfaces          # 0..1, back face = 1
    return 0.3 + 0.7 * z                              # monotonic increase

def distribute_damage(layup_deg: List[float], target_dpa_mm2: float,
                      dent_depth_mm: float, fiber_break_radius_mm: float,
                      centroid_mm: tuple[float, float] = (0.0, 0.0)) -> list[DelaminationEllipse]:
    n_plies = len(layup_deg)
    n_interfaces = n_plies - 1
    if n_interfaces <= 0 or target_dpa_mm2 <= 0:
        return []

    templates = []
    for i in range(n_interfaces):
        dtheta = layup_deg[i+1] - layup_deg[i]
        ar = _aspect_ratio(dtheta)
        orient = _orientation_deg(layup_deg[i], layup_deg[i+1])
        rel = _relative_size(i, n_interfaces)
        # At scalar=1, ellipse has minor_mm=rel, major_mm=ar*rel  → area = pi*ar*rel^2
        templates.append({"i": i, "ar": ar, "orient": orient, "rel": rel})

    def _union_area(scalar: float) -> float:
        es = [DelaminationEllipse(
                  interface_index=t["i"],
                  centroid_mm=centroid_mm,
                  major_mm=scalar * t["ar"] * t["rel"],
                  minor_mm=scalar * t["rel"],
                  orientation_deg=t["orient"],
              ) for t in templates]
        return DamageState(es, dent_depth_mm, fiber_break_radius_mm).projected_damage_area_mm2

    # Bracket: at s=0.1 union tiny, at s=50 union >> target for any sane target.
    # Expand upper bracket if needed.
    lo, hi = 0.1, 50.0
    while _union_area(hi) < target_dpa_mm2:
        hi *= 2
        if hi > 1e6:
            raise RuntimeError("shape_templates: cannot bracket target DPA")
    scalar = brentq(lambda s: _union_area(s) - target_dpa_mm2, lo, hi, xtol=1e-3)

    return [DelaminationEllipse(
                interface_index=t["i"],
                centroid_mm=centroid_mm,
                major_mm=scalar * t["ar"] * t["rel"],
                minor_mm=scalar * t["rel"],
                orientation_deg=t["orient"],
            ) for t in templates]
```

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add peanut-template DPA distribution with Brent scaling for conservation"
```

### Task 4.3: `impact/dent_model.py`

**Files:**
- Create: `src/bvidfe/impact/dent_model.py`
- Create: `tests/impact/test_dent_model.py`

- [ ] **Step 1: Tests**

```python
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.dent_model import dent_depth_mm, fiber_break_radius_mm

def test_zero_below_threshold():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert dent_depth_mm(m, E_impact_J=1.0, E_onset_J=5.0, h_mm=2.0) == 0.0

def test_monotonic_in_energy():
    m = MATERIAL_LIBRARY["IM7/8552"]
    d1 = dent_depth_mm(m, E_impact_J=10.0, E_onset_J=5.0, h_mm=2.0)
    d2 = dent_depth_mm(m, E_impact_J=30.0, E_onset_J=5.0, h_mm=2.0)
    assert d2 > d1 > 0

def test_fiber_break_radius_zero_when_eta_zero():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert m.fiber_break_eta == 0.0
    assert fiber_break_radius_mm(m, E_impact_J=100.0) == 0.0
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement**

```python
import math
from bvidfe.core.material import OrthotropicMaterial

def dent_depth_mm(m: OrthotropicMaterial, E_impact_J: float, E_onset_J: float, h_mm: float) -> float:
    if E_impact_J <= E_onset_J:
        return 0.0
    # d/h = beta * ((E - E_onset) / (G_Ic * h^2))^gamma
    # G_Ic in N/mm; E in J = N*m = 1e3 N*mm
    numerator = (E_impact_J - E_onset_J) * 1e3
    denom = m.G_Ic * h_mm**2
    return h_mm * m.dent_beta * (numerator / denom) ** m.dent_gamma

def fiber_break_radius_mm(m: OrthotropicMaterial, E_impact_J: float) -> float:
    excess = max(0.0, E_impact_J - m.fiber_break_E_threshold)
    return m.fiber_break_eta * math.sqrt(excess)
```

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add dent-depth and fiber-break-radius empirical models"
```

### Task 4.4: `impact/mapping.py` — orchestrator

**Files:**
- Create: `src/bvidfe/impact/mapping.py`
- Create: `tests/impact/test_mapping.py`

- [ ] **Step 1: Tests**

```python
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.core.laminate import Laminate
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage

def test_below_threshold_returns_empty_damage():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0,45,-45,90]*4, 0.152)
    pan = PanelGeometry(150, 100)
    ev = ImpactEvent(energy_J=0.01, impactor=ImpactorGeometry(), mass_kg=5.5)
    ds = impact_to_damage(ev, lam, pan)
    assert ds.delaminations == []
    assert ds.dent_depth_mm == 0.0

def test_above_threshold_produces_n_minus_1_ellipses():
    m = MATERIAL_LIBRARY["IM7/8552"]
    layup = [0,45,-45,90]*4
    lam = Laminate(m, layup, 0.152)
    pan = PanelGeometry(150, 100)
    ev = ImpactEvent(energy_J=30.0, impactor=ImpactorGeometry(), mass_kg=5.5)
    ds = impact_to_damage(ev, lam, pan)
    assert len({e.interface_index for e in ds.delaminations}) == len(layup) - 1
    assert ds.dent_depth_mm > 0
    assert ds.projected_damage_area_mm2 > 0
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement `ImpactEvent` dataclass + `impact_to_damage(event, lam, pan) -> DamageState`**

Pipeline inside `impact_to_damage`:
1. Compute `E_onset` via `onset_energy(...)`
2. If `event.energy_J <= E_onset`, return empty `DamageState`
3. Compute `DPA_target = material.olsson_alpha * (E_impact - E_onset) * 1e3 / (G_IIc * h_total)` (mm²). Guard negative.
4. Call `distribute_damage(layup, DPA_target, dent_depth, fiber_break_radius, centroid=event.location_xy_mm)`
5. Compute dent and fiber-break radius via phase 4.3 helpers
6. Return `DamageState(ellipses, dent_depth_mm, fiber_break_radius_mm)`

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add impact_to_damage orchestrator tying Olsson + templates + dent"
```

---

## Phase 5 — Failure Criteria

### Task 5.1: `failure/tsai_wu.py`

**Files:**
- Create: `src/bvidfe/failure/__init__.py`
- Create: `src/bvidfe/failure/tsai_wu.py`
- Create: `tests/failure/__init__.py`
- Create: `tests/failure/test_tsai_wu.py`

- [ ] **Step 1: Tests — pristine strength and simple uniaxial recovery**

```python
import math
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.tsai_wu import tsai_wu_index, tsai_wu_strength_uniaxial

def test_uniaxial_Xt_gives_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    # Apply sigma_1 = Xt, index should be ~1
    s = [m.Xt, 0, 0, 0, 0, 0]
    assert abs(tsai_wu_index(m, s) - 1.0) < 0.02

def test_zero_stress_gives_zero_index():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert tsai_wu_index(m, [0]*6) == 0

def test_strength_uniaxial_tension_matches_Xt_within_5pct():
    m = MATERIAL_LIBRARY["IM7/8552"]
    sig = tsai_wu_strength_uniaxial(m, direction=1, sign=+1)
    assert 0.95*m.Xt <= sig <= 1.05*m.Xt
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement 3D Tsai-Wu** — 6-vector Voigt stress. Coefficients F1..F6 and F11..F66 per classical formulation; F12 uses the Tsai-Hahn coupling.

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add 3D Tsai-Wu failure criterion"
```

### Task 5.2: `failure/larc05.py`

**Files:**
- Create: `src/bvidfe/failure/larc05.py`
- Create: `tests/failure/test_larc05.py`

- [ ] **Step 1: Tests — port validation vectors from WrinkleFE `tests/test_failure_larc05.py`**

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement LaRC05** by porting `src/wrinklefe/failure/larc05.py` from WrinkleFE, adjusting types to match `OrthotropicMaterial`.

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Port LaRC05 failure criterion from WrinkleFE"
```

### Task 5.3: `failure/soutis_openhole.py`

**Files:**
- Create: `src/bvidfe/failure/soutis_openhole.py`
- Create: `tests/failure/test_soutis.py`

- [ ] **Step 1: Tests**

```python
import math
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.soutis_openhole import (
    soutis_cai, whitney_nuismer_tai,
)

def test_soutis_monotone_decreasing_in_dpa():
    m = MATERIAL_LIBRARY["IM7/8552"]
    A_panel = 150*100
    s1 = soutis_cai(m, dpa_mm2=100, A_panel_mm2=A_panel, sigma_pristine_MPa=500)
    s2 = soutis_cai(m, dpa_mm2=500, A_panel_mm2=A_panel, sigma_pristine_MPa=500)
    assert s1 > s2

def test_soutis_zero_dpa_returns_pristine():
    m = MATERIAL_LIBRARY["IM7/8552"]
    s = soutis_cai(m, dpa_mm2=0, A_panel_mm2=15000, sigma_pristine_MPa=500)
    assert math.isclose(s, 500, rel_tol=1e-6)

def test_wn_tai_monotone_decreasing_in_dpa():
    m = MATERIAL_LIBRARY["IM7/8552"]
    t1 = whitney_nuismer_tai(m, dpa_mm2=100, sigma_pristine_MPa=800)
    t2 = whitney_nuismer_tai(m, dpa_mm2=400, sigma_pristine_MPa=800)
    assert t1 > t2
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement**

```python
import math
from bvidfe.core.material import OrthotropicMaterial

def soutis_cai(m, dpa_mm2, A_panel_mm2, sigma_pristine_MPa):
    """σ_CAI / σ_0 = 1 / (1 + k_s * (DPA/A_panel)^m)"""
    if dpa_mm2 <= 0:
        return sigma_pristine_MPa
    kd = 1.0 / (1.0 + m.soutis_k_s * (dpa_mm2 / A_panel_mm2) ** m.soutis_m)
    return kd * sigma_pristine_MPa

def whitney_nuismer_tai(m, dpa_mm2, sigma_pristine_MPa):
    """Point-stress on equivalent circular hole of diameter 2*sqrt(DPA/π)."""
    if dpa_mm2 <= 0:
        return sigma_pristine_MPa
    R = math.sqrt(dpa_mm2 / math.pi)
    d0 = m.wn_d0_mm
    xi = R / (R + d0)
    Kt_factor = (1 + 0.5*xi**2 + 1.5*xi**4 - (3 - 1)*(5*xi**6 - 7*xi**8))   # Whitney-Nuismer
    return sigma_pristine_MPa / Kt_factor
```

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add Soutis CAI and Whitney-Nuismer TAI open-hole-equivalent models"
```

### Task 5.4: `failure/evaluator.py`

**Files:**
- Create: `src/bvidfe/failure/evaluator.py`
- Create: `tests/failure/test_evaluator.py`

- [ ] **Step 1: Tests** — `FailureEvaluator` applies a chosen criterion across a stress field and returns max index + critical element/gauss-point.

- [ ] **Step 2–5:** Implement, port patterns from WrinkleFE `failure/evaluator.py`, test, commit.

```bash
git commit -am "Add FailureEvaluator with criterion dispatch"
```

---

## Phase 6 — Empirical Tier Wired End-to-End

### Task 6.1: `analysis/config.py` + `analysis/bvid.py` (empirical path only)

**Files:**
- Create: `src/bvidfe/analysis/__init__.py`
- Create: `src/bvidfe/analysis/config.py`
- Create: `src/bvidfe/analysis/results.py`
- Create: `src/bvidfe/analysis/bvid.py`
- Create: `tests/analysis/__init__.py`
- Create: `tests/analysis/test_empirical_path.py`

- [ ] **Step 1: Tests for the empirical end-to-end path** — both workflow paths (impact + damage) yield sensible CAI and TAI.

```python
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.analysis import AnalysisConfig, BvidAnalysis

def test_empirical_cai_impact_path():
    cfg = AnalysisConfig(
        material="IM7/8552", layup_deg=[0,45,-45,90]*4, ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5),
        loading="compression", tier="empirical",
    )
    r = BvidAnalysis(cfg).run()
    assert 0.3 < r.knockdown < 1.0
    assert r.residual_strength_MPa < r.pristine_strength_MPa
    assert r.damage.projected_damage_area_mm2 > 0
    assert r.tier_used == "empirical"

def test_empirical_damage_path_equivalent_to_impact_path():
    # Run impact path, capture damage, re-run with damage path, assert equal KD.
    cfg_imp = AnalysisConfig(
        material="IM7/8552", layup_deg=[0,45,-45,90]*4, ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5),
        loading="compression", tier="empirical",
    )
    r_imp = BvidAnalysis(cfg_imp).run()
    cfg_dmg = AnalysisConfig(
        material="IM7/8552", layup_deg=[0,45,-45,90]*4, ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        damage=r_imp.damage,
        loading="compression", tier="empirical",
    )
    r_dmg = BvidAnalysis(cfg_dmg).run()
    assert abs(r_imp.knockdown - r_dmg.knockdown) < 1e-6

def test_tension_path_runs():
    cfg = AnalysisConfig(
        material="IM7/8552", layup_deg=[0,45,-45,90]*4, ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        damage=DamageState([DelaminationEllipse(3,(75,50),20,12,45)], 0.4),
        loading="tension", tier="empirical",
    )
    r = BvidAnalysis(cfg).run()
    assert r.knockdown < 1.0
```

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement**

`src/bvidfe/analysis/config.py` — `AnalysisConfig` + `MeshParams` per spec Section 3.1.

`src/bvidfe/analysis/results.py` — `AnalysisResults` + `FieldResults` per spec. `FieldResults` fields optional/None now; populated by phase 10.

`src/bvidfe/analysis/bvid.py`:
```python
from copy import deepcopy
from dataclasses import asdict
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.core.laminate import Laminate
from bvidfe.impact.mapping import impact_to_damage
from bvidfe.failure.soutis_openhole import soutis_cai, whitney_nuismer_tai
from .config import AnalysisConfig
from .results import AnalysisResults

def _pristine_strength(lam: Laminate, loading: str) -> float:
    """Rough laminate-level pristine strength estimate via rule-of-mixtures on fiber-aligned plies."""
    # compression: min over plies of Xc*cos^2 + Yc*sin^2  (CLT-weighted)
    # tension: analogous with Xt, Yt
    import math
    m = lam.material
    total_t = 0.0
    num = 0.0
    for th in lam.layup_deg:
        c2, s2 = math.cos(math.radians(th))**2, math.sin(math.radians(th))**2
        if loading == "compression":
            num += lam.ply_thickness_mm * (m.Xc*c2 + m.Yc*s2)
        else:
            num += lam.ply_thickness_mm * (m.Xt*c2 + m.Yt*s2)
        total_t += lam.ply_thickness_mm
    return num / total_t

class BvidAnalysis:
    def __init__(self, config: AnalysisConfig):
        self.config = config

    def run(self) -> AnalysisResults:
        mat = (self.config.material if not isinstance(self.config.material, str)
               else MATERIAL_LIBRARY[self.config.material])
        lam = Laminate(mat, self.config.layup_deg, self.config.ply_thickness_mm)
        damage = self._resolve_damage(lam)
        sigma_0 = _pristine_strength(lam, self.config.loading)
        if self.config.tier == "empirical":
            sigma = self._empirical(lam, damage, sigma_0)
            field_results = None
            eigs = None
            critical_iface = None
        else:
            # wired in phases 9 and 10
            raise NotImplementedError(f"tier '{self.config.tier}' not yet wired")

        return AnalysisResults(
            residual_strength_MPa=sigma,
            pristine_strength_MPa=sigma_0,
            knockdown=sigma / sigma_0,
            damage=damage,
            dpa_mm2=damage.projected_damage_area_mm2,
            buckling_eigenvalues=eigs,
            critical_sublaminate=critical_iface,
            field_results=field_results,
            tier_used=self.config.tier,
            config_snapshot=deepcopy(asdict(self.config)),
        )

    def _resolve_damage(self, lam):
        if self.config.damage is not None:
            return self.config.damage
        return impact_to_damage(self.config.impact, lam, self.config.panel)

    def _empirical(self, lam, damage, sigma_0):
        A_panel = self.config.panel.Lx_mm * self.config.panel.Ly_mm
        dpa = damage.projected_damage_area_mm2
        if self.config.loading == "compression":
            return soutis_cai(lam.material, dpa, A_panel, sigma_0)
        return whitney_nuismer_tai(lam.material, dpa, sigma_0)
```

Expose `AnalysisConfig`, `BvidAnalysis`, `AnalysisResults`, `MeshParams` from `src/bvidfe/analysis/__init__.py`.

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Wire empirical tier end-to-end: BvidAnalysis.run() for CAI and TAI"
```

### Task 6.2: Smoke-test the CLI

**Files:**
- Create: `src/bvidfe/cli.py`
- Create: `tests/test_cli.py`

- [ ] **Step 1: Tests** — `bvidfe --tier empirical --energy 30 --material IM7/8552 --layup 0,45,-45,90,0,45,-45,90 --thickness 0.152 --panel 150x100 --loading compression` prints JSON with `knockdown`.

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement** using `argparse`. Parse layup from comma-separated, panel from `<Lx>x<Ly>`, print `json.dumps(result.to_dict(), indent=2)`.

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add bvidfe CLI entry point for empirical tier"
```

---

## Phase 7 — Elements

### Task 7.1: `elements/gauss.py`

**Files:**
- Create: `src/bvidfe/elements/__init__.py`
- Create: `src/bvidfe/elements/gauss.py`
- Create: `tests/elements/__init__.py`
- Create: `tests/elements/test_gauss.py`

- [ ] **Step 1: Test** — Gauss points integrate `∫_{-1}^{1} x^2 dx = 2/3` exactly for n=2.

- [ ] **Step 2–5:** Port WrinkleFE `elements/gauss.py`; run; commit.

```bash
git commit -am "Port Gauss-Legendre quadrature utilities"
```

### Task 7.2: `elements/hex8.py` (port)

**Files:**
- Create: `src/bvidfe/elements/hex8.py`
- Create: `tests/elements/test_hex8.py`

- [ ] **Steps 1–5:** Port WrinkleFE `elements/hex8.py` (shape funcs, Jacobian, B matrix, stiffness via 2×2×2 Gauss). Port tests. Validate against WrinkleFE's golden values where applicable. Commit.

```bash
git commit -am "Port Hex8 element with 2x2x2 Gauss stiffness integration"
```

### Task 7.3: `elements/hex8i.py` (port)

Same pattern. Hex8 with incompatible modes for improved bending response. Commit as `Port Hex8i incompatible-modes element`.

### Task 7.4: `elements/cohesive.py` — zero-thickness cohesive surface

**Files:**
- Create: `src/bvidfe/elements/cohesive.py`
- Create: `tests/elements/test_cohesive.py`

- [ ] **Step 1: Tests** — (a) element stiffness for pristine (undamaged) state produces near-rigid connection in normal direction; (b) traction-separation softens after δ_0, reaches zero at δ_f; (c) ERR under simple mode-II shear integrates to G_IIc.

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement** 8-node zero-thickness cohesive surface element (4 nodes on each side of an interface). Bilinear traction-separation law parameterized by `(G_Ic, G_IIc, sigma_n_max, tau_max)`. Expose `stiffness_matrix(u_rel)`, `traction(u_rel)`, and an `is_precracked` flag (when set, traction ≡ 0 for all δ).

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add zero-thickness cohesive surface element with bilinear traction-separation"
```

---

## Phase 8 — Solver

### Task 8.1: `solver/assembler.py`

**Files:**
- Create: `src/bvidfe/solver/__init__.py`
- Create: `src/bvidfe/solver/assembler.py`
- Create: `tests/solver/__init__.py`
- Create: `tests/solver/test_assembler.py`

- [ ] **Step 1: Test** — assembly of K on a 2-element cube of 1-mm hex8s produces a 24-DOF symmetric sparse matrix whose free-free rigid-body null space has 6 zero eigenvalues.

- [ ] **Step 2–5:** Port WrinkleFE `solver/assembler.py`; COO→CSC; test; commit.

```bash
git commit -am "Port sparse K assembler for hex8 meshes"
```

### Task 8.2: `solver/boundary.py`

**Files:**
- Create: `src/bvidfe/solver/boundary.py`
- Create: `tests/solver/test_boundary.py`

Port WrinkleFE patterns for penalty BC application, compression/tension BCs. Commit.

### Task 8.3: `solver/static.py`

**Files:**
- Create: `src/bvidfe/solver/static.py`
- Create: `tests/solver/test_static.py`

Port and adapt WrinkleFE `solver/static.py`. Validates against analytical axial-strain benchmark. Commit.

### Task 8.4: `solver/buckling.py` — linear eigenvalue buckling with `eigsh`

**Files:**
- Create: `src/bvidfe/solver/buckling.py`
- Create: `tests/solver/test_buckling.py`

- [ ] **Step 1: Test** — SSSS isotropic square plate in uniaxial compression: first eigenvalue matches classical `π²·D/(a²·σ_applied·h)` within 5%.

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement**

```python
from scipy.sparse.linalg import eigsh
def linear_buckling(K, Kg, n_modes: int = 3):
    # K φ = λ Kg φ, symmetric generalized, smallest positive λ via shift-invert
    vals, vecs = eigsh(K, M=Kg, k=n_modes, sigma=0.0, which="LM")
    # Sort ascending and discard negative/spurious
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    return vals[order], vecs[:, order]
```

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add linear buckling eigenvalue solver using scipy eigsh with shift-invert"
```

---

## Phase 9 — Semi-Analytical Tier

### Task 9.1: Rayleigh-Ritz sublaminate buckling over an ellipse

**Files:**
- Create: `src/bvidfe/analysis/semi_analytical.py`
- Create: `tests/analysis/test_semi_analytical.py`

- [ ] **Step 1: Tests** — (a) circular-limit reduces to known closed-form for orthotropic plate; (b) sublaminate buckling load decreases monotonically with ellipse area; (c) semi_analytical CAI < pristine and > empirical CAI for small DPAs and sometimes < empirical for large DPAs (acceptance: just check <1 and >0).

- [ ] **Step 2: Run — fail**

- [ ] **Step 3: Implement** `sublaminate_buckling_load(material, layup_slice, ellipse)`:
  1. Build sublaminate from plies *above* the critical interface (i.e., between outer surface and the delaminated interface).
  2. Compute sublaminate ABD via CLT.
  3. Map ellipse to local (η, ξ) with `η = x/a`, `ξ = y/b`.
  4. Build 5×5 basis of clamped-beam eigenfunction products `φ_mn(η, ξ) = w_m(η) * w_n(ξ)`.
  5. Assemble K (stiffness) and K_g (in-plane unit compressive load) via Gauss quadrature (10×10).
  6. Solve `eigsh(K, M=K_g, k=1, sigma=0)`; return smallest positive eigenvalue as critical buckling load.
  Also implement `find_critical_interface(damage, laminate)` by scoring each interface as `area * |z_centroid|` (maximum through-thickness asymmetry × largest ellipse).

- [ ] **Step 4: Run — pass**

- [ ] **Step 5: Commit**

```bash
git commit -am "Add Rayleigh-Ritz sublaminate buckling over elliptical delamination"
```

### Task 9.2: Post-buckling residual strength (Soutis envelope) + TAI notch

**Files:**
- Modify: `src/bvidfe/analysis/semi_analytical.py`
- Modify: `tests/analysis/test_semi_analytical.py`

- [ ] Test: `semi_analytical_cai` combines the buckling load and Soutis envelope to give σ_CAI; for a pristine case returns σ_0.

- [ ] Implement function `semi_analytical_cai(lam, damage, sigma_0, A_panel)` and `semi_analytical_tai(lam, damage, sigma_0)` (the latter already reusable via Soutis notch: equivalent circular hole + cohesive-zone strength from in-situ transverse strength).

- [ ] Commit: `Add post-buckling CAI envelope and Soutis notch TAI`

### Task 9.3: Wire into `BvidAnalysis`

**Files:**
- Modify: `src/bvidfe/analysis/bvid.py`
- Modify: `tests/analysis/test_empirical_path.py` (extend, or add a new file)

- [ ] Add `elif self.config.tier == "semi_analytical":` branch calling 9.1/9.2.
- [ ] Populate `buckling_eigenvalues`, `critical_sublaminate` in `AnalysisResults`.
- [ ] Tests: round-trip both loading modes through the tier; assert monotonicity vs energy on a parametric sweep.
- [ ] Commit: `Wire semi_analytical tier into BvidAnalysis`

---

## Phase 10 — 3D FE Tier

### Task 10.1: Structured hex mesh with delamination interfaces

**Files:**
- Create: `src/bvidfe/analysis/fe_mesh.py`
- Create: `tests/analysis/test_fe_mesh.py`

- [ ] Test: mesh generator returns `n_nodes = (nx+1)(ny+1)(nz+1)` and `n_elems = nx·ny·nz`; pre-cracked cohesive interface element count equals number of 4-node quads covered by the ellipse footprint.

- [ ] Implement `build_fe_mesh(config, damage) -> FeMesh`. Honors `MeshParams.elements_per_ply`, `MeshParams.in_plane_size_mm`. Assigns ply-index and ply-angle per element. Inserts cohesive surface elements only at interfaces appearing in `damage.delaminations`; marks those inside ellipse footprint as `is_precracked=True`. Assigns near-zero isotropic stiffness (~1 MPa) to elements inside `fiber_break_radius_mm` of the centroid.

- [ ] Commit: `Add FE mesh builder with pre-cracked cohesive surfaces and fiber-break core`

### Task 10.2: 3D FE CAI solve (linear buckling + LaRC05 eval)

**Files:**
- Create: `src/bvidfe/analysis/fe_tier.py`
- Create: `tests/analysis/test_fe_tier.py`

- [ ] Test: for a small 6×4×4 panel with one small delamination, `fe3d_cai` returns a knockdown < 1 and > 0.3; runs in under 60 s on CI; matches semi_analytical CAI within 40% (loose, since the two models are different fidelity).

- [ ] Implement `fe3d_cai(config, damage, lam, sigma_0)`:
  1. Build mesh (task 10.1)
  2. Assemble K under uniaxial pre-stress → K_g via element-by-element geometric stiffness integration
  3. `linear_buckling(K, Kg, n_modes=3)`
  4. Scale mode 1 to unit amplitude; recover stresses from displacement via Hex8 element stress recovery
  5. Evaluate LaRC05 at every Gauss point; σ_CAI = σ_applied × λ_1 × (1 / max_failure_index_of_scaled_state)
  Likewise `fe3d_tai(config, damage, lam, sigma_0)` as a static uniaxial tension solve with Tsai-Wu failure search (bisection on applied strain).

- [ ] Commit: `Add 3D FE CAI and TAI tiers with cohesive pre-cracked interfaces`

### Task 10.3: Wire `fe3d` into `BvidAnalysis` + populate `FieldResults`

**Files:**
- Modify: `src/bvidfe/analysis/bvid.py`
- Modify: `tests/analysis/test_empirical_path.py` (add fe3d cases, possibly marked `@pytest.mark.slow`)

- [ ] Add `elif self.config.tier == "fe3d":` branch calling `fe_tier`. Populate `field_results=FieldResults(...)` with nodal displacement, element stresses, buckling mode.

- [ ] Commit: `Wire fe3d tier into BvidAnalysis and populate FieldResults`

---

## Phase 11 — Visualization

### Task 11.1: `viz/style.py` + `viz/plots_2d.py`

**Files:**
- Create: `src/bvidfe/viz/__init__.py`
- Create: `src/bvidfe/viz/style.py`
- Create: `src/bvidfe/viz/plots_2d.py`
- Create: `tests/viz/__init__.py`
- Create: `tests/viz/test_plots_2d.py`

- [ ] Tests: smoke — each plot function returns a Matplotlib Figure with expected number of axes and saves a PNG larger than 10 KB to `tmp_path`.

- [ ] Implement:
  - `plot_damage_map(damage: DamageState, panel: PanelGeometry) -> Figure` — top-down ellipse overlay, color by `interface_index`, plus side-view dent profile subplot.
  - `plot_knockdown_curve(energies_J, knockdowns, tier_label) -> Figure`
  - `plot_tier_comparison(energies_J, results_per_tier) -> Figure`

- [ ] Commit: `Add 2D visualization (damage map, knockdown curves, tier comparison)`

### Task 11.2: `viz/plots_3d.py`

**Files:**
- Create: `src/bvidfe/viz/plots_3d.py`
- Create: `tests/viz/test_plots_3d.py`

- [ ] Tests: smoke — `plot_mesh_with_delams(field_results, damage)` returns a `pyvista.Plotter` (off-screen), `plot_buckling_mode(field_results)` likewise, `plot_stress_field(field_results, component="s11")` likewise.

- [ ] Implement using pyvista `UnstructuredGrid` from mesh + cohesive surface polygons overlaid.

- [ ] Commit: `Add 3D visualization with PyVista`

---

## Phase 12 — Parametric Sweeps

### Task 12.1: `sweep/parametric_sweep.py`

**Files:**
- Create: `src/bvidfe/sweep/__init__.py`
- Create: `src/bvidfe/sweep/parametric_sweep.py`
- Create: `tests/sweep/__init__.py`
- Create: `tests/sweep/test_sweep.py`

- [ ] Tests: `sweep_energies(cfg, energies_J=[5,10,20,40]) -> pandas.DataFrame` returns one row per energy with `knockdown, residual_MPa, dpa_mm2, dent_mm`.

- [ ] Implement `sweep_energies`, `sweep_layups`, `sweep_thicknesses`; all return DataFrames and write CSV if path given.

- [ ] Commit: `Add parametric sweep utilities with DataFrame + CSV output`

Note: add `pandas` to dependencies in this task.

---

## Phase 13 — GUI

### Task 13.1: GUI scaffolding + material panel

**Files:**
- Create: `src/bvidfe/gui/__init__.py`
- Create: `src/bvidfe/gui/app.py` (entry point)
- Create: `src/bvidfe/gui/main_window.py`
- Create: `src/bvidfe/gui/panels/material_panel.py`
- Create: `tests/gui/__init__.py`
- Create: `tests/gui/test_material_panel.py`

- [ ] Tests (pytest-qt): constructing `MaterialPanel` populates dropdown with 4 presets; selecting a preset emits `configChanged` signal.

- [ ] Implement `QMainWindow` skeleton with empty dock + stub panels. Port threading pattern from WrinkleFE `gui/main_window.py`.

- [ ] Commit: `Scaffold GUI with main window and material panel`

### Task 13.2: Remaining input panels

**Files:**
- Create: `src/bvidfe/gui/panels/panel_panel.py`
- Create: `src/bvidfe/gui/panels/input_mode_panel.py`
- Create: `src/bvidfe/gui/panels/impact_panel.py`
- Create: `src/bvidfe/gui/panels/damage_panel.py`
- Create: `src/bvidfe/gui/panels/analysis_panel.py`
- Create: `src/bvidfe/gui/panels/sweep_panel.py`
- Create: `tests/gui/test_panels_smoke.py`

- [ ] For each panel: smoke test (constructs without error; emits `configChanged` when value edited). Implement. Commit.

**Impact panel** must include: energy, impactor diameter, impactor **shape** dropdown (hemispherical / flat / conical), impactor mass, location (Lx/2, Ly/2 default). Show live-computed `E_onset` in a read-only label (throttle recomputation on a 200ms timer).

**Damage panel** must include: QTableWidget for delaminations (interface, major, minor, orientation, centroid_x, centroid_y), dent-depth spinbox, fiber-break-radius spinbox, "Import C-scan…" button opening a QFileDialog and calling `damage.io.load_cscan_json`.

**Analysis panel** must include: tier dropdown, loading dropdown, mesh-params group (visible only if tier=fe3d), Run button that enqueues an `AnalysisWorker`.

- [ ] Commit after each panel finishes its tests.

### Task 13.3: Result tabs

**Files:**
- Create: `src/bvidfe/gui/tabs/summary_tab.py`
- Create: `src/bvidfe/gui/tabs/damage_map_tab.py`
- Create: `src/bvidfe/gui/tabs/knockdown_tab.py`
- Create: `src/bvidfe/gui/tabs/mesh_tab.py` (fe3d)
- Create: `src/bvidfe/gui/tabs/buckling_tab.py`
- Create: `src/bvidfe/gui/tabs/stress_tab.py` (fe3d)
- Create: `tests/gui/test_tabs.py`

Each tab: QWidget embedding a Matplotlib `FigureCanvasQTAgg` (or `pyvistaqt.QtInteractor`) that accepts an `AnalysisResults` and draws the corresponding `viz.plots_*` output. Smoke-test construction + `update(results)` call.

- [ ] Commit each.

### Task 13.4: `AnalysisWorker` + `SweepWorker`

**Files:**
- Create: `src/bvidfe/gui/workers.py`
- Create: `tests/gui/test_workers.py`

- [ ] `AnalysisWorker(QThread)` emits `progress(int)`, `resultReady(AnalysisResults)`, `error(str)`. Takes `AnalysisConfig` at construction, calls `BvidAnalysis(config).run()` in `run()`.

- [ ] Tests use `qtbot.waitSignal` to assert `resultReady` fires within 5s on an empirical-tier mock.

- [ ] Commit: `Add AnalysisWorker and SweepWorker with Qt signals`

### Task 13.5: Config save/load + results export

**Files:**
- Create: `src/bvidfe/gui/dialogs/export_dialog.py`
- Modify: `src/bvidfe/gui/main_window.py`
- Create: `tests/gui/test_io_menu.py`

- [ ] Implement File menu: Save Config (`AnalysisConfig → JSON`), Load Config, Export Results (PNG/JSON/CSV).

- [ ] Commit: `Add GUI file I/O menu with config and results export`

---

## Phase 14 — Validation

### Task 14.1: Validation harness

**Files:**
- Create: `validation/__init__.py` (empty marker)
- Create: `validation/reference_data.json`
- Create: `validation/validate_bvid_public.py`
- Create: `tests/test_validation_smoke.py`

- [ ] Test: `tests/test_validation_smoke.py` imports the validator module, constructs a single `DatasetCase`, runs the tier dispatcher, asserts a dict with keys `{sigma_pred, sigma_test, error_pct}` is produced.

- [ ] Implement the harness: `DatasetCase` dataclass (material, layup, ply_thickness, panel_size, impact_energy, measured_dent, measured_dpa, measured_cai_MPa). `run_dataset(cases, tier) -> pd.DataFrame`. Compute MAE% on `residual_strength_MPa`.

- [ ] Commit: `Add validation harness with DatasetCase and MAE% metric`

### Task 14.2: Digitize + ingest published datasets

**Files:**
- Create: `validation/datasets/soutis_as4_3501-6.json`
- Create: `validation/datasets/caprino_as4_epoxy.json`
- Create: `validation/datasets/sanchez-saez_im7_8552.json`
- Create: `validation/datasets/nasa_cai_round_robin.json`

For each dataset:
- [ ] Digitize test data from the cited source into the JSON schema (`interface_index` map, energies, CAI/TAI strengths).
- [ ] Manually tune `olsson_alpha`, `soutis_k_s`, `soutis_m`, `wn_d0_mm`, `dent_beta`, `dent_gamma` in the matching `MATERIAL_LIBRARY` entry until the MAE% for the empirical tier is within the spec's target.
- [ ] Commit per dataset: `Calibrate <dataset>: MAE <x>%`

### Task 14.3: Regression gate

**Files:**
- Modify: `.github/workflows/tests.yml` (or create `validation.yml`)

- [ ] Add CI job running `python validation/validate_bvid_public.py --gate` that exits non-zero if any dataset MAE exceeds 1.25× the spec target.

- [ ] Commit: `Add validation regression gate to CI`

---

## Phase 15 — Packaging

### Task 15.1: PyInstaller spec

**Files:**
- Create: `BvidFE.spec`

Adapt from WrinkleFE `wrinklefe.spec` and PorosityFE `PorosityFE.spec`. Bundle pyvista QT plugins per those templates.

- [ ] Build locally: `pyinstaller BvidFE.spec --noconfirm --clean`. Confirm `dist/BvidFE/BvidFE.app` (macOS) or `dist/BvidFE/BvidFE.exe` (Windows) launches.

- [ ] Commit: `Add PyInstaller spec for standalone app bundling`

### Task 15.2: Release workflow

**Files:**
- Create: `.github/workflows/release.yml`

- [ ] On tag push `v*`, build PyInstaller app on macOS + Windows runners, upload as GitHub Release artifacts.

- [ ] Commit: `Add release workflow producing signed-app artifacts on tag push`

---

## Phase 16 — Documentation

### Task 16.1: README

**Files:**
- Create: `README.md`
- Create: `screenshots/` (run GUI, capture PNGs)

- [ ] Include: logo/badges, one-paragraph purpose, Features, Screenshots, Installation (from source + PyPI), Quick Start (GUI, CLI, Python API), Physics Models (summary of Section 4), Validation table (from phase 14 results), Citation BibTeX stub, MIT License line, Contributing link.

- [ ] Commit: `Add README with features, validation, and usage`

### Task 16.2: ARCHITECTURE.md, CLAUDE.md, CONTRIBUTING.md, CHANGELOG.md, CITATION.cff, LICENSE

**Files:**
- Create: `ARCHITECTURE.md` (port structure from WrinkleFE; update module table to match this plan)
- Create: `CLAUDE.md` (contributor / future-Claude guidance)
- Create: `CONTRIBUTING.md`
- Create: `CHANGELOG.md`
- Create: `CITATION.cff`
- Create: `LICENSE` (MIT, same as siblings)

- [ ] Commit each: `Add <doc>.md`

### Task 16.3: Final verification + v0.1.0 tag

- [ ] Run the full test suite: `pytest -v`
- [ ] Run validation gate: `python validation/validate_bvid_public.py --gate`
- [ ] Run `ruff check src tests && black --check src tests`
- [ ] Manually walk through the GUI: impact-driven path + damage-driven path + sweep + config save/load + PNG export.

- [ ] Commit: `Prepare v0.1.0 release`
- [ ] Tag: `git tag v0.1.0 -m "BVID-FE v0.1.0: initial release"`

---

## Notes for the Executing Agent

- **Phase 6 is the first functional milestone.** After phase 6 the package has a working empirical pipeline and CLI. Stop and verify the pipeline end-to-end before starting phase 7 — it's the cheapest place to catch integration issues.
- **Phases 7–10 are the most expensive.** The FE tier is only useful once all of core, elements, solver, failure are in place. Do not cut corners on element / solver tests — they prevent debugging hell in phase 10.
- **Phase 14 (validation) is where calibration happens.** The `MATERIAL_LIBRARY` constants in `core/material.py` are placeholders until this phase. Do not over-tune before there is data to tune against.
- **GUI tests** should all run headless using `pytest-qt` with `QT_QPA_PLATFORM=offscreen` on CI.
- **After every phase**, run `pytest -v` and confirm all phases up to that point still pass. Keep the test suite green; fix fast when something breaks.
