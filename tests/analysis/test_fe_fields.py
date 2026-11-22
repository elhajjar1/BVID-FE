"""fe3d fields for the CAI load-ramp animation (``bvidfe.analysis.fe_fields``).

``fe3d_field_history`` returns the reference displacement, the per-element
failure multipliers and the first-ply-failure and residual strengths of an
fe3d config. They must equal what ``BvidAnalysis.run()`` reports. ``run()``
must not call the module, and the helpers extracted from ``fe_tier`` for it
must leave ``run()`` bit-identical.

The main-suite cases are 150 x 100 mm panels at 10 mm and 40 x 30 mm panels
at 2.5 mm, where the 2.5 mm averaging disk spans neighbouring elements.
The 150 x 100 mm panel at 2.5 and 5 mm is ``slow``.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from dataclasses import replace
from functools import lru_cache
from pathlib import Path
from unittest import mock

import numpy as np
import pytest
from scipy.ndimage import convolve, minimum_filter

from bvidfe import AnalysisConfig, BvidAnalysis, DamageState, DelaminationEllipse, MeshParams
from bvidfe.analysis import bvid, fe_fields, fe_tier
from bvidfe.analysis.fe_fields import Fe3dFields, Fe3dPanelFields, fe3d_field_history
from bvidfe.analysis.fe_mesh import element_grid_shape
from bvidfe.analysis.fe_tier import (
    FPF_AVERAGING_RADIUS_MM,
    _build_elements,
    _disk_kernel,
    _failure_stress_material,
    _index_coefficients,
    _resolve_material,
    _smallest_positive_root,
)
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.failure.larc05 import larc05_index_batch
from bvidfe.failure.puck import puck_index_batch
from bvidfe.failure.tsai_wu import tsai_wu_index_batch
from bvidfe.solver.assembler import assemble_global_stiffness
from tests.analysis._fe3d_fixtures import QI8, impact_cfg

REPO = Path(__file__).resolve().parents[2]
UNDAMAGED = DamageState(delaminations=[], dent_depth_mm=0.0)
# Two delaminations under the centre of a 40 x 30 mm panel.
SMALL_DAMAGE = DamageState(
    delaminations=[
        DelaminationEllipse(2, (20.0, 15.0), 9.0, 6.0, 0.0),
        DelaminationEllipse(5, (20.0, 15.0), 8.0, 5.0, 90.0),
    ],
    dent_depth_mm=0.0,
)


def _small_cfg(loading: str, beta: float = 0.0) -> AnalysisConfig:
    """8-ply QI IM7/8552, 40 x 30 mm, 2.5 mm mesh, ``SMALL_DAMAGE``."""
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(40.0, 30.0),
        loading=loading,
        tier="fe3d",
        damage=SMALL_DAMAGE,
        transverse_stress_ratio=beta,
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.5),
    )


CONFIGS: dict[str, AnalysisConfig] = {
    "cai": impact_cfg("compression", 10.0),
    "tai": impact_cfg("tension", 20.0),
    "cai_biaxial": replace(impact_cfg("compression", 10.0), transverse_stress_ratio=0.5),
    "tai_biaxial": replace(impact_cfg("tension", 20.0), transverse_stress_ratio=-0.5),
    # Delamination growth governs (see test_fe3d_knockdown_normalization).
    "cai_growth": replace(impact_cfg("compression", 5.0), layup_deg=[45, -45, 90, 90, -45, 45]),
    "cai_2p5mm": _small_cfg("compression"),
    "tai_2p5mm": _small_cfg("tension"),
    "cai_2p5mm_biaxial": _small_cfg("compression", beta=-0.5),
}
GOVERNING = {name: "first_ply_failure" for name in CONFIGS} | {"cai_growth": "delamination_growth"}

# run() of the base code (main at caec56b, before fe_fields and the helpers
# extracted for it), captured on linux x86_64 (numpy 2.4.6, scipy 1.17.1,
# OMP_NUM_THREADS=1): (residual_strength_MPa, knockdown, E11_knockdown,
# critical_sublaminate, solves). There the current code reproduces every
# float bit for bit (float.hex), and every note. The 1e-9 tolerance leaves
# room for other platforms' round-off: solving these systems by direct LU
# instead of ILU-GMRES moves them by at most 1e-11.
PINNED_RUN = {
    "cai": (225.81067241022077, 0.3755457816861775, 0.8990954628082054, 5, 2),
    "tai": (509.70724222756286, 0.6766850972446213, 0.9855173443918511, None, 2),
    "cai_biaxial": (246.13533775698113, 0.3497476837491107, 0.8990954628082054, 5, 4),
    "tai_biaxial": (431.8784258661742, 0.6692645496230041, 0.9855173443918511, 4, 4),
    "cai_growth": (266.1697359926951, 0.5188528426604728, 0.9584092906188425, 4, 2),
    "cai_2p5mm": (198.4584599290917, 0.33005631076156455, 0.7524524945761623, 5, 2),
    "tai_2p5mm": (670.1792531123286, 0.8897270344476976, 0.9831826517027581, None, 2),
    "cai_2p5mm_biaxial": (193.45651234621198, 0.3779751912800412, 0.7524524945761623, 5, 4),
}


def _laminate(cfg: AnalysisConfig) -> Laminate:
    return Laminate(_resolve_material(cfg), cfg.layup_deg, cfg.ply_thickness_mm)


def _counting_solves(counter: list):
    solve = fe_tier.solve_linear_static

    def counted(*args, **kwargs):
        counter.append(1)
        return solve(*args, **kwargs)

    return mock.patch.object(fe_tier, "solve_linear_static", counted)


@lru_cache(maxsize=None)
def _computed(name: str):
    """``(run() result, its solve count, fields, fields' solve count, the
    fe3d calls run() made)``. The calls are ``{(damaged, beta):
    (first-ply failure stress, E_x)}`` from ``fe3d_cai_channels_and_modulus``
    or ``fe3d_tai_and_modulus``."""
    cfg = CONFIGS[name]
    calls: dict = {}
    real_cai, real_tai = bvid.fe3d_cai_channels_and_modulus, bvid.fe3d_tai_and_modulus

    def cai(cfg_, damage, lam, cap):
        out = real_cai(cfg_, damage, lam, cap)
        calls[(bool(damage.delaminations), cfg_.transverse_stress_ratio)] = out[1:]
        return out

    def tai(cfg_, damage, lam, cap):
        out = real_tai(cfg_, damage, lam, cap)
        calls[(bool(damage.delaminations), cfg_.transverse_stress_ratio)] = out
        return out

    run_solves: list = []
    with (
        mock.patch.object(bvid, "fe3d_cai_channels_and_modulus", cai),
        mock.patch.object(bvid, "fe3d_tai_and_modulus", tai),
        _counting_solves(run_solves),
    ):
        result = BvidAnalysis(cfg).run()
    field_solves: list = []
    with _counting_solves(field_solves):
        fields = fe3d_field_history(cfg)
    return result, len(run_solves), fields, len(field_solves), calls


def _panels(fields: Fe3dFields) -> tuple[Fe3dPanelFields, Fe3dPanelFields]:
    return fields.damaged, fields.undamaged


# ---------------------------------------------------------------------------
# Equal to run()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", CONFIGS)
def test_fields_match_run(name):
    """First-ply failure (damaged and undamaged, fe3d's own stresses), the
    growth channel, the governing one and every strength and knockdown
    are run()'s, to 1e-9 (on the biaxial configs, the load ray's)."""
    cfg = CONFIGS[name]
    r, _, f, _, calls = _computed(name)
    beta = cfg.transverse_stress_ratio
    for panel, damaged in ((f.damaged, True), (f.undamaged, False)):
        fpf, modulus = calls[(damaged, beta)]
        assert panel.fpf_stress_MPa == pytest.approx(fpf, rel=1e-9)
        assert panel.secant_modulus_MPa == pytest.approx(modulus, rel=1e-9)
    assert f.pristine_strength_MPa == pytest.approx(r.pristine_strength_MPa, rel=1e-9)
    assert f.residual_strength_MPa == pytest.approx(r.residual_strength_MPa, rel=1e-9)
    assert f.knockdown == pytest.approx(r.knockdown, rel=1e-9)
    assert f.E11_knockdown == pytest.approx(r.E11_knockdown, rel=1e-9)
    assert f.critical_sublaminate == r.critical_sublaminate
    assert f.governing == GOVERNING[name]
    assert f.criterion == ("tsai_wu" if cfg.loading == "tension" else "larc05")
    tier = "TAI" if cfg.loading == "tension" else "CAI"
    if math.isfinite(f.growth_stress_MPa):
        assert f.critical_sublaminate == f.growth.interface_index
        assert any(
            n.startswith(f"fe3d delamination growth: interface {f.critical_sublaminate} grows at ")
            and f"{f.growth_stress_MPa:.1f} MPa" in n
            for n in r.notes
        ), r.notes
        channel = f.governing.replace("_", " ").replace("first ply", "first-ply")
        assert any(n.startswith(f"fe3d {tier}: {channel} governs") for n in r.notes), r.notes
    else:
        assert f.governing == "first_ply_failure"
        assert f.growth is None and tier == "TAI" and beta == 0
    assert f.fpf_residual_MPa == pytest.approx(f.pristine_strength_MPa * f.fpf_knockdown)
    assert f.residual_strength_MPa == min(f.growth_stress_MPa, f.fpf_residual_MPa)
    if beta != 0:
        # The biaxial fields are the load ray's, not the uniaxial residual.
        assert f.damaged.sigma_y_ref_MPa != 0.0
        assert calls[(True, 0.0)][0] != pytest.approx(f.damaged.fpf_stress_MPa, rel=1e-3)


@pytest.mark.parametrize("name", CONFIGS)
def test_run_matches_the_pinned_base_code(name):
    """run() is unchanged by fe_fields and the helpers extracted for it."""
    r, run_solves, _, _, _ = _computed(name)
    residual, knockdown, e11, critical, solves = PINNED_RUN[name]
    assert r.residual_strength_MPa == pytest.approx(residual, rel=1e-9)
    assert r.knockdown == pytest.approx(knockdown, rel=1e-9)
    assert r.E11_knockdown == pytest.approx(e11, rel=1e-9)
    assert r.critical_sublaminate == critical
    assert r.field_results is None
    assert run_solves == solves


def test_run_does_not_call_fe_fields(monkeypatch):
    """With fe3d_field_history and its per-mesh worker replaced by functions
    that raise, run() still runs, and solves the same number of systems as
    before (2, or 4 biaxial); fe3d_field_history solves as many."""
    names = ("cai", "tai_biaxial")
    field_solves = {name: _computed(name)[3] for name in names}

    def boom(*args, **kwargs):
        raise AssertionError("run() called fe_fields")

    for attr in ("fe3d_field_history", "_panel_fields"):
        monkeypatch.setattr(fe_fields, attr, boom)
    for name in names:
        solves: list = []
        with _counting_solves(solves):
            r = BvidAnalysis(CONFIGS[name]).run()
        assert len(solves) == PINNED_RUN[name][4] == field_solves[name]
        assert r.residual_strength_MPa == pytest.approx(PINNED_RUN[name][0], rel=1e-9)


def test_run_does_not_import_fe_fields():
    """In a fresh interpreter, an fe3d run() leaves fe_fields unimported."""
    code = (
        "import sys\n"
        "from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, MeshParams, "
        "PanelGeometry\n"
        "cfg = AnalysisConfig(material='IM7/8552', layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],"
        " ply_thickness_mm=0.152, panel=PanelGeometry(150.0, 100.0), tier='fe3d',"
        " impact=ImpactEvent(10.0), mesh=MeshParams(in_plane_size_mm=10.0))\n"
        "BvidAnalysis(cfg).run()\n"
        "assert 'bvidfe.analysis.fe_tier' in sys.modules\n"
        "print('bvidfe.analysis.fe_fields' in sys.modules)\n"
    )
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO / "src"), env.get("PYTHONPATH", "")])
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
    )
    assert out.stdout.strip() == "False", out.stderr


# ---------------------------------------------------------------------------
# Where first-ply failure is set
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", CONFIGS)
def test_fpf_element_is_where_first_ply_failure_is_set(name):
    """At fpf_stress_MPa the disk-averaged load ratio reaches 1 at
    fpf_element, the element of least c_averaged, and nowhere earlier. On
    the damaged panels that element is in the damage zone."""
    _, _, f, _, _ = _computed(name)
    for panel in _panels(f):
        c_min = panel.c_averaged.min()
        assert panel.c_averaged[panel.fpf_element] == c_min
        assert not panel.fpf_capped and panel.fpf_multiplier == c_min
        assert panel.fpf_stress_MPa == pytest.approx(c_min * panel.sigma_ref_MPa, rel=1e-12)
        ratio = panel.load_ratio(panel.fpf_stress_MPa)
        assert ratio[panel.fpf_element] == pytest.approx(1.0, rel=1e-12)
        assert ratio.max() == pytest.approx(1.0, rel=1e-12)
        assert panel.load_ratio(0.999 * panel.fpf_stress_MPa).max() < 1.0
        # The pointwise check fails first (see test_pointwise_and_averaged...).
        assert panel.load_ratio(panel.fpf_stress_MPa, pointwise=True).max() >= 1.0
    assert f.damaged.mesh.in_plane_damage_factors[f.damaged.fpf_element] < 1.0


@pytest.mark.parametrize("name", CONFIGS)
def test_pointwise_and_averaged_multipliers(name):
    """Each element's disk-averaged multiplier is at least the least
    pointwise one in its disk, or 1 (with the TAI notch, both times the
    notch): below that every Gauss point of the disk is below failure. At
    10 mm the disk is the element alone, so the relation holds element by
    element. Globally, pointwise failure comes first."""
    cfg = CONFIGS[name]
    _, _, f, _, _ = _computed(name)
    for panel in _panels(f):
        K = np.ones(panel.mesh.n_elements) if panel.notch is None else panel.notch
        nz, ny, nx = panel.grid_shape
        kernel = _disk_kernel(cfg.panel.Lx_mm / nx, cfg.panel.Ly_mm / ny, FPF_AVERAGING_RADIUS_MM)
        raw = (panel.c_pointwise * K).reshape(nz, ny, nx)
        disk_min = np.stack(
            [
                minimum_filter(layer, footprint=kernel > 0, mode="constant", cval=np.inf)
                for layer in raw
            ]
        ).ravel()
        lower = np.minimum(1.0, disk_min)
        assert np.all(panel.c_averaged * K >= lower * (1 - 1e-12))
        if np.count_nonzero(kernel) == 1:
            assert np.all(panel.c_averaged >= np.minimum(1 / K, panel.c_pointwise) * (1 - 1e-12))
        if panel.notch is None:
            assert panel.c_averaged.min() >= panel.c_pointwise.min()
        assert np.all(panel.c_pointwise > 0) and np.all(panel.c_averaged > 0)
    # On the undamaged panel pointwise failure is at most 18% earlier (at
    # 2.5 mm, from the free-edge and corner stresses the disk averages out).
    gap = f.undamaged.c_averaged.min() / f.undamaged.c_pointwise.min()
    assert 1.0 <= gap < 1.25


def test_disk_spans_neighbours_on_the_2p5mm_mesh():
    """The 2.5 mm cases exercise the disk over neighbouring elements: there
    c_averaged >= c_pointwise fails for some elements, as it may."""
    _, _, f, _, _ = _computed("cai_2p5mm")
    kernel = _disk_kernel(2.5, 2.5, FPF_AVERAGING_RADIUS_MM)
    assert np.count_nonzero(kernel) > 1
    assert np.any(f.damaged.c_averaged < f.damaged.c_pointwise)


# ---------------------------------------------------------------------------
# Undamaged panel
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_undamaged_panel(loading):
    """With no damage both panels are the same: knockdowns exactly 1, no
    growth, and c_averaged uniform in each layer away from the edges (the
    y faces are free edges; the disk is clipped at every edge)."""
    cfg = impact_cfg(loading, 10.0)
    f = fe3d_field_history(cfg, damage=UNDAMAGED)
    assert f.damage is UNDAMAGED and f.growth is None
    assert f.knockdown == 1.0 and f.fpf_knockdown == 1.0 and f.E11_knockdown == 1.0
    assert f.residual_strength_MPa == f.pristine_strength_MPa
    assert f.governing == "first_ply_failure"
    assert f.damaged.fpf_stress_MPa == f.undamaged.fpf_stress_MPa
    np.testing.assert_array_equal(f.damaged.c_averaged, f.undamaged.c_averaged)
    assert f.residual_applied_stress_MPa == f.undamaged.fpf_stress_MPa
    for name in ("cai", "tai", "cai_2p5mm", "tai_2p5mm"):
        if CONFIGS[name].loading != loading:
            continue
        panels = [_computed(name)[2].undamaged]
        if name in ("cai", "tai"):
            panels.append(f.undamaged)
        for panel in panels:
            nz, ny, nx = panel.grid_shape
            interior = panel.c_averaged.reshape(nz, ny, nx)[:, 2 : ny - 2, 2 : nx - 2]
            lo, hi = interior.min(axis=(1, 2)), interior.max(axis=(1, 2))
            assert np.all((hi - lo) / lo < 5e-3), (hi - lo) / lo
            assert np.all(panel.mesh.in_plane_damage_factors == 1.0)
            assert panel.notch is None


def test_undamaged_reference_is_the_undamaged_config():
    """A damaged config's undamaged panel is that of the same config without
    damage."""
    damaged = _computed("cai")[2]
    clean = fe3d_field_history(CONFIGS["cai"], damage=UNDAMAGED)
    assert damaged.undamaged.fpf_stress_MPa == pytest.approx(clean.damaged.fpf_stress_MPa)
    np.testing.assert_allclose(damaged.undamaged.u_ref, clean.damaged.u_ref, rtol=0, atol=1e-12)


# ---------------------------------------------------------------------------
# Reference displacement, reactions and the load-displacement line
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["cai", "tai", "cai_biaxial", "tai_biaxial"])
def test_u_ref_reaction_and_scaling(name):
    """The reaction K u_ref on the x_max face over the gross section is
    sigma_ref_MPa (signed as the reference strain), and on the y_max face
    over Lx h it is sigma_y_ref_MPa = beta sigma_x; the platen moves by the
    reference end shortening. At any stress the scaled displacement carries
    that stress, and the load-displacement line goes through it."""
    cfg = CONFIGS[name]
    beta = cfg.transverse_stress_ratio
    _, _, f, _, _ = _computed(name)
    lam = _laminate(cfg)
    h = lam.thickness_mm
    for panel in _panels(f):
        mesh = panel.mesh
        elements = _build_elements(mesh, lam, tension=cfg.loading == "tension")
        K = assemble_global_stiffness(elements, mesh.element_dof_maps, mesh.n_dof)
        x, y = mesh.node_coords[:, 0], mesh.node_coords[:, 1]
        x_max, x_min = np.flatnonzero(np.isclose(x, x.max())), np.flatnonzero(np.isclose(x, 0))
        y_max = np.flatnonzero(np.isclose(y, y.max()))
        R = K @ panel.u_ref
        sigma_x = R[3 * x_max].sum() / (cfg.panel.Ly_mm * h)
        sigma_y = R[3 * y_max + 1].sum() / (cfg.panel.Lx_mm * h)
        assert panel.gross_section_mm2 == pytest.approx(cfg.panel.Ly_mm * h)
        assert abs(sigma_x) == pytest.approx(panel.sigma_ref_MPa, rel=1e-9)
        assert np.sign(sigma_x) == np.sign(panel.reference_strain)
        assert panel.reaction_force_ref_N == pytest.approx(abs(R[3 * x_max].sum()), rel=1e-9)
        tol = 1e-9 * panel.sigma_ref_MPa
        assert sigma_y == pytest.approx(panel.sigma_y_ref_MPa, rel=1e-9, abs=tol)
        assert sigma_y == pytest.approx(beta * sigma_x, rel=1e-9, abs=tol)
        u = panel.u_ref.reshape(-1, 3)
        assert panel.reference_strain == (0.05 if cfg.loading == "tension" else -0.05)
        np.testing.assert_allclose(u[x_max, 0], panel.reference_strain * cfg.panel.Lx_mm)
        np.testing.assert_array_equal(u[x_min, 0], 0.0)
        assert panel.end_shortening_ref_mm == pytest.approx(0.05 * cfg.panel.Lx_mm)

        stress = 0.37 * panel.fpf_stress_MPa
        u_s = panel.displacement(stress)
        assert u_s.shape == (mesh.n_nodes, 3)
        R_s = K @ u_s.ravel()
        assert abs(R_s[3 * x_max].sum()) / panel.gross_section_mm2 == pytest.approx(stress)
        assert panel.reaction_force_N(stress) == pytest.approx(abs(R_s[3 * x_max].sum()))
        assert panel.end_shortening_mm(stress) == pytest.approx(abs(u_s[x_max, 0]).max())
        slope = panel.reaction_force_ref_N / panel.end_shortening_ref_mm
        assert panel.reaction_force_N(stress) == pytest.approx(
            slope * panel.end_shortening_mm(stress)
        )
        assert not np.any(panel.displacement(0.0))
        if beta == 0:
            # The uniaxial line's slope is the secant modulus.
            assert slope * cfg.panel.Lx_mm / panel.gross_section_mm2 == pytest.approx(
                panel.secant_modulus_MPa
            )


def test_load_ramp_end_and_stress_scales():
    """The ramp ends at the damaged first-ply failure where it governs; where
    growth governs, at the growth stress on the applied scale, when no
    element has failed yet. reported_stress and applied_stress invert each
    other, and map undamaged first-ply failure to the pristine strength."""
    fpf = _computed("cai")[2]
    assert fpf.residual_applied_stress_MPa == fpf.damaged.fpf_stress_MPa
    assert fpf.reported_stress(fpf.damaged.fpf_stress_MPa) == pytest.approx(
        fpf.residual_strength_MPa, rel=1e-12
    )
    growth = _computed("cai_growth")[2]
    end = growth.residual_applied_stress_MPa
    assert end == pytest.approx(growth.applied_stress(growth.growth_stress_MPa), rel=1e-12)
    assert end < growth.damaged.fpf_stress_MPa
    assert growth.damaged.load_ratio(end).max() < 1.0
    for f in (fpf, growth):
        s = 123.4
        assert f.applied_stress(f.reported_stress(s)) == pytest.approx(s, rel=1e-12)
        assert f.reported_stress(f.undamaged.fpf_stress_MPa) == pytest.approx(
            f.pristine_strength_MPa, rel=1e-12
        )


def test_fields_layout():
    """Shapes, read-only arrays, and the boundary faces (4 nodes and the
    parent element) of the structured grid, shared by the two panels."""
    cfg = CONFIGS["tai"]
    _, _, f, _, _ = _computed("tai")
    nz, ny, nx = element_grid_shape(cfg)
    for panel in _panels(f):
        n = panel.mesh.n_elements
        assert panel.grid_shape == (nz, ny, nx) and n == nz * ny * nx
        assert panel.u_ref.shape == (panel.mesh.n_dof,)
        for a in (panel.u_ref, panel.c_pointwise, panel.c_averaged, panel.boundary_faces):
            assert not a.flags.writeable
        assert panel.c_pointwise.shape == panel.c_averaged.shape == (n,)
        faces = panel.boundary_faces
        assert faces.shape == (2 * (nx * ny + nx * nz + ny * nz), 5)
        assert set(np.unique(faces[:, 4])) <= set(range(n))
        assert np.all(faces[:, :4] < panel.mesh.n_nodes)
    assert f.undamaged.boundary_faces is f.damaged.boundary_faces
    # The TAI fibre-damage hole notches the elements under it, in every ply.
    assert f.undamaged.notch is None
    assert f.damaged.notch is not None
    under = (f.damaged.notch > 1.0).reshape(nz, ny * nx)
    assert np.any(under) and not np.all(under)
    assert np.all(under == under[0])


# ---------------------------------------------------------------------------
# The extracted fe_tier helpers are bit-identical to the code they came from
# ---------------------------------------------------------------------------


def _base_critical_strain_multiplier(cfg, mesh, elements, u_ref, criterion, notch=None):
    """fe_tier._critical_strain_multiplier of the base code (caec56b),
    verbatim but for the pointwise branch and the grid-size check."""
    a, b = _index_coefficients(cfg, mesh, elements, u_ref, criterion)
    nz, ny, nx = element_grid_shape(cfg)
    kernel = _disk_kernel(cfg.panel.Lx_mm / nx, cfg.panel.Ly_mm / ny, FPF_AVERAGING_RADIUS_MM)
    weight = convolve(np.ones((ny, nx)), kernel, mode="constant")
    a_avg = np.stack(
        [convolve(layer, kernel, mode="constant") / weight for layer in a.reshape(nz, ny, nx)]
    )
    b_avg = np.stack(
        [convolve(layer, kernel, mode="constant") / weight for layer in b.reshape(nz, ny, nx)]
    )
    if notch is None:
        return float(_smallest_positive_root(a_avg.ravel(), b_avg.ravel()).min())
    K = np.asarray(notch, dtype=float)
    return float(_smallest_positive_root(K * a_avg.ravel(), K * K * b_avg.ravel()).min())


def _base_critical_strain_multiplier_pointwise(cfg, mesh, elements, u_ref, criterion):
    """fe_tier._critical_strain_multiplier_pointwise of the base code,
    verbatim."""
    material = _resolve_material(cfg)
    in_plane = mesh.in_plane_damage_factors
    c_crit_min = np.inf
    for eidx, elem in enumerate(elements):
        dof_map = mesh.element_dof_maps[eidx]
        sigma_ref = _failure_stress_material(elem, u_ref[dof_map], in_plane[eidx])
        if criterion == "larc05":
            idx_ref = larc05_index_batch(material, sigma_ref)
            valid = idx_ref > 0
            if not bool(valid.any()):
                continue
            c_crit = np.where(valid, 1.0 / np.sqrt(np.where(valid, idx_ref, 1.0)), np.inf)
            c_crit_elem_min = float(c_crit.min())
        elif criterion == "puck":
            idx_ref = puck_index_batch(material, sigma_ref)
            valid = idx_ref > 0
            if not bool(valid.any()):
                continue
            c_crit = np.where(valid, 1.0 / np.where(valid, idx_ref, 1.0), np.inf)
            c_crit_elem_min = float(c_crit.min())
        else:
            idx_ref = tsai_wu_index_batch(material, sigma_ref)
            sigma_2 = 2.0 * sigma_ref
            idx_2 = tsai_wu_index_batch(material, sigma_2)
            b = (idx_2 - 2.0 * idx_ref) / 2.0
            a = idx_ref - b
            valid = idx_ref > 0
            tiny = 1e-14
            c_crit = np.full_like(idx_ref, np.inf, dtype=float)
            mask_lin = valid & (np.abs(b) < tiny) & (np.abs(a) >= tiny)
            c_lin = 1.0 / np.where(mask_lin, a, 1.0)
            c_crit = np.where(mask_lin, c_lin, c_crit)
            mask_quad = valid & (np.abs(b) >= tiny)
            disc = a * a + 4.0 * b
            mask_quad_real = mask_quad & (disc >= 0)
            sqrt_disc = np.sqrt(np.where(mask_quad_real, disc, 0.0))
            c_quad = (-a + sqrt_disc) / np.where(mask_quad_real, 2.0 * b, 1.0)
            c_crit = np.where(mask_quad_real, c_quad, c_crit)
            c_crit = np.where(c_crit > 0, c_crit, np.inf)
            c_crit_elem_min = float(c_crit.min())
        if c_crit_elem_min < c_crit_min:
            c_crit_min = c_crit_elem_min
    return c_crit_min


def _base_basis_failure_stress_and_modulus(cfg, basis, lam, criterion, notch, beta):
    """fe_tier._basis_failure_stress_and_modulus of the base code, verbatim
    but for the log lines."""
    u_ref, force_x, _ = fe_tier._ratio_displacement(cfg, basis, lam, beta)
    c_crit = _base_critical_strain_multiplier(
        cfg, basis.mesh, basis.elements, u_ref, criterion, notch
    )
    if not np.isfinite(c_crit) or c_crit <= 0:
        c_crit = 1.0
    c_crit = min(c_crit, 1.0)
    area = cfg.panel.Ly_mm * lam.thickness_mm
    force = abs(force_x)
    sigma = c_crit * force / area
    return sigma, abs(basis.force_x) / area / basis.strain_cap


@pytest.mark.parametrize(
    "name, criteria",
    [
        ("cai_2p5mm", ("larc05", "tsai_wu", "puck")),
        ("tai_2p5mm", ("tsai_wu",)),
        ("cai_2p5mm_biaxial", ("larc05",)),
    ],
)
def test_extracted_helpers_are_bit_identical_to_the_base_code(name, criteria):
    """On one solve (the same reference displacement), the refactored fe_tier
    functions give ``==`` the base code's results: the averaged multiplier
    (the minimum of _averaged_critical_multipliers) with and without the
    TAI notch, the pointwise one (of _pointwise_critical_multipliers) and
    the first-ply failure stress and modulus (_fpf_stress_and_modulus)."""
    cfg = CONFIGS[name]
    lam = _laminate(cfg)
    tension = cfg.loading == "tension"
    basis, _ = fe_tier._fe3d_basis(
        cfg, SMALL_DAMAGE, lam, label="test", tension=tension, strain_sign=1 if tension else -1
    )
    mesh, elements = basis.mesh, basis.elements
    notch = fe_tier._tai_fibre_damage_notch(cfg, SMALL_DAMAGE, lam, mesh) if tension else None
    u_ref, _, _ = fe_tier._ratio_displacement(cfg, basis, lam, cfg.transverse_stress_ratio)
    for criterion in criteria:
        if notch is not None:
            assert fe_tier._critical_strain_multiplier(
                cfg, mesh, elements, u_ref, criterion
            ) == _base_critical_strain_multiplier(cfg, mesh, elements, u_ref, criterion)
        assert fe_tier._critical_strain_multiplier(
            cfg, mesh, elements, u_ref, criterion, None
        ) == _base_critical_strain_multiplier_pointwise(cfg, mesh, elements, u_ref, criterion)
        args = (cfg, basis, lam, criterion, notch, cfg.transverse_stress_ratio)
        assert fe_tier._basis_failure_stress_and_modulus(
            *args
        ) == _base_basis_failure_stress_and_modulus(*args)


def test_extracted_growth_channels_are_the_base_calls():
    """_fe3d_cai_growth and _fe3d_tai_growth make the calls fe3d CAI and TAI
    made inline before (``==`` on the dataclass)."""
    from bvidfe.analysis.semi_analytical import weakest_sublaminate_growth

    cfg = CONFIGS["cai_2p5mm"]
    lam = _laminate(cfg)
    boundary = cfg.panel.boundary
    assert fe_tier._fe3d_cai_growth(cfg, SMALL_DAMAGE, lam) == weakest_sublaminate_growth(
        lam, SMALL_DAMAGE, boundary=boundary
    )
    biax = replace(cfg, transverse_stress_ratio=-0.5)
    assert fe_tier._fe3d_cai_growth(biax, SMALL_DAMAGE, lam) == weakest_sublaminate_growth(
        lam, SMALL_DAMAGE, boundary=boundary, transverse_stress_ratio=-0.5
    )
    tai = replace(cfg, loading="tension")
    assert fe_tier._fe3d_tai_growth(tai, SMALL_DAMAGE, lam) is None
    tai_biax = replace(tai, transverse_stress_ratio=-0.5)
    assert fe_tier._fe3d_tai_growth(tai_biax, SMALL_DAMAGE, lam) == weakest_sublaminate_growth(
        lam, SMALL_DAMAGE, boundary=boundary, loading="tension", transverse_stress_ratio=-0.5
    )


# ---------------------------------------------------------------------------
# The 150 x 100 mm reference panel at production meshes
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("loading, size_mm", [("compression", 2.5), ("tension", 5.0)])
def test_fields_match_run_on_the_reference_panel(loading, size_mm):
    """8-ply QI 150 x 100 mm at 2.5 mm (CAI, 10 J) and 5 mm (TAI, 20 J):
    the fields equal run(), first-ply failure is set at fpf_element, and the
    disk relation of test_pointwise_and_averaged_multipliers holds."""
    cfg = replace(
        impact_cfg(loading, 10.0 if loading == "compression" else 20.0),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=size_mm),
    )
    r = BvidAnalysis(cfg).run()
    f = fe3d_field_history(cfg)
    assert f.residual_strength_MPa == pytest.approx(r.residual_strength_MPa, rel=1e-9)
    assert f.knockdown == pytest.approx(r.knockdown, rel=1e-9)
    assert f.E11_knockdown == pytest.approx(r.E11_knockdown, rel=1e-9)
    assert f.critical_sublaminate == r.critical_sublaminate
    assert f.governing == "first_ply_failure"
    for panel in _panels(f):
        ratio = panel.load_ratio(panel.fpf_stress_MPa)
        assert ratio[panel.fpf_element] == pytest.approx(1.0, rel=1e-12)
        assert ratio.max() == pytest.approx(1.0, rel=1e-12)
        K = np.ones(panel.mesh.n_elements) if panel.notch is None else panel.notch
        nz, ny, nx = panel.grid_shape
        kernel = _disk_kernel(cfg.panel.Lx_mm / nx, cfg.panel.Ly_mm / ny, FPF_AVERAGING_RADIUS_MM)
        disk_min = np.stack(
            [
                minimum_filter(layer, footprint=kernel > 0, mode="constant", cval=np.inf)
                for layer in (panel.c_pointwise * K).reshape(nz, ny, nx)
            ]
        ).ravel()
        assert np.all(panel.c_averaged * K >= np.minimum(1.0, disk_min) * (1 - 1e-12))
