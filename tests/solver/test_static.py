"""Tests for the linear static FE solver."""

import re

import numpy as np
import pytest

from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.elements.hex8 import Hex8Element
from bvidfe.solver import static
from bvidfe.solver.boundary import BoundaryCondition, uniaxial_x_bcs
from bvidfe.solver.static import (
    dirichlet_dofs_and_values,
    solve_linear_static,
)


def _single_cube_model():
    """Create a simple single-element cube for testing."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    node_coords = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    elem = Hex8Element(node_coords, m)
    # global DOF map: node i -> DOFs 3*i, 3*i+1, 3*i+2
    dof_map = np.arange(24)
    return node_coords, [elem], [dof_map], m


def test_solve_compression_returns_prescribed_displacement_at_xmax():
    """Test that prescribed displacement is achieved at x_max."""
    node_coords, elems, dofs, m = _single_cube_model()
    n_dof = 3 * node_coords.shape[0]
    bcs = uniaxial_x_bcs(node_coords, applied_strain=-0.01)
    u = solve_linear_static(elems, dofs, n_dof, bcs)
    assert u.shape == (n_dof,)
    # Node 1 is at x=1 (x_max); u_x should be ≈ -0.01
    assert abs(u[3 * 1 + 0] - (-0.01)) < 1e-3


def test_solve_preserves_rigid_y_displacement_at_xmin():
    """Test that clamped boundary conditions are preserved at x_min."""
    node_coords, elems, dofs, m = _single_cube_model()
    n_dof = 3 * node_coords.shape[0]
    bcs = uniaxial_x_bcs(node_coords, applied_strain=-0.01)
    u = solve_linear_static(elems, dofs, n_dof, bcs)
    # Node 0 at x=y=z=0: u_x=u_y=u_z=0
    for k in range(3):
        assert abs(u[0 + k]) < 1e-6


def test_solve_axial_strain_matches_prescribed_for_single_cube():
    """Axial strain in a single cube under uniform compression should match the prescribed value."""
    node_coords, elems, dofs, m = _single_cube_model()
    n_dof = 3 * node_coords.shape[0]
    applied_strain = -0.01
    bcs = uniaxial_x_bcs(node_coords, applied_strain=applied_strain)
    u = solve_linear_static(elems, dofs, n_dof, bcs)
    # Node 2 at x=1: u_x ≈ applied_strain * Lx
    assert abs(u[3 * 2 + 0] - applied_strain) < 1e-3


# --- Dirichlet elimination + ILU-GMRES ----------------------------------------


def _small_panel_model(delaminated: bool = True):
    """A 20 x 10 mm, 4-ply fe3d mesh (990 DOFs) with the production BCs."""
    from bvidfe.analysis.config import AnalysisConfig, MeshParams
    from bvidfe.analysis.fe_mesh import build_fe_mesh
    from bvidfe.analysis.fe_tier import _build_elements
    from bvidfe.core.geometry import PanelGeometry
    from bvidfe.core.laminate import Laminate
    from bvidfe.damage.state import DamageState, DelaminationEllipse

    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0),
        tier="fe3d",
        damage=DamageState(delaminations=[], dent_depth_mm=0.0),
        mesh=MeshParams(in_plane_size_mm=2.0),
    )
    delams = [DelaminationEllipse(1, (10, 5), 6, 3, 0)] if delaminated else []
    damage = DamageState(delaminations=delams, dent_depth_mm=0.2 if delaminated else 0.0)
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], cfg.layup_deg, cfg.ply_thickness_mm)
    mesh = build_fe_mesh(cfg, damage)
    elements = _build_elements(mesh, lam)
    bcs = uniaxial_x_bcs(mesh.node_coords, -0.01, boundary="simply_supported")
    return mesh, elements, bcs


def _dense_reference(mesh, elements, bcs):
    """Dense elimination solve with numpy, as an independent reference."""
    from bvidfe.solver.assembler import assemble_global_stiffness

    K = assemble_global_stiffness(elements, mesh.element_dof_maps, mesh.n_dof).toarray()
    prescribed = {bc.dof: bc.value for bc in bcs}
    dofs = np.array(sorted(prescribed))
    values = np.array([prescribed[d] for d in dofs])
    free = np.setdiff1d(np.arange(mesh.n_dof), dofs)
    u = np.zeros(mesh.n_dof)
    u[dofs] = values
    u[free] = np.linalg.solve(K[np.ix_(free, free)], -K[np.ix_(free, dofs)] @ values)
    return u


@pytest.mark.parametrize("method", ["auto", "iterative", "direct"])
def test_elimination_solve_matches_dense_reference(method):
    mesh, elements, bcs = _small_panel_model()
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs, method=method)
    u_ref = _dense_reference(mesh, elements, bcs)
    assert np.linalg.norm(u - u_ref) <= 1e-8 * np.linalg.norm(u_ref)
    # Constrained DOFs carry their prescribed values exactly.
    for bc in bcs:
        assert u[bc.dof] == bc.value


def test_elimination_agrees_with_legacy_penalty_solve():
    """The 1e10 penalty only enforces the BCs to ~K/penalty; the two
    solves agree to that level."""
    mesh, elements, bcs = _small_panel_model()
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs)
    u_pen = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs, method="penalty")
    assert np.linalg.norm(u - u_pen) <= 1e-6 * np.linalg.norm(u)


def test_external_force_with_elimination():
    """A load on a free DOF enters the reduced right-hand side."""
    mesh, elements, _ = _small_panel_model(delaminated=False)
    x = mesh.node_coords[:, 0]
    bcs = [
        BoundaryCondition(dof=3 * int(n) + k, value=0.0)
        for n in np.flatnonzero(x == 0)
        for k in range(3)
    ]
    F = np.zeros(mesh.n_dof)
    tip = np.flatnonzero(x == x.max())
    F[3 * tip + 2] = 1.0  # transverse tip load, cantilever
    u_it = solve_linear_static(
        elements, mesh.element_dof_maps, mesh.n_dof, bcs, F, method="iterative"
    )
    u_dir = solve_linear_static(
        elements, mesh.element_dof_maps, mesh.n_dof, bcs, F, method="direct"
    )
    assert u_it[3 * tip + 2].min() > 0
    assert np.linalg.norm(u_it - u_dir) <= 1e-8 * np.linalg.norm(u_dir)


def test_duplicate_bcs_are_merged():
    """A DOF listed twice with the same value (as the builder once listed the
    loaded-edge bottom nodes' u_z = 0) collapses to one constraint and does
    not shift the solution."""
    mesh, elements, unique = _small_panel_model()
    bcs = unique + unique[::7]
    dofs, values = dirichlet_dofs_and_values(bcs, mesh.n_dof)
    assert len(bcs) > len(dofs) == len({bc.dof for bc in bcs})
    assert np.all(np.diff(dofs) > 0)
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs)
    u_unique = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, unique)
    np.testing.assert_array_equal(u, u_unique)


def test_conflicting_bcs_raise():
    bcs = [BoundaryCondition(dof=3, value=0.0), BoundaryCondition(dof=3, value=0.1)]
    with pytest.raises(ValueError, match="conflicting boundary conditions on DOF 3"):
        dirichlet_dofs_and_values(bcs, 24)
    # Agreeing duplicates (to round-off) are fine.
    dofs, values = dirichlet_dofs_and_values(
        [BoundaryCondition(3, 0.1), BoundaryCondition(3, 0.1 + 1e-17)], 24
    )
    assert dofs.tolist() == [3] and values.tolist() == [0.1]
    with pytest.raises(ValueError, match="out of range"):
        dirichlet_dofs_and_values([BoundaryCondition(24, 0.0)], 24)


def test_gmres_non_convergence_falls_back_to_direct(caplog):
    mesh, elements, bcs = _small_panel_model()
    u_ref = _dense_reference(mesh, elements, bcs)
    with caplog.at_level("INFO", logger="bvidfe.fe3d.solver"):
        u = solve_linear_static(
            elements, mesh.element_dof_maps, mesh.n_dof, bcs, method="iterative", maxiter=1
        )
    assert "ILU-GMRES did not reach rtol=1e-10 in 1 iterations" in caplog.text
    assert "using the direct solver" in caplog.text
    assert np.linalg.norm(u - u_ref) <= 1e-8 * np.linalg.norm(u_ref)


def test_iterative_solve_is_ilu_gmres_not_cg(monkeypatch, caplog):
    """The ILU preconditioner is not symmetric, so CG is not used at all
    (it ran 200 iterations without converging before GMRES on every 5 mm
    panel): the iterative solve is ILU-GMRES, which converges on the first
    attempt, without a fallback."""

    def _no_cg(*args, **kwargs):
        raise AssertionError("CG must not be called")

    monkeypatch.setattr(static.spla, "cg", _no_cg)
    mesh, elements, bcs = _small_panel_model()
    with caplog.at_level("INFO", logger="bvidfe.fe3d.solver"):
        u = solve_linear_static(
            elements, mesh.element_dof_maps, mesh.n_dof, bcs, method="iterative"
        )
    assert "did not reach" not in caplog.text
    u_ref = _dense_reference(mesh, elements, bcs)
    assert np.linalg.norm(u - u_ref) <= 1e-8 * np.linalg.norm(u_ref)


def test_gmres_breakdown_falls_back_without_a_numpy_warning(monkeypatch):
    """A Krylov solver breaks down when the preconditioner maps the
    residual to zero (in CG, r . M r = 0, seen with the real ILU on macOS /
    Python 3.12 for this cantilever): GMRES then normalises a zero vector,
    dividing by zero, which numpy reports as a RuntimeWarning ("divide by
    zero encountered in scalar divide"). The solver already falls back on
    the non-finite result, so the warning is noise and must not escape
    (pytest turns it into an error). The preconditioner below projects out
    the right-hand side, so M b = 0 for the single-entry tip load."""
    mesh, elements, _ = _small_panel_model(delaminated=False)
    x = mesh.node_coords[:, 0]
    bcs = [
        BoundaryCondition(dof=3 * int(n) + k, value=0.0)
        for n in np.flatnonzero(x == 0)
        for k in range(3)
    ]
    F = np.zeros(mesh.n_dof)
    F[3 * int(np.flatnonzero(x == x.max())[-1]) + 2] = 1.0
    fixed, _ = dirichlet_dofs_and_values(bcs, mesh.n_dof)
    b = np.delete(F, fixed)  # the reduced right-hand side (u_c = 0)

    class _AnnihilatesLoad:
        def solve(self, v):
            return v - b * (b @ v) / (b @ b)

    monkeypatch.setattr(static, "_ilu", lambda _K: _AnnihilatesLoad())
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs, F, method="iterative")
    u_dir = solve_linear_static(
        elements, mesh.element_dof_maps, mesh.n_dof, bcs, F, method="direct"
    )
    assert np.linalg.norm(u - u_dir) <= 1e-8 * np.linalg.norm(u_dir)


def test_ilu_failure_falls_back_to_direct(monkeypatch):
    mesh, elements, bcs = _small_panel_model()

    def _broken_ilu(_K):
        raise RuntimeError("Factor is exactly singular")

    monkeypatch.setattr(static, "_ilu", _broken_ilu)
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs, method="iterative")
    u_ref = _dense_reference(mesh, elements, bcs)
    assert np.linalg.norm(u - u_ref) <= 1e-8 * np.linalg.norm(u_ref)


def test_unknown_method_raises():
    node_coords, elems, dofs, m = _single_cube_model()
    bcs = uniaxial_x_bcs(node_coords, applied_strain=-0.01)
    with pytest.raises(ValueError, match="unknown method"):
        solve_linear_static(elems, dofs, 24, bcs, method="cholesky")


# ------------------------------------------------------- stage log (#110)


def _stage_messages(caplog) -> list[str]:
    return [r.getMessage() for r in caplog.records if r.name == "bvidfe.fe3d.solver"]


@pytest.mark.parametrize(
    "method, how",
    [
        ("direct", r"solve done: direct LU \(\d+\.\d\ds\)"),
        ("iterative", r"solve done: ILU-GMRES, \d+ iterations \(\d+\.\d\ds\)"),
        ("penalty", r"solve done: spsolve \(\d+\.\d\ds\)"),
    ],
)
def test_solve_logs_each_stage_with_dofs_method_and_time(caplog, method, how):
    """One INFO line at the start and end of assembly, BCs and the solve,
    with DOF counts, the solver used, its iterations and the wall time."""

    mesh, elements, bcs = _small_panel_model()
    with caplog.at_level("INFO", logger="bvidfe.fe3d.solver"):
        solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs, method=method)
    msgs = _stage_messages(caplog)
    stages = [m.split(":")[0] for m in msgs]
    assert stages == [
        "assembly start",
        "assembly done",
        "BCs start",
        "BCs done",
        "solve start",
        "solve done",
    ], msgs
    assert msgs[0] == f"assembly start: {len(elements)} elements, {mesh.n_dof} DOFs"
    assert re.fullmatch(r"assembly done: \d+ nonzeros \(\d+\.\d\ds\)", msgs[1])
    if method != "penalty":
        n_fixed = len({bc.dof for bc in bcs})
        assert msgs[3].startswith(
            f"BCs done: {n_fixed} constrained DOFs eliminated, {mesh.n_dof - n_fixed} free"
        )
        assert msgs[4].startswith(f"solve start: method={method}, {mesh.n_dof - n_fixed} free")
    assert re.fullmatch(how, msgs[5]), msgs[5]


def test_solve_log_names_the_fallback_that_produced_the_answer(caplog):
    mesh, elements, bcs = _small_panel_model()
    with caplog.at_level("INFO", logger="bvidfe.fe3d.solver"):
        solve_linear_static(
            elements, mesh.element_dof_maps, mesh.n_dof, bcs, method="iterative", maxiter=1
        )
    (done,) = [m for m in _stage_messages(caplog) if m.startswith("solve done")]
    after = r"\(after 1 ILU-GMRES iterations\)"
    assert re.fullmatch(rf"solve done: direct LU {after} \(\d+\.\d\ds\)", done), done
