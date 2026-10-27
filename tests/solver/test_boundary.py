import numpy as np
import pytest
import scipy.sparse as sp

from bvidfe.damage.state import DelaminationEllipse
from bvidfe.solver.assembler import assemble_global_stiffness
from bvidfe.solver.boundary import (
    BoundaryCondition,
    apply_dirichlet_penalty,
    mid_surface_z,
    uniaxial_x_bcs,
)
from bvidfe.solver.static import solve_linear_static


def test_apply_dirichlet_penalty_enforces_prescribed_value():
    # Small 3-DOF system: u = [u0, u1, u2]
    K = sp.csc_matrix(np.array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]], dtype=float))
    F = np.zeros(3)
    bcs = [BoundaryCondition(dof=0, value=1.5), BoundaryCondition(dof=2, value=-0.5)]
    Kmod, Fmod = apply_dirichlet_penalty(K, F, bcs, penalty=1e10)
    u = sp.linalg.spsolve(Kmod, Fmod)
    assert abs(u[0] - 1.5) < 1e-6
    assert abs(u[2] - (-0.5)) < 1e-6


def test_uniaxial_x_bcs_compression_returns_xmin_and_xmax_bcs():
    # 8 nodes of a unit cube
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
    bcs = uniaxial_x_bcs(node_coords, applied_strain=-0.01)
    # 4 nodes on x_min get u_x=0, 4 nodes on x_max get u_x=-0.01*Lx=-0.01
    xmin_bcs = [b for b in bcs if abs(b.value) < 1e-12]
    xmax_bcs = [b for b in bcs if abs(b.value - (-0.01)) < 1e-6]
    assert len(xmin_bcs) >= 4  # includes symmetry y and z
    assert len(xmax_bcs) == 4


def test_uniaxial_x_bcs_tension_positive_displacement():
    node_coords = np.array(
        [
            [0, 0, 0],
            [2, 0, 0],
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [2, 0, 1],
            [2, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    bcs = uniaxial_x_bcs(node_coords, applied_strain=0.005)
    xmax = [b for b in bcs if abs(b.value - (0.005 * 2)) < 1e-6]
    assert len(xmax) == 4


def _grid(xs, ys, zs):
    """Structured node grid, ordered like the fe3d mesh (z outer, then y, x)."""
    return np.array([[x, y, z] for z in zs for y in ys for x in xs], dtype=float)


def _constrained(bcs, component):
    """{node: value} of the BCs on displacement component 0/1/2 (x/y/z)."""
    return {b.dof // 3: b.value for b in bcs if b.dof % 3 == component}


def test_uniaxial_x_bcs_platens_mid_surface_layer_and_one_y_node():
    """u_x on every loaded-face node, u_z = 0 on the mid-surface node layer
    only (so no face is held against its through-thickness Poisson
    expansion), u_y = 0 at one node, and no DOF constrained twice."""
    coords = _grid(np.linspace(0, 4, 5), np.linspace(0, 4, 5), np.linspace(0, 0.8, 5))
    bcs = uniaxial_x_bcs(coords, applied_strain=-0.01)
    assert len({b.dof for b in bcs}) == len(bcs)
    x, y, z = coords.T

    ux = _constrained(bcs, 0)
    assert set(ux) == set(np.flatnonzero((x == 0) | (x == 4)))
    for n, value in ux.items():
        assert value == pytest.approx(-0.04 if x[n] == 4 else 0.0)

    uz = _constrained(bcs, 2)
    assert set(uz) == set(np.flatnonzero(np.isclose(z, 0.4)))
    assert set(uz.values()) == {0.0}

    uy = _constrained(bcs, 1)
    assert len(uy) == 1
    (node,) = uy
    np.testing.assert_allclose(coords[node], [0.0, 2.0, 0.4])
    assert uy[node] == 0.0


def test_mid_surface_layer_for_odd_element_layer_count():
    """Three element layers: no node layer on the mid-plane (z = 0.3), so
    the lower of the two nearest layers is restrained."""
    zs = np.array([0.0, 0.2, 0.4, 0.6])
    coords = _grid([0.0, 1.0], [0.0, 1.0], zs)
    assert mid_surface_z(coords) == pytest.approx(0.2)
    uz = _constrained(uniaxial_x_bcs(coords, -0.01), 2)
    assert set(uz) == set(np.flatnonzero(np.isclose(coords[:, 2], 0.2)))


def test_boundary_does_not_change_the_bcs():
    """The membrane solve does not use ``boundary``: out-of-plane edge
    supports carry no load in it (see ``uniaxial_x_bcs``)."""
    coords = _grid(np.linspace(0, 2, 3), np.linspace(0, 2, 3), np.linspace(0, 1, 3))
    reference = uniaxial_x_bcs(coords, -0.01)
    for boundary in ("simply_supported", "clamped", "free", "bogus"):
        assert uniaxial_x_bcs(coords, -0.01, boundary=boundary) == reference


def _small_fe3d_model(layup, delaminations=()):
    """A 20 x 10 mm fe3d mesh at 5 mm (4 x 2 elements, one per ply)."""
    from bvidfe.analysis import AnalysisConfig, MeshParams
    from bvidfe.analysis.fe_mesh import build_fe_mesh
    from bvidfe.analysis.fe_tier import _build_elements
    from bvidfe.core.geometry import PanelGeometry
    from bvidfe.core.laminate import Laminate
    from bvidfe.core.material import MATERIAL_LIBRARY
    from bvidfe.damage.state import DamageState

    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=list(layup),
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0),
        tier="fe3d",
        damage=DamageState(delaminations=[], dent_depth_mm=0.0),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=5.0),
    )
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], cfg.layup_deg, cfg.ply_thickness_mm)
    mesh = build_fe_mesh(cfg, DamageState(delaminations=list(delaminations), dent_depth_mm=0.0))
    elements = _build_elements(mesh, lam)
    K = assemble_global_stiffness(elements, mesh.element_dof_maps, mesh.n_dof)
    return mesh, elements, K.toarray()


def _rigid_body_modes(coords):
    """(n_dof, 6) translations and infinitesimal rotations of the node set."""
    n = coords.shape[0]
    modes = np.zeros((3 * n, 6))
    c = coords - coords.mean(axis=0)
    for k in range(3):
        modes[k::3, k] = 1.0
    modes[1::3, 3], modes[2::3, 3] = -c[:, 2], c[:, 1]  # about x
    modes[0::3, 4], modes[2::3, 4] = c[:, 2], -c[:, 0]  # about y
    modes[0::3, 5], modes[1::3, 5] = -c[:, 1], c[:, 0]  # about z
    return modes


def test_bcs_remove_every_rigid_body_mode():
    """Every rigid-body mode moves a constrained DOF, and the reduced
    stiffness is positive definite. Without the single u_y node exactly one
    mode (the y translation) is left."""
    mesh, _, K = _small_fe3d_model([0.0, 45.0, -45.0, 90.0])
    bcs = uniaxial_x_bcs(mesh.node_coords, -0.01)
    dofs = np.array(sorted({b.dof for b in bcs}))
    assert np.linalg.matrix_rank(_rigid_body_modes(mesh.node_coords)[dofs]) == 6

    def reduced_eigenvalues(constrained):
        free = np.setdiff1d(np.arange(mesh.n_dof), constrained)
        return np.linalg.eigvalsh(K[np.ix_(free, free)])

    eig = reduced_eigenvalues(dofs)
    assert eig.min() > 1e-8 * eig.max()
    y_dof = [b.dof for b in bcs if b.dof % 3 == 1]
    eig = reduced_eigenvalues(np.setdiff1d(dofs, y_dof))
    assert np.sum(eig < 1e-10 * eig.max()) == 1


def test_mid_surface_restraint_carries_no_load_on_a_symmetric_laminate():
    """On a laminate symmetric about its mid-plane, damage included (here a
    delamination on the mid-plane interface, and the through-thickness soft
    zone), u_z = 0 on the mid-plane is a symmetry condition: its reactions,
    and the u_y node's, vanish."""
    delam = DelaminationEllipse(1, (10.0, 5.0), 6.0, 3.0, 0.0)
    mesh, elements, K = _small_fe3d_model([0.0, 90.0, 90.0, 0.0], [delam])
    assert mesh.in_plane_damage_factors.min() < 1.0
    bcs = uniaxial_x_bcs(mesh.node_coords, -0.01)
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs)
    reaction = K @ u
    x_max = mesh.node_coords[:, 0] == mesh.node_coords[:, 0].max()
    axial = abs(reaction[3 * np.flatnonzero(x_max)].sum())
    support = [b.dof for b in bcs if b.dof % 3 != 0]
    assert np.abs(reaction[support]).max() < 1e-9 * axial
