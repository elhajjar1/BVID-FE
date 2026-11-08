import warnings

import numpy as np
import pytest

from bvidfe.analysis.config import MeshParams
from bvidfe.analysis.fe_mesh import FiberBreakCoreUnmeshedWarning, build_fe_mesh
from bvidfe.damage.state import DamageState, DelaminationEllipse


def _simple_config():
    from bvidfe.analysis.config import AnalysisConfig
    from bvidfe.core.geometry import PanelGeometry
    from bvidfe.impact.mapping import ImpactEvent
    from bvidfe.core.geometry import ImpactorGeometry

    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],  # 4 plies
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0),
        loading="compression",
        tier="fe3d",
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=2, in_plane_size_mm=2.0),
    )


def test_mesh_node_and_element_counts_are_consistent():
    cfg = _simple_config()
    mesh = build_fe_mesh(cfg, DamageState([], dent_depth_mm=0.0))
    # nx = 20/2 = 10, ny = 10/2 = 5, nz = 4 plies * 2 = 8
    assert mesh.n_elements == 10 * 5 * 8
    assert mesh.n_nodes == 11 * 6 * 9
    # DOF count = 3 * nodes
    assert mesh.n_dof == 3 * mesh.n_nodes


def test_mesh_ply_assignment_consistent():
    cfg = _simple_config()
    mesh = build_fe_mesh(cfg, DamageState([], dent_depth_mm=0.0))
    # All elements should have ply_index in [0, n_plies-1]
    assert mesh.ply_indices.min() == 0
    assert mesh.ply_indices.max() == 3


def test_mesh_damage_factor_defaults_to_one_for_pristine():
    cfg = _simple_config()
    mesh = build_fe_mesh(cfg, DamageState([], dent_depth_mm=0.0))
    assert np.all(mesh.damage_factors == 1.0)
    assert np.all(mesh.in_plane_damage_factors == 1.0)
    assert np.all(mesh.fiber_tension_factors == 1.0)


def test_mesh_damage_factor_reduced_inside_ellipse():
    cfg = _simple_config()
    # Damage at interface 1 (between ply 1 and ply 2), large ellipse covering panel center
    ds = DamageState(
        [DelaminationEllipse(1, (10, 5), 8, 4, 0)],
        dent_depth_mm=0.3,
    )
    mesh = build_fe_mesh(cfg, ds)
    # Some elements should have reduced stiffness
    n_damaged = (mesh.damage_factors < 1.0).sum()
    assert n_damaged > 0


def test_mesh_delamination_only_softens_damage_zone_in_plane():
    """Pure delamination (no fiber-break radius): OOP reduced at the
    delaminated interface, in-plane reduced to the damage-zone value in the
    footprint (area-weighted where an element straddles its edge, see
    test_fe3d_inplane_damage.py), 1.0 elsewhere."""
    from bvidfe.analysis.fe_mesh import DAMAGE_OOP_FACTOR, DAMAGE_ZONE_INPLANE_FACTOR

    cfg = _simple_config()
    ds = DamageState(
        [DelaminationEllipse(1, (10, 5), 8, 4, 0)],
        dent_depth_mm=0.3,
        fiber_break_radius_mm=0.0,
    )
    mesh = build_fe_mesh(cfg, ds)
    damaged = mesh.damage_factors < 1.0
    assert damaged.any(), "expected at least one damaged element"
    assert np.allclose(mesh.damage_factors[damaged], DAMAGE_OOP_FACTOR)
    softened = mesh.in_plane_damage_factors < 1.0
    assert softened[damaged].all(), "delaminated elements lie in the damage zone"
    assert np.all(mesh.in_plane_damage_factors >= DAMAGE_ZONE_INPLANE_FACTOR - 1e-12)
    assert np.isclose(mesh.in_plane_damage_factors, DAMAGE_ZONE_INPLANE_FACTOR).any()
    assert not softened.all(), "elements outside the footprint keep full in-plane stiffness"
    # Delaminations and matrix cracks cut no fibres: under tension the zone
    # keeps its fibre-direction stiffness.
    assert np.all(mesh.fiber_tension_factors == 1.0)


def test_mesh_damage_zone_takes_the_material_factor():
    """The damage zone's in-plane factor is the material's
    fe3d_zone_inplane_factor (a preset by name or a custom card);
    apply_damage defaults to DAMAGE_ZONE_INPLANE_FACTOR."""
    import dataclasses

    from bvidfe.analysis.fe_mesh import (
        DAMAGE_ZONE_INPLANE_FACTOR,
        apply_damage,
        build_fe_mesh_skeleton,
    )
    from bvidfe.core.material import MATERIAL_LIBRARY

    ds = DamageState([DelaminationEllipse(1, (10, 5), 8, 4, 0)], dent_depth_mm=0.3)
    custom = dataclasses.replace(
        MATERIAL_LIBRARY["IM7/8552"], name="t", fe3d_zone_inplane_factor=0.4
    )
    for material, f in (("IM7/8552", 0.15), ("AS4/8552", 0.219), (custom, 0.4)):
        cfg = dataclasses.replace(_simple_config(), material=material)
        factors = build_fe_mesh(cfg, ds).in_plane_damage_factors
        assert factors.min() == pytest.approx(f)
        assert np.all((factors >= f - 1e-12) & (factors <= 1.0))
    skeleton = build_fe_mesh_skeleton(_simple_config())
    assert apply_damage(skeleton, ds).in_plane_damage_factors.min() == pytest.approx(
        DAMAGE_ZONE_INPLANE_FACTOR
    )
    for bad in (0.0, 1.5):
        with pytest.raises(ValueError, match="zone_inplane_factor"):
            apply_damage(skeleton, ds, bad)


@pytest.mark.parametrize("interface", [0, 1, 2])
def test_oop_factor_takes_the_two_plies_bounding_the_interface(interface):
    """Issue #167: the interface lies on an element face, so selecting the
    elements whose z-extent reaches it took one element layer on each side,
    two plies at elements_per_ply = 1 but one ply's thickness at 2. The two
    bounding plies are softened at any through-thickness refinement."""
    import dataclasses

    from bvidfe.analysis.fe_mesh import DAMAGE_OOP_FACTOR

    ds = DamageState([DelaminationEllipse(interface, (10, 5), 8, 4, 0)], dent_depth_mm=0.3)
    thickness = {}
    for epp in (1, 2, 3):
        cfg = dataclasses.replace(
            _simple_config(), mesh=MeshParams(elements_per_ply=epp, in_plane_size_mm=2.0)
        )
        mesh = build_fe_mesh(cfg, ds)
        hit = mesh.damage_factors < 1.0
        assert np.all(mesh.damage_factors[hit] == DAMAGE_OOP_FACTOR)
        assert set(mesh.ply_indices[hit].tolist()) == {interface, interface + 1}
        z = mesh.node_coords[mesh.element_connectivity[hit]][:, :, 2]
        columns = len({tuple(c) for c in mesh.node_coords[mesh.element_connectivity[hit, 0], :2]})
        thickness[epp] = float((z.max(axis=1) - z.min(axis=1)).sum()) / columns
    # 2 plies x 0.2 mm under each softened column, whatever the refinement.
    assert thickness[1] == pytest.approx(0.4)
    assert thickness[2] == pytest.approx(0.4)
    assert thickness[3] == pytest.approx(0.4)


def test_mesh_fiber_break_core_reduces_in_plane_factor():
    """Inside the fiber-break radius, both factors are reduced; the fibre
    tension factor is not (the core's broken fibres are fe3d TAI's hole,
    issue #164)."""
    from bvidfe.analysis.fe_mesh import (
        DAMAGE_FIBER_BREAK_INPLANE_FACTOR,
        DAMAGE_OOP_FACTOR,
    )

    cfg = _simple_config()
    ds = DamageState(
        [DelaminationEllipse(1, (10, 5), 8, 4, 0)],
        dent_depth_mm=0.3,
        fiber_break_radius_mm=3.0,  # carve out a fiber-break core at the centroid
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", FiberBreakCoreUnmeshedWarning)
        mesh = build_fe_mesh(cfg, ds)
    # Identified by distance: the surrounding damage zone is softened too.
    centroids = np.array([mesh.node_coords[c].mean(axis=0) for c in mesh.element_connectivity])
    fiber_break = np.hypot(centroids[:, 0] - 10.0, centroids[:, 1] - 5.0) <= 3.0
    assert fiber_break.any(), "expected at least one fiber-break-core element"
    assert np.allclose(
        mesh.in_plane_damage_factors[fiber_break],
        DAMAGE_FIBER_BREAK_INPLANE_FACTOR,
    )
    # Fiber-break-core elements must also have OOP reduced
    assert np.all(mesh.damage_factors[fiber_break] == DAMAGE_OOP_FACTOR)
    # The broken fibres are not a softened core: fe3d TAI carries them as
    # its sub-mesh hole, whose radius is at least the core's.
    assert np.all(mesh.fiber_tension_factors == 1.0)


def test_mesh_fiber_break_core_smaller_than_an_element_warns():
    """A core that holds no element centre (the centroid sits on an
    element edge here, 1 mm from the nearest centres of the 2 mm mesh)
    softens nothing and says so; the damage zone around it is unchanged."""
    cfg = _simple_config()
    ell = DelaminationEllipse(1, (10, 5), 8, 4, 0)
    with pytest.warns(FiberBreakCoreUnmeshedWarning, match="0.90 mm.*2 mm mesh.*1.3 mm"):
        mesh = build_fe_mesh(cfg, DamageState([ell], dent_depth_mm=0.3, fiber_break_radius_mm=0.9))
    no_core = build_fe_mesh(cfg, DamageState([ell], dent_depth_mm=0.3))
    np.testing.assert_array_equal(mesh.in_plane_damage_factors, no_core.in_plane_damage_factors)
    np.testing.assert_array_equal(mesh.damage_factors, no_core.damage_factors)
    assert np.all(mesh.fiber_tension_factors == 1.0)
    # A core with no delamination has no centre to sit on: no warning.
    with warnings.catch_warnings():
        warnings.simplefilter("error", FiberBreakCoreUnmeshedWarning)
        build_fe_mesh(cfg, DamageState([], dent_depth_mm=0.0, fiber_break_radius_mm=0.9))


def test_mesh_element_dof_maps_cover_24_dofs_each():
    cfg = _simple_config()
    mesh = build_fe_mesh(cfg, DamageState([], dent_depth_mm=0.0))
    for dof_map in mesh.element_dof_maps:
        assert len(dof_map) == 24
        assert all(0 <= d < mesh.n_dof for d in dof_map)


def test_mesh_node_coordinates_span_panel():
    cfg = _simple_config()
    mesh = build_fe_mesh(cfg, DamageState([], dent_depth_mm=0.0))
    assert mesh.node_coords[:, 0].min() == 0.0
    assert abs(mesh.node_coords[:, 0].max() - 20.0) < 1e-9
    assert abs(mesh.node_coords[:, 1].max() - 10.0) < 1e-9


def test_estimate_fe_mesh_size_returns_sensible_counts():
    from bvidfe.analysis.config import AnalysisConfig, MeshParams
    from bvidfe.analysis.fe_mesh import estimate_fe_mesh_size
    from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
    from bvidfe.impact.mapping import ImpactEvent

    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20, 10),
        loading="compression",
        tier="fe3d",
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=2, in_plane_size_mm=2.0),
    )
    stats = estimate_fe_mesh_size(cfg)
    assert stats["n_elements"] == 10 * 5 * 8
    assert stats["n_dof"] == 3 * 11 * 6 * 9


def test_default_mesh_params_are_conservative():
    from bvidfe.analysis.config import MeshParams

    mp = MeshParams()
    # At these defaults a 150x100 panel with 8 plies must stay below 10k elements
    # (so non-expert users don't hang their GUI).
    assert mp.elements_per_ply == 1
    assert mp.in_plane_size_mm >= 5.0
