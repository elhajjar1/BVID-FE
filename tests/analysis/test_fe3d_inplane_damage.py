"""fe3d first-ply failure must see impact damage.

Delaminations used to soften only the out-of-plane stiffness, and the
fibre-break core is disabled in every material preset, so under in-plane
load the damaged and undamaged stress fields were the same and the fe3d
tension knockdown was 1.000 at 5-30 J. The damage zone (the projected
footprint of the delaminations, where matrix cracks, delaminations and
fibre damage overlap) is now a soft inclusion: in-plane stiffness x0.3
through the thickness.

Damage-zone elements fail at the strain intact material would (strain
equivalence): the criterion sees their effective stress, so their strength
drops with their stiffness. Keeping full strengths made a larger zone look
stronger, up to an intact strength for a fully softened panel.

Under tension the zone keeps its fibre-direction stiffness, and so its
fibre tension strength: delaminations and matrix cracks cut no fibres.
With the whole in-plane stiffness softened in tension too, the zone failed
as if its fibres were cut (8-ply QI 15 J knockdown 0.24 at 5 mm, about one
third of the closed-form TAI; 0.335 against 0.889 measured on Korbelin et
al.'s coupons). Compression is unchanged.

The failure stress was also the failure strain times the pristine CLT
modulus, so a softened panel was credited with load it cannot carry (a
uniformly softened panel reported 3.3x its intact strength). It is now the
reaction force on the loaded edge at failure over the gross area.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.analysis import fe_tier as ft
from bvidfe.analysis.fe_mesh import (
    DAMAGE_OOP_FACTOR,
    DAMAGE_ZONE_INPLANE_FACTOR,
    build_fe_mesh,
)
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import ImpactEvent

_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
_UNDAMAGED = DamageState(delaminations=[], dent_depth_mm=0.0)


def _small_cfg(loading: str = "tension", damage: DamageState = _UNDAMAGED) -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0),
        loading=loading,
        tier="fe3d",
        damage=damage,
        mesh=MeshParams(elements_per_ply=2, in_plane_size_mm=2.0),
    )


def test_damage_zone_softens_in_plane_stiffness_through_the_thickness():
    ell = DelaminationEllipse(1, (10.0, 5.0), 6.0, 3.0, 0.0)
    mesh = build_fe_mesh(_small_cfg(), DamageState([ell], dent_depth_mm=0.3))
    centroids = np.array([mesh.node_coords[c].mean(axis=0) for c in mesh.element_connectivity])
    inside = ((centroids[:, 0] - 10.0) / 6.0) ** 2 + ((centroids[:, 1] - 5.0) / 3.0) ** 2 <= 1.0
    assert inside.any() and (~inside).any()
    # Every ply in the footprint's columns, not just the delaminated interface.
    assert set(mesh.ply_indices[inside]) == {0, 1, 2, 3}
    f = mesh.in_plane_damage_factors
    # Elements wholly inside the footprint take the zone factor, those
    # straddling its edge the area-weighted value between it and 1.
    assert f.min() == pytest.approx(DAMAGE_ZONE_INPLANE_FACTOR)
    assert np.all((f >= DAMAGE_ZONE_INPLANE_FACTOR - 1e-12) & (f <= 1.0))
    assert np.all(f[inside] < 1.0)
    assert ((f > DAMAGE_ZONE_INPLANE_FACTOR + 1e-12) & (f < 1.0)).any()
    # Same in-plane factor in every ply of a column.
    assert np.ptp(f.reshape(-1, 5 * 10), axis=0).max() == 0.0
    # Out-of-plane softening stays at the delaminated interface only.
    assert np.all(mesh.damage_factors[mesh.damage_factors < 1.0] == DAMAGE_OOP_FACTOR)
    assert (mesh.damage_factors[inside] == 1.0).any()


def test_damage_zone_factor_is_the_soft_inclusion_value():
    """Fitted to the NCAMP and Lovejoy & Scotti quasi-isotropic and hard CAI
    coupons (see the comment in fe_mesh.py)."""
    assert DAMAGE_ZONE_INPLANE_FACTOR == pytest.approx(0.15)


@pytest.mark.parametrize("loading", ["tension", "compression"])
def test_first_ply_failure_stress_of_uniformly_softened_panel(loading, monkeypatch):
    """A panel softened in-plane x0.3 everywhere (fibres included: a
    fiber-break core over the whole panel in tension) fails at the intact
    failure strain (strain equivalence) and so carries ~0.3x the load.
    Failure strain x pristine modulus with full strengths reported it 3.3x
    stronger; the reaction force with full strengths reported it as strong
    as intact."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 90, 0, 90], 0.2)
    cfg = _small_cfg(loading)
    fpf = ft._fe3d_cai_first_ply_failure if loading == "compression" else ft.fe3d_tai
    intact = fpf(cfg, _UNDAMAGED, lam, math.inf)

    real_build = ft.build_fe_mesh

    def uniformly_softened(config, damage):
        mesh = real_build(config, damage)
        mesh.in_plane_damage_factors[:] = 0.3
        mesh.fiber_tension_factors[:] = 0.3
        return mesh

    monkeypatch.setattr(ft, "build_fe_mesh", uniformly_softened)
    softened = fpf(cfg, _UNDAMAGED, lam, math.inf)
    assert softened / intact == pytest.approx(0.3, abs=0.05)


def _zone_everywhere(monkeypatch, factor: float = DAMAGE_ZONE_INPLANE_FACTOR) -> None:
    """Make build_fe_mesh put the whole panel inside the damage zone."""
    real_build = ft.build_fe_mesh

    def whole_panel_zone(config, damage):
        mesh = real_build(config, damage)
        mesh.in_plane_damage_factors[:] = factor
        return mesh

    monkeypatch.setattr(ft, "build_fe_mesh", whole_panel_zone)


def test_damage_zone_keeps_fibre_tension_strength(monkeypatch):
    """A [0]_4 coupon inside the damage zone fails in tension where the
    intact one does, at Xt, and carries the same load: the zone's fibres
    are intact. Softening the whole in-plane stiffness in tension (the
    earlier model) made it fail at the intact strain with 0.15x the load."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 0, 0, 0], 0.2)
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 0, 0, 0],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0),
        loading="tension",
        tier="fe3d",
        damage=_UNDAMAGED,
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.0),
    )
    intact = ft.fe3d_tai(cfg, _UNDAMAGED, lam, math.inf)
    assert intact == pytest.approx(MATERIAL_LIBRARY["IM7/8552"].Xt, rel=0.01)
    _zone_everywhere(monkeypatch)
    in_zone = ft.fe3d_tai(cfg, _UNDAMAGED, lam, math.inf)
    assert in_zone / intact == pytest.approx(1.0, abs=0.01)


def test_damage_zone_softens_only_the_matrix_terms_in_tension():
    """Ply-frame compliance of a zone element under tension: S22 and S66
    divided by DAMAGE_ZONE_INPLANE_FACTOR, S11 and the Poisson terms
    pristine. In compression every in-plane stiffness term takes the
    factor, as before tension was told apart."""
    from bvidfe.elements.hex8 import _T_sigma_z

    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], _QI8, 0.152)
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(20.0, 10.0),
        loading="tension",
        tier="fe3d",
        damage=_UNDAMAGED,
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.0),
    )
    ell = DelaminationEllipse(2, (10.0, 5.0), 8.0, 4.0, 0.0)
    mesh = build_fe_mesh(cfg, DamageState([ell], dent_depth_mm=0.2))
    C_pristine = lam.material.get_stiffness_matrix()
    S_pristine = lam.material.get_compliance_matrix()
    f = DAMAGE_ZONE_INPLANE_FACTOR
    S_expected = S_pristine.copy()
    S_expected[1, 1] /= f
    S_expected[5, 5] /= f
    tension = ft._build_elements(mesh, lam, tension=True)
    compression = ft._build_elements(mesh, lam)
    zone = np.flatnonzero(
        np.isclose(mesh.in_plane_damage_factors, f) & (mesh.damage_factors == 1.0)
    )
    assert {float(mesh.ply_angles_deg[e]) for e in zone} == {0.0, 45.0, -45.0, 90.0}
    for e in zone:
        T = _T_sigma_z(np.radians(mesh.ply_angles_deg[e]))
        C_t = T @ tension[e]._C_global @ T.T
        C_c = T @ compression[e]._C_global @ T.T
        S_t = np.linalg.inv(C_t)
        np.testing.assert_allclose(S_t, S_expected, rtol=1e-9, atol=1e-12 * S_pristine[0, 0])
        mask = np.ones((6, 6))
        mask[np.ix_([0, 1, 5], [0, 1, 5])] = f
        np.testing.assert_allclose(C_c, C_pristine * mask, rtol=1e-9, atol=1e-6)
    # Outside the zone both loadings keep the pristine element untouched.
    outside = np.flatnonzero(mesh.in_plane_damage_factors == 1.0)
    assert outside.size
    for e in outside:
        np.testing.assert_array_equal(tension[e]._C_global, compression[e]._C_global)


# Values on the code before tension was told apart (the soft zone's whole
# in-plane stiffness softened in both loadings), same small panel and
# delamination. Compression must reproduce them; tension now keeps the
# zone's fibres.
_SMALL_ELLIPSE = DelaminationEllipse(1, (10.0, 5.0), 6.0, 3.0, 0.0)
_CAI_MULTIPLIER_BEFORE = 0.09863416263950059
_CAI_FPF_BEFORE_MPA = 272.2799510816287
_CAI_KNOCKDOWN_BEFORE = 0.3260979695149508
_TAI_KNOCKDOWN_BEFORE = 0.3157040508823482


def test_compression_path_is_unchanged_by_the_tension_zone_model():
    damage = DamageState([_SMALL_ELLIPSE], dent_depth_mm=0.3)
    cfg = _small_cfg("compression", damage)
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 90, 0, 90], 0.2)
    mesh, elements, _ = ft._fe3d_preflight(cfg, damage, lam)
    u_ref = ft._reference_displacement(cfg, mesh, elements, -1, 0.05)
    c = ft._critical_strain_multiplier(cfg, mesh, elements, u_ref, "larc05")
    assert c == pytest.approx(_CAI_MULTIPLIER_BEFORE, rel=1e-9)
    assert ft._fe3d_cai_first_ply_failure(cfg, damage, lam, math.inf) == pytest.approx(
        _CAI_FPF_BEFORE_MPA, rel=1e-9
    )
    assert BvidAnalysis(cfg).run().knockdown == pytest.approx(_CAI_KNOCKDOWN_BEFORE, rel=1e-9)


def test_tension_knockdown_of_small_panel_rises_with_intact_zone_fibres():
    """0.316 when the zone failed as if its fibres were cut; 0.968 with
    them intact, the 90 deg plies' matrix cracking governing; 0.964 with
    the closed-form fibre-damage hole (0.3 mm across here) at the impact
    site (test_fe3d_tai_fibre_damage.py)."""
    r = BvidAnalysis(_small_cfg("tension", DamageState([_SMALL_ELLIPSE], 0.3))).run()
    assert r.knockdown > _TAI_KNOCKDOWN_BEFORE
    assert r.knockdown == pytest.approx(0.9642, abs=2e-3)


@pytest.mark.parametrize("f_oop", [1.0, DAMAGE_OOP_FACTOR])
def test_damage_zone_failure_stress_is_the_effective_stress(f_oop):
    """In-plane components: what the pristine ply carries at the same strain.
    Out-of-plane components: the actual (damaged) tractions."""
    from bvidfe.elements.hex8 import Hex8Element

    corners = np.array(
        [
            [0, 0, 0],
            [2, 0, 0],
            [2, 2, 0],
            [0, 2, 0],
            [0, 0, 0.2],
            [2, 0, 0.2],
            [2, 2, 0.2],
            [0, 2, 0.2],
        ],
        dtype=float,
    )
    pristine = Hex8Element(corners, MATERIAL_LIBRARY["IM7/8552"], ply_angle_deg=45.0)
    damaged = Hex8Element(corners, MATERIAL_LIBRARY["IM7/8552"], ply_angle_deg=45.0)
    mask = np.full((6, 6), f_oop)
    mask[np.ix_([0, 1, 5], [0, 1, 5])] = DAMAGE_ZONE_INPLANE_FACTOR
    damaged._C_global = damaged._C_global * mask
    u = np.random.default_rng(0).normal(scale=1e-3, size=24)

    effective = ft._failure_stress_material(damaged, u, DAMAGE_ZONE_INPLANE_FACTOR)
    in_plane, out_of_plane = [0, 1, 5], [2, 3, 4]
    np.testing.assert_allclose(
        effective[:, in_plane], pristine.stress_at_gauss_points_material(u)[:, in_plane], rtol=1e-12
    )
    np.testing.assert_allclose(
        effective[:, out_of_plane],
        damaged.stress_at_gauss_points_material(u)[:, out_of_plane],
        rtol=1e-12,
    )
    # Undamaged elements see their actual stress.
    np.testing.assert_array_equal(
        ft._failure_stress_material(pristine, u, 1.0), pristine.stress_at_gauss_points_material(u)
    )


def test_first_ply_failure_stress_of_intact_panel_matches_clt_modulus():
    """For an undamaged panel the reaction-force stress stays close to the
    failure strain times the CLT modulus it replaces."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 90, 0, 90], 0.2)
    cfg = _small_cfg("tension")
    mesh, elements, _ = ft._fe3d_preflight(cfg, _UNDAMAGED, lam)
    strain = ft._solve_failure_strain_analytic(cfg, mesh, elements, +1, "tsai_wu")
    Ex = lam.effective_engineering_constants()[0]
    assert ft.fe3d_tai(cfg, _UNDAMAGED, lam, math.inf) == pytest.approx(strain * Ex, rel=0.1)


def _impact_cfg(loading: str, energy_J: float) -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        loading=loading,
        tier="fe3d",
        impact=ImpactEvent(energy_J, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0),
    )


def test_fe3d_tension_knockdown_responds_to_impact_damage():
    """Was 1.000: first-ply failure sat at a panel corner in both runs.
    With the zone's fibres cut in tension too it was 0.286 on this mesh;
    with them intact the matrix-cracked zone cost about 9% (0.909); the
    closed-form fibre-damage hole at the impact site brings it to 0.714
    (closed form 0.771, see test_fe3d_tai_fibre_damage.py)."""
    r = BvidAnalysis(_impact_cfg("tension", 15.0)).run()
    assert r.knockdown == pytest.approx(0.7144, abs=2e-3), r.notes


def test_fe3d_tension_knockdown_falls_with_impact_energy():
    kd = [BvidAnalysis(_impact_cfg("tension", E)).run().knockdown for E in (5.0, 15.0)]
    assert kd[1] < kd[0] < 1.0, kd


def test_fe3d_tension_knockdown_does_not_rise_as_the_damage_zone_spans_the_panel():
    """With full strengths in the damage zone, a zone spanning most of the
    panel lost its stress concentration and the panel looked stronger:
    cross-ply tension 0.49 at 5 J -> 0.64 at 15 J on this mesh."""
    cross_ply = [0, 90, 0, 90, 90, 0, 90, 0]
    kd = []
    for energy_J in (5.0, 15.0):
        cfg = AnalysisConfig(
            material="IM7/8552",
            layup_deg=cross_ply,
            ply_thickness_mm=0.152,
            panel=PanelGeometry(150.0, 100.0),
            loading="tension",
            tier="fe3d",
            impact=ImpactEvent(energy_J, ImpactorGeometry(), mass_kg=5.5),
            mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=10.0),
        )
        kd.append(BvidAnalysis(cfg).run().knockdown)
    assert kd[1] <= kd[0] < 1.0, kd
