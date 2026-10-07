"""fe3d first-ply failure must see impact damage.

Delaminations used to soften only the out-of-plane stiffness, and the
fibre-break core is disabled in every material preset, so under in-plane
load the damaged and undamaged stress fields were the same and the fe3d
tension knockdown was 1.000 at 5-30 J. The damage zone (the projected
footprint of the delaminations, where matrix cracks, delaminations and
fibre damage overlap) is now a soft inclusion: in-plane stiffness x0.3
through the thickness.

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
    assert np.all(mesh.in_plane_damage_factors[inside] == DAMAGE_ZONE_INPLANE_FACTOR)
    assert np.all(mesh.in_plane_damage_factors[~inside] == 1.0)
    # Out-of-plane softening stays at the delaminated interface only.
    assert np.all(mesh.damage_factors[mesh.damage_factors < 1.0] == DAMAGE_OOP_FACTOR)
    assert (mesh.damage_factors[inside] == 1.0).any()


def test_damage_zone_factor_is_the_soft_inclusion_value():
    assert DAMAGE_ZONE_INPLANE_FACTOR == pytest.approx(0.3)


@pytest.mark.parametrize("loading", ["tension", "compression"])
def test_first_ply_failure_stress_counts_the_load_a_softened_panel_carries(loading, monkeypatch):
    """Softening the in-plane stiffness of every element by x0.3, with the
    ply strengths unchanged, raises the strain at failure ~3.3x but leaves
    the load at failure the same. Failure strain x pristine modulus
    reported the softened panel as 3.3x stronger."""
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], [0, 90, 0, 90], 0.2)
    cfg = _small_cfg(loading)
    fpf = ft._fe3d_cai_first_ply_failure if loading == "compression" else ft.fe3d_tai
    intact = fpf(cfg, _UNDAMAGED, lam, math.inf)

    real_build = ft.build_fe_mesh

    def uniformly_softened(config, damage):
        mesh = real_build(config, damage)
        mesh.in_plane_damage_factors[:] = 0.3
        return mesh

    monkeypatch.setattr(ft, "build_fe_mesh", uniformly_softened)
    softened = fpf(cfg, _UNDAMAGED, lam, math.inf)
    assert softened / intact == pytest.approx(1.0, abs=0.05)


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
    """Was 1.000: first-ply failure sat at a panel corner in both runs."""
    r = BvidAnalysis(_impact_cfg("tension", 15.0)).run()
    assert r.knockdown < 0.8, r.notes


def test_fe3d_tension_knockdown_falls_with_impact_energy():
    kd = [BvidAnalysis(_impact_cfg("tension", E)).run().knockdown for E in (5.0, 15.0)]
    assert kd[1] < kd[0] < 1.0, kd
