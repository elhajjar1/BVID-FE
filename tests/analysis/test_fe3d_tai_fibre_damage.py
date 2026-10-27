"""fe3d TAI must see the impact's broken fibres.

fe3d's tension damage zone keeps its fibres (delaminations and matrix
cracks cut none), and the fibre-break core is off in every material
preset, so fe3d TAI counted only the matrix-cracked zone: nearly flat in
energy and above the closed-form TAI from about 15 J (8-ply QI IM7/8552
150 x 100 mm, 5 mm mesh: 0.928 / 0.912 / 0.894 at 5 / 15 / 30 J against
0.892 / 0.771 / 0.684).

The fibre damage is now the closed-form tiers' equivalent hole
(``tai_hole_radius_mm``, 1.4-3.6 mm across here), the same size and the
same Whitney-Nuismer point-stress notch factor. The hole is far below the
mesh (elements of 2-5 mm, failure index averaged over 2.5 mm, against
d0 = 1 mm), so it is not meshed: its notch factor multiplies the stress the
failure criterion sees in the elements under it. Compression is unchanged.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.analysis import fe_tier as ft
from bvidfe.analysis.fe_mesh import build_fe_mesh, fibre_damage_hole_elements
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.soutis_openhole import (
    fibre_damage_tai,
    lekhnitskii_kt_infinity,
    tai_hole_radius_mm,
)
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage

sys.path.insert(0, str(Path(__file__).parent))
from test_fe3d_inplane_damage import _QI8, _impact_cfg  # noqa: E402

_IM7 = MATERIAL_LIBRARY["IM7/8552"]
_LAM = Laminate(_IM7, _QI8, 0.152)


def _cfg(loading: str, energy_J: float, size_mm: float = 10.0, tier: str = "fe3d"):
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
        loading=loading,
        tier=tier,
        impact=ImpactEvent(energy_J, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=size_mm) if tier == "fe3d" else None,
    )


def test_hole_elements_are_those_under_the_hole():
    """Every ply of each element whose footprint overlaps the hole: the
    4 columns around the impact point (a mesh node at 5 mm) for a hole
    smaller than an element, none for no hole."""
    cfg = _cfg("tension", 15.0, size_mm=5.0)
    damage = impact_to_damage(cfg.impact, _LAM, cfg.panel)
    mesh = build_fe_mesh(cfg, damage)
    R = tai_hole_radius_mm(_IM7, damage, 100.0)
    assert 1.0 < R < 2.5
    hole = fibre_damage_hole_elements(mesh, damage, R)
    assert hole.sum() == 4 * len(_QI8)
    centroids = np.array([mesh.node_coords[c].mean(axis=0) for c in mesh.element_connectivity])
    np.testing.assert_allclose(np.abs(centroids[hole, :2] - (75.0, 50.0)), 2.5)
    assert set(mesh.ply_indices[hole]) == set(range(len(_QI8)))
    assert not fibre_damage_hole_elements(mesh, damage, 0.0).any()
    # A hole larger than an element reaches the next ring.
    assert fibre_damage_hole_elements(mesh, damage, 6.0).sum() > 4 * len(_QI8)


def test_tai_notch_is_the_closed_form_hole():
    """Notch factor = 1 / the closed-form fibre-damage knockdown, at the
    hole elements; 1 elsewhere; none without damage."""
    cfg = _cfg("tension", 15.0)
    damage = impact_to_damage(cfg.impact, _LAM, cfg.panel)
    mesh = build_fe_mesh(cfg, damage)
    notch = ft._tai_fibre_damage_notch(cfg, damage, _LAM, mesh)
    kd = fibre_damage_tai(_IM7, damage, 1.0, lekhnitskii_kt_infinity(_LAM), 100.0)
    hole = fibre_damage_hole_elements(mesh, damage, tai_hole_radius_mm(_IM7, damage, 100.0))
    np.testing.assert_allclose(notch, np.where(hole, 1.0 / kd, 1.0))
    assert ft._tai_fibre_damage_notch(cfg, DamageState(), _LAM, mesh) is None


def test_explicit_fibre_break_core_is_meshed_not_notched():
    """A damage state that carries its own fibre-break core keeps the meshed
    core (fiber_tension_factors) and gets no extra notch."""
    ell = DelaminationEllipse(3, (75.0, 50.0), 20.0, 12.0, 0.0)
    damage = DamageState([ell], dent_depth_mm=0.3, fiber_break_radius_mm=6.0)
    cfg = _cfg("tension", 15.0)
    mesh = build_fe_mesh(cfg, damage)
    assert (mesh.fiber_tension_factors < 1.0).any()
    assert ft._tai_fibre_damage_notch(cfg, damage, _LAM, mesh) is None


def test_fe3d_tai_with_only_the_hole_is_the_closed_form(monkeypatch):
    """Without the damage zone (mesh factors pristine) fe3d sees only the
    hole, under the uniform far-field stress, and returns the closed-form
    fibre-damage knockdown."""
    real_build = ft.build_fe_mesh

    def pristine_factors(config, damage):
        mesh = real_build(config, damage)
        mesh.damage_factors[:] = 1.0
        mesh.in_plane_damage_factors[:] = 1.0
        mesh.fiber_tension_factors[:] = 1.0
        return mesh

    monkeypatch.setattr(ft, "build_fe_mesh", pristine_factors)
    for energy_J in (5.0, 30.0):
        fe3d = BvidAnalysis(_cfg("tension", energy_J)).run().knockdown
        closed = BvidAnalysis(_cfg("tension", energy_J, tier="semi_analytical")).run().knockdown
        assert fe3d == pytest.approx(closed, rel=0.01), energy_J


def test_fe3d_tai_falls_with_energy_with_the_closed_form():
    """Broken fibres set tension after impact: fe3d TAI falls with energy
    like the closed form (0.892 / 0.771 / 0.684 here) and, adding the
    matrix-cracked zone around the hole, stays at or below it. Without the
    hole it was 0.93 / 0.91 / 0.89 on this mesh."""
    fe3d = [BvidAnalysis(_cfg("tension", E)).run().knockdown for E in (5.0, 15.0, 30.0)]
    closed = [
        BvidAnalysis(_cfg("tension", E, tier="semi_analytical")).run().knockdown
        for E in (5.0, 15.0, 30.0)
    ]
    assert fe3d[2] < fe3d[1] < fe3d[0] < 1.0, fe3d
    assert fe3d[0] - fe3d[2] > 0.15, fe3d
    for f, c in zip(fe3d, closed):
        assert c - 0.1 < f <= c + 1e-9, (fe3d, closed)


# fe3d CAI on the code before the TAI fibre damage (10 mm mesh).
_CAI_KNOCKDOWN_BEFORE = {5.0: 0.490301045563204, 15.0: 0.301015353370238, 30.0: 0.23529319659770684}


@pytest.mark.parametrize("energy_J", [5.0, 15.0, 30.0])
def test_fe3d_cai_is_unchanged_by_the_tai_fibre_damage(energy_J):
    r = BvidAnalysis(_impact_cfg("compression", energy_J)).run()
    assert r.knockdown == pytest.approx(_CAI_KNOCKDOWN_BEFORE[energy_J], rel=1e-12)


def test_fe3d_tai_notes_the_hole():
    r = BvidAnalysis(_cfg("tension", 15.0)).run()
    assert any(n.startswith("fe3d TAI: fibre damage") for n in r.notes), r.notes
    assert any(n.startswith("TAI: equivalent hole") for n in r.notes), r.notes


def test_fe3d_tai_strength_is_capped_and_finite():
    damage = impact_to_damage(
        ImpactEvent(15.0, ImpactorGeometry(), 5.5), _LAM, PanelGeometry(150, 100)
    )
    s = ft.fe3d_tai(_cfg("tension", 15.0), damage, _LAM, 400.0)
    assert math.isfinite(s) and 0.0 < s <= 400.0
