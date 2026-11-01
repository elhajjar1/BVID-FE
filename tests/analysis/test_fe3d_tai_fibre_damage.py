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

A damage state's own fibre-break core (a C-scan's ``fiber_break_radius_mm``)
is part of that hole, which is at least the core's radius, as in the
closed-form tiers. It used to switch the hole off in favour of the meshed
core, which takes elements by their centre and held none at 5 mm for
cores under 3.5 mm: a 1.5 mm core raised the 8-ply QI 15 J knockdown
from 0.714 to 0.912 against 0.727 closed-form (issue #164). The meshed
core keeps its compression factors and no longer softens the fibres in
tension, which would count the hole's broken fibres twice.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import dataclasses

import numpy as np
import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.analysis import fe_tier as ft
from bvidfe.analysis.fe_mesh import (
    DAMAGE_FIBER_BREAK_INPLANE_FACTOR,
    FiberBreakCoreUnmeshedWarning,
    build_fe_mesh,
    fibre_damage_hole_elements,
)
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
_PANEL = PanelGeometry(150.0, 100.0)


def _cfg(loading: str, energy_J: float, size_mm: float = 10.0, tier: str = "fe3d"):
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=_PANEL,
        loading=loading,
        tier=tier,
        impact=ImpactEvent(energy_J, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=size_mm) if tier == "fe3d" else None,
    )


def _cscan(energy_J: float, fiber_break_radius_mm: float) -> DamageState:
    """The impact's damage state as a C-scan would give it, with a measured
    fibre-break core."""
    damage = impact_to_damage(ImpactEvent(energy_J, ImpactorGeometry(), 5.5), _LAM, _PANEL)
    return dataclasses.replace(damage, fiber_break_radius_mm=fiber_break_radius_mm)


def _cscan_cfg(loading: str, damage: DamageState, size_mm: float = 10.0, tier: str = "fe3d"):
    return dataclasses.replace(_cfg(loading, 15.0, size_mm, tier), impact=None, damage=damage)


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


def test_explicit_fibre_break_core_is_the_hole_not_a_softened_core():
    """A damage state that carries its own fibre-break core gets the notch,
    on the elements under a hole of the core's radius (6 mm here, above
    the footprint's hole), and no fibre softening in the meshed core; the
    core's compression factors are as before."""
    ell = DelaminationEllipse(3, (75.0, 50.0), 20.0, 12.0, 0.0)
    damage = DamageState([ell], dent_depth_mm=0.3, fiber_break_radius_mm=6.0)
    cfg = _cfg("tension", 15.0)
    mesh = build_fe_mesh(cfg, damage)
    assert np.all(mesh.fiber_tension_factors == 1.0)
    core = mesh.in_plane_damage_factors == DAMAGE_FIBER_BREAK_INPLANE_FACTOR
    assert core.any() and set(mesh.ply_indices[core]) == set(range(len(_QI8)))
    R, kd = ft.tai_fibre_damage_hole(cfg, damage, _LAM)
    assert (
        R == 6.0 > tai_hole_radius_mm(_IM7, dataclasses.replace(damage, fiber_break_radius_mm=0.0))
    )
    assert 0.0 < kd < 1.0
    notch = ft._tai_fibre_damage_notch(cfg, damage, _LAM, mesh)
    hole = fibre_damage_hole_elements(mesh, damage, 6.0)
    np.testing.assert_allclose(notch, np.where(hole, 1.0 / kd, 1.0))
    # Every meshed core element is under the hole, and the hole reaches
    # further (overlap, not centre-in-disk).
    assert core.sum() < hole.sum() and not (core & ~hole).any()


@pytest.mark.parametrize("fiber_break_radius_mm", [0.0, 1.5, 4.0])
def test_fe3d_hole_is_the_closed_form_hole(fiber_break_radius_mm):
    """fe3d and the closed-form tiers share one hole radius, max(footprint
    hole, fibre-break core), and one knockdown."""
    damage = _cscan(15.0, fiber_break_radius_mm)
    cfg = _cscan_cfg("tension", damage)
    R, kd = ft.tai_fibre_damage_hole(cfg, damage, _LAM)
    assert R == tai_hole_radius_mm(_IM7, damage, 100.0)
    assert R == pytest.approx(max(1.247, fiber_break_radius_mm), abs=0.01)
    assert kd == fibre_damage_tai(_IM7, damage, 1.0, lekhnitskii_kt_infinity(_LAM), 100.0)


def test_sub_element_fibre_break_core_lowers_fe3d_tai_on_the_default_mesh():
    """On the 5 mm mesh a 1.5 mm C-scan core holds no element centre (the
    nearest is 3.54 mm from the impact point). It used to switch the hole
    off and raise TAI from 0.714 to 0.912 (no core / core); now the hole is
    at least the core, so the strength is at or below the no-core one
    (closed form 0.771 -> 0.727)."""
    no_core = _cscan(15.0, 0.0)
    core = _cscan(15.0, 1.5)
    s_no_core = ft.fe3d_tai(_cscan_cfg("tension", no_core, 5.0), no_core, _LAM, math.inf)
    with pytest.warns(FiberBreakCoreUnmeshedWarning):
        s_core = ft.fe3d_tai(_cscan_cfg("tension", core, 5.0), core, _LAM, math.inf)
    assert 0.0 < s_core <= s_no_core, (s_core, s_no_core)
    assert s_core < 0.98 * s_no_core, (s_core, s_no_core)


def test_unmeshed_core_is_recorded_in_the_results():
    r = BvidAnalysis(_cscan_cfg("tension", _cscan(15.0, 1.5))).run()
    assert any(w.startswith("FiberBreakCoreUnmeshedWarning") for w in r.warnings), r.warnings
    assert any(n.startswith("fe3d TAI: fibre damage") for n in r.notes), r.notes


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


# fe3d CAI on the code before the TAI fibre damage (10 mm mesh). The 5 J
# value is recaptured for #165: fe3d CAI takes min(growth, FPF) and, with
# the sublaminate buckling rectangle following the ellipse orientation, the
# growth stress at 5 J fell from 300.5 to 281.8 MPa and now governs there
# (0.490 -> 0.469); 15 and 30 J stay at first-ply failure.
_CAI_KNOCKDOWN_BEFORE = {
    5.0: 0.4691595567346301,
    15.0: 0.301015353370238,
    30.0: 0.23529319659770684,
}


@pytest.mark.parametrize("energy_J", [5.0, 15.0, 30.0])
def test_fe3d_cai_is_unchanged_by_the_tai_fibre_damage(energy_J):
    r = BvidAnalysis(_impact_cfg("compression", energy_J)).run()
    assert r.knockdown == pytest.approx(_CAI_KNOCKDOWN_BEFORE[energy_J], rel=1e-12)


# fe3d CAI with a measured fibre-break core on the code before #164 (10 mm
# mesh, 15 J): a 1.5 mm core holds no element centre, an 8 mm one four columns.
_CAI_CORE_KNOCKDOWN_BEFORE = {1.5: 0.301015353370238, 8.0: 0.2884425155093483}


@pytest.mark.parametrize("fiber_break_radius_mm", [1.5, 8.0])
def test_fe3d_cai_with_a_fibre_break_core_is_unchanged(fiber_break_radius_mm):
    cfg = _cscan_cfg("compression", _cscan(15.0, fiber_break_radius_mm))
    r = BvidAnalysis(cfg).run()
    expected = _CAI_CORE_KNOCKDOWN_BEFORE[fiber_break_radius_mm]
    assert r.knockdown == pytest.approx(expected, rel=1e-12)


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
