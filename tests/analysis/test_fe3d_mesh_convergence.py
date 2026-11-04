"""fe3d first-ply failure converges under mesh refinement.

With the damage zone assigned per element centroid (binary) and failure
checked at every Gauss point, the cross-ply knockdown moved ~25% between
5 and 2.5 mm meshes: elements dropped in and out of the zone, and the
stress at the zone's re-entrant corners is singular, so the pointwise
failure load has no mesh-independent limit. The zone is now weighted by
element area fraction and the failure index averaged over a disk of
``FPF_AVERAGING_RADIUS_MM`` (``fe_mesh._zone_area_fraction``,
``fe_tier._critical_strain_multiplier``).

The crossed-ellipse coupon below is the hard case: its arms are 6 mm
wide, with re-entrant corners where they cross. 2.5 and 1.25 mm agree,
but finer than that the first-ply failure ratio in ``BvidAnalysis``
scatters with no trend (0.481 / 0.447 / 0.476 / 0.494 / 0.496 / 0.487 /
0.463 at 1.67 / 1.5 / 1.36 / 1.25 / 1.15 / 1.07 / 1.0 mm), as the
corners fall differently on each grid.

The h-refinement tests at the end (issue #89, which proposed them as
``tests/analysis/test_fe3d_h_refinement.py``) run ``BvidAnalysis`` itself,
so they cover the knockdown and E11 knockdown users see.
"""

from __future__ import annotations

import numpy as np
import pytest

from bvidfe import AnalysisConfig, AnalysisResults, BvidAnalysis, MeshParams
from bvidfe.analysis import fe_tier as ft
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.impact.mapping import ImpactEvent

_LAYUP = [0, 90, 90, 0]
# Crossed ellipses: their union has re-entrant corners.
_DAMAGE = DamageState(
    [
        DelaminationEllipse(0, (30, 15), 9, 3, 0),
        DelaminationEllipse(1, (30, 15), 3, 9, 90),
        DelaminationEllipse(2, (30, 15), 9, 3, 0),
    ],
    dent_depth_mm=0.2,
)
_LAM = Laminate(MATERIAL_LIBRARY["IM7/8552"], _LAYUP, 0.2)


def _cfg(h: float) -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_LAYUP,
        ply_thickness_mm=0.2,
        panel=PanelGeometry(60.0, 30.0),
        loading="compression",
        tier="fe3d",
        damage=_DAMAGE,
        mesh=MeshParams(in_plane_size_mm=h),
    )


def _fpf_load(cfg: AnalysisConfig, damage: DamageState, radius) -> float:
    mesh = build_fe_mesh(cfg, damage)
    elements = ft._build_elements(mesh, _LAM)
    u = ft._reference_displacement(cfg, mesh, elements, -1, 0.05)
    c = ft._critical_strain_multiplier(cfg, mesh, elements, u, "larc05", radius)
    return c * ft._loaded_edge_force(mesh, elements, u)


def _knockdown(h: float, radius=ft.FPF_AVERAGING_RADIUS_MM) -> float:
    cfg = _cfg(h)
    undamaged = DamageState(delaminations=[], dent_depth_mm=0.0)
    return _fpf_load(cfg, _DAMAGE, radius) / _fpf_load(cfg, undamaged, radius)


def test_cross_ply_knockdown_is_mesh_converged():
    coarse, fine = _knockdown(2.5), _knockdown(1.25)
    assert 0.0 < fine < 1.0
    # 0.523 vs 0.508; with the binary zone and pointwise check, 0.460 vs
    # 0.358 vs 0.448 at 2.5 / 1.875 / 1.25 mm.
    assert fine == pytest.approx(coarse, rel=0.06)


def test_averaging_lowers_the_peak_index():
    """The averaged index never exceeds the peak pointwise one, so the
    averaged failure load is never below the pointwise one."""
    cfg = _cfg(2.5)
    assert _fpf_load(cfg, _DAMAGE, None) < _fpf_load(cfg, _DAMAGE, ft.FPF_AVERAGING_RADIUS_MM)


@pytest.mark.parametrize(
    "a, b, expected",
    [
        (0.0, 4.0, 0.5),  # LaRC05: idx = b c^2
        (2.0, 0.0, 0.5),  # Puck: idx = a c
        (1.0, 2.0, 0.5),  # Tsai-Wu: 2 c^2 + c = 1
        (0.0, -1.0, np.inf),  # never reaches 1
        (-1.0, 0.0, np.inf),
    ],
)
def test_smallest_positive_root(a, b, expected):
    c = ft._smallest_positive_root(np.array([a]), np.array([b]))
    assert c[0] == pytest.approx(expected)


@pytest.mark.parametrize("h", [2.5, 2.0, 1.25, 1.0, 0.5])
def test_disk_kernel_keeps_the_disk_area(h):
    """Weighted by overlap the averaging disk keeps its area on any mesh;
    counting whole elements by their centres shrank it from 31 to 20 mm^2
    between 2.5 and 2 mm."""
    k = ft._disk_kernel(h, h, 2.5)
    assert k.sum() * h * h == pytest.approx(np.pi * 2.5**2, rel=0.03)
    assert k[k.shape[0] // 2, k.shape[1] // 2] == 1.0


def test_disk_kernel_at_the_default_mesh_is_the_element_alone():
    """At 5 mm the 2.5 mm disk fits inside its element, so the default-mesh
    results (and the damage-zone factor fitted on them) do not move."""
    k = ft._disk_kernel(5.0, 5.0, 2.5)
    assert np.count_nonzero(k) == 1 and k[1, 1] > 0


# --- h-refinement of the public results (issue #89) ---------------------
#
# fe3d's strength is a disk-averaged first-ply failure load, a local peak
# that no energy argument bounds, so it need not move monotonically with
# the mesh. Once the zone is resolved it lies in a band set by where the
# zone's edge falls on each grid (its elements are weighted by area
# fraction): within about 4% of the finest mesh at every resolved size
# measured, but with no shrinking step-to-step trend (150 x 100 mm, 25 J,
# compression: 0.2467 / 0.2487 / 0.2353 / 0.2412 / 0.2391 at
# 5 / 3.33 / 3.0 / 2.5 / 2.0 mm). A circular zone (r = 6 mm) in the
# coupon below gives 0.3277 / 0.3254 / 0.3265 / 0.3149 at 2.5 / 1.67 /
# 1.25 / 1.0 mm, 3.6% on the last step, so _PLATEAU holds for the
# ellipse's alignment, not for every zone. The E11 knockdown integrates
# the whole field and in compression falls at each step of both sweeps,
# by less each time, as a displacement-based model softens when refined
# (between meshes that are not refinements of each other, 3.0 and 2.5 mm
# on the crossed-ellipse coupon, it can rise by 0.1%); in tension the
# zone keeps its fibres and E11 barely moves.
#
# Zones only a few elements across at 2.5 mm are not resolved at these
# sizes: a 3 J impact on a 60 x 40 mm [0/90]s coupon (304 mm^2) gives a
# compression knockdown of 0.447 / 0.402 / 0.247 / 0.235 at 5 / 2.5 /
# 1.25 / 1.0 mm, still moving 5% on the last step.

#: Every resolved mesh within this of the finest one (the alignment band).
_BAND = 0.05
#: The two finest meshes within this of each other.
_PLATEAU = 0.02


def _rel_change(coarse: float, fine: float) -> float:
    return abs(coarse - fine) / fine


def _check_h_refinement(
    runs: dict[float, AnalysisResults], loading: str, resolved: list[float]
) -> None:
    """The properties both sweeps share. ``runs`` maps the mesh size
    (coarse to fine) to its result; ``resolved`` lists the sizes that
    resolve the damage zone."""
    hs = list(runs)
    kd = [runs[h].knockdown for h in hs]
    e11 = []
    for h in hs:
        e = runs[h].E11_knockdown
        assert e is not None
        e11.append(e)
    assert all(np.isfinite(k) and 0.0 < k < 1.0 for k in kd)
    assert all(np.isfinite(e) and 0.0 < e < 1.0 for e in e11)
    # Strength: in the band once resolved, the finest pair plateaued.
    for h in resolved:
        assert _rel_change(runs[h].knockdown, kd[-1]) < _BAND, h
    assert _rel_change(kd[-2], kd[-1]) < _PLATEAU
    stiffness = np.array(e11)
    if loading == "compression":
        # Stiffness: falls with every refinement, by less each time.
        steps = -np.diff(stiffness) / stiffness[1:]
        assert np.all(steps > 0.0)
        assert np.all(np.diff(steps) < 0.0)
    else:
        assert stiffness.max() / stiffness.min() - 1.0 < 1e-3


# 4-ply [0/90]s 40 x 30 mm coupon, one elliptical zone (16 x 10 mm,
# 126 mm^2) through its three interfaces: 3 x 2 elements at 5 mm.
_ELLIPSE = DamageState(
    [DelaminationEllipse(i, (20, 15), 8, 5, 0) for i in range(3)],
    dent_depth_mm=0.2,
)
_COUPON_H = (5.0, 2.5, 1.25, 1.0)


@pytest.fixture(scope="module", params=["compression", "tension"])
def coupon_sweep(request) -> tuple[str, dict[float, AnalysisResults]]:
    """``BvidAnalysis`` on the coupon at each of ``_COUPON_H``, about 5 s
    per loading."""
    runs = {
        h: BvidAnalysis(
            AnalysisConfig(
                material="IM7/8552",
                layup_deg=_LAYUP,
                ply_thickness_mm=0.2,
                panel=PanelGeometry(40.0, 30.0),
                loading=request.param,
                tier="fe3d",
                damage=_ELLIPSE,
                mesh=MeshParams(in_plane_size_mm=h),
            )
        ).run()
        for h in _COUPON_H
    }
    return request.param, runs


def test_coupon_knockdown_h_refines(coupon_sweep):
    """At 5 / 2.5 / 1.25 / 1.0 mm the knockdown is 0.4255 / 0.3829 /
    0.3882 / 0.3884 in compression (0.3742 at 1.67 mm) and 0.9556 / 0.9464
    / 0.9503 / 0.9514 in tension; the E11 knockdown 0.8504 / 0.8393 /
    0.8323 / 0.8307 in compression and 0.9952 at every size in tension."""
    loading, runs = coupon_sweep
    _check_h_refinement(runs, loading, resolved=[2.5, 1.25, 1.0])
    # The default 5 mm mesh under-resolves the zone (the averaging disk is
    # the element alone), so its step is the largest.
    kd = [runs[h].knockdown for h in _COUPON_H]
    steps = [_rel_change(a, b) for a, b in zip(kd, kd[1:])]
    assert steps[0] == max(steps)


# The issue's reference config: 8-ply QI IM7/8552 150 x 100 mm, 25 J
# (2,615 mm^2). 3.4 mm gives 3.33 mm elements (45 x 30).
_PANEL_H = (5.0, 3.4, 2.5, 2.0)


@pytest.mark.slow
@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_reference_panel_knockdown_h_refines(loading):
    """At 5 / 3.33 / 2.5 / 2.0 mm (about 75 s per loading) the knockdown
    is 0.2467 / 0.2487 / 0.2412 / 0.2391 in compression and 0.6511 /
    0.6558 / 0.6570 / 0.6566 in tension; the E11 knockdown 0.7025 / 0.6951
    / 0.6910 / 0.6896 in compression and 0.9813 at every size in tension.
    The zone is resolved at 5 mm: the default mesh is within 3.2% of
    2.0 mm."""
    runs = {
        h: BvidAnalysis(
            AnalysisConfig(
                material="IM7/8552",
                layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
                ply_thickness_mm=0.152,
                panel=PanelGeometry(150.0, 100.0),
                loading=loading,
                tier="fe3d",
                impact=ImpactEvent(25.0),
                mesh=MeshParams(in_plane_size_mm=h),
            )
        ).run()
        for h in _PANEL_H
    }
    _check_h_refinement(runs, loading, resolved=list(_PANEL_H))
