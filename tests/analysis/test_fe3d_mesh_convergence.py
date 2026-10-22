"""fe3d first-ply failure converges under mesh refinement.

With the damage zone assigned per element centroid (binary) and failure
checked at every Gauss point, the cross-ply knockdown moved ~25% between
5 and 2.5 mm meshes: elements dropped in and out of the zone, and the
stress at the zone's re-entrant corners is singular, so the pointwise
failure load has no mesh-independent limit. The zone is now weighted by
element area fraction and the failure index averaged over a disk of
``FPF_AVERAGING_RADIUS_MM`` (``fe_mesh._zone_area_fraction``,
``fe_tier._critical_strain_multiplier``).
"""

from __future__ import annotations

import numpy as np
import pytest

from bvidfe import AnalysisConfig, MeshParams
from bvidfe.analysis import fe_tier as ft
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse

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
