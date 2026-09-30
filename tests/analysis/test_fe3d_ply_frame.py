"""fe3d first-ply failure evaluates the criteria in the ply material frame.

Unidirectional coupons under uniaxial x-strain have a uniform, closed-form
stress state, so the failure strain is known exactly for every criterion:

* ``[0]_4``  in compression: s1 = E11 * eps fails at s1 = -Xc  ->  eps = Xc / E11
* ``[90]_4`` in compression: s2 = E22 * eps fails at s2 = -Yc  ->  eps = Yc / E22

Tsai-Wu, LaRC05 and Puck all reduce to these uniaxial strengths. Before the
fix, fe3d passed global-frame stress to the criteria, so a 90 deg ply's
transverse stress was compared with Xc and ``[90]_4`` never failed below the
5% strain cap.
"""

from __future__ import annotations

import pytest

from bvidfe.analysis import AnalysisConfig, MeshParams
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.analysis.fe_tier import _build_elements, _solve_failure_strain_analytic
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState
from bvidfe.impact.mapping import ImpactEvent

_MATERIAL = MATERIAL_LIBRARY["IM7/8552"]


def _pristine_setup(ply_angle: float):
    layup = [ply_angle] * 4
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=layup,
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0),
        loading="compression",
        tier="fe3d",
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.5),
    )
    lam = Laminate(material=_MATERIAL, layup_deg=layup, ply_thickness_mm=0.2)
    mesh = build_fe_mesh(cfg, DamageState(delaminations=[], dent_depth_mm=0.0))
    return cfg, mesh, _build_elements(mesh, lam)


@pytest.mark.parametrize("criterion", ["tsai_wu", "larc05", "puck"])
def test_ninety_degree_laminate_fails_at_transverse_compressive_strength(criterion):
    cfg, mesh, elements = _pristine_setup(90.0)
    eps = _solve_failure_strain_analytic(cfg, mesh, elements, strain_sign=-1, criterion=criterion)
    assert eps < 0.05, "never failed below the strain cap: transverse stress checked against Xc"
    assert eps == pytest.approx(_MATERIAL.Yc / _MATERIAL.E22, rel=0.01)


# Tsai-Wu's linear F3 * s3 term (Zt falls back to Yt = 73 MPa) picks up the
# ~5 MPa through-thickness tension that the pinned z_min face and loaded-edge
# u_z restraint induce, lowering the failure strain ~4.7%. LaRC05 and Puck
# ignore small s3 in fiber mode. The effect is the same before and after the
# frame fix, since a 0 deg ply's material frame is the global frame.
_ZERO_DEG_REL_TOL = {"tsai_wu": 0.06, "larc05": 0.01, "puck": 0.01}


@pytest.mark.parametrize("criterion", ["tsai_wu", "larc05", "puck"])
def test_zero_degree_laminate_fails_at_fiber_compressive_strength(criterion):
    cfg, mesh, elements = _pristine_setup(0.0)
    eps = _solve_failure_strain_analytic(cfg, mesh, elements, strain_sign=-1, criterion=criterion)
    assert eps == pytest.approx(_MATERIAL.Xc / _MATERIAL.E11, rel=_ZERO_DEG_REL_TOL[criterion])
