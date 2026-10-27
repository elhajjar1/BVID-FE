"""fe3d first-ply failure evaluates the criteria in the ply material frame.

Unidirectional coupons under uniaxial x-strain have a uniform, closed-form
stress state, so the failure strain is known exactly for every criterion:

* ``[0]_4``  in compression: s1 = E11 * eps fails at s1 = -Xc  ->  eps = Xc / E11
* ``[90]_4`` in compression: s2 = E22 * eps fails at s2 = -Yc  ->  eps = Yc / E22

(and Xt / E11, Yt / E22 in tension). Tsai-Wu, LaRC05 and Puck all reduce to
these uniaxial strengths. Before the frame fix, fe3d passed global-frame
stress to the criteria, so a 90 deg ply's transverse stress was compared
with Xc and ``[90]_4`` never failed below the 5% strain cap.

The uniform state also needs boundary conditions that let the plies expand
freely by Poisson: with ``u_z = 0`` on the whole loaded faces (the
previous simply-supported and clamped restraints) and the whole ``z_min``
face, ``[90]_4`` failed 3-4% (LaRC05 under clamped: 10.7%) early in
compression and ``[0]_4`` up to 1.7% early under Tsai-Wu, from the
through-thickness stress next to the edges (``nu23`` = 0.5 for IM7/8552).
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


def _pristine_setup(ply_angle: float, boundary: str = "simply_supported"):
    layup = [ply_angle] * 4
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=layup,
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0, boundary=boundary),
        loading="compression",
        tier="fe3d",
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.5),
    )
    lam = Laminate(material=_MATERIAL, layup_deg=layup, ply_thickness_mm=0.2)
    mesh = build_fe_mesh(cfg, DamageState(delaminations=[], dent_depth_mm=0.0))
    return cfg, mesh, _build_elements(mesh, lam)


_BOUNDARIES = ["simply_supported", "clamped", "free"]
_CRITERIA = ["tsai_wu", "larc05", "puck"]
# Puck's inter-fibre mode C does not reduce exactly to Yc: the index under
# pure transverse compression is 1 - 4e-4 for IM7/8552, a property of the
# criterion, not of the stress field.
_REL = 1e-3


@pytest.mark.parametrize("boundary", _BOUNDARIES)
@pytest.mark.parametrize("criterion", _CRITERIA)
def test_ninety_degree_laminate_fails_at_transverse_compressive_strength(criterion, boundary):
    cfg, mesh, elements = _pristine_setup(90.0, boundary=boundary)
    eps = _solve_failure_strain_analytic(cfg, mesh, elements, strain_sign=-1, criterion=criterion)
    assert eps < 0.05, "never failed below the strain cap: transverse stress checked against Xc"
    assert eps == pytest.approx(_MATERIAL.Yc / _MATERIAL.E22, rel=_REL)


@pytest.mark.parametrize("boundary", _BOUNDARIES)
@pytest.mark.parametrize("criterion", _CRITERIA)
def test_zero_degree_laminate_fails_at_fiber_compressive_strength(criterion, boundary):
    cfg, mesh, elements = _pristine_setup(0.0, boundary=boundary)
    eps = _solve_failure_strain_analytic(cfg, mesh, elements, strain_sign=-1, criterion=criterion)
    assert eps == pytest.approx(_MATERIAL.Xc / _MATERIAL.E11, rel=_REL)


@pytest.mark.parametrize("ply_angle, strength, modulus", [(0.0, "Xt", "E11"), (90.0, "Yt", "E22")])
def test_unidirectional_laminate_fails_at_tensile_strength(ply_angle, strength, modulus):
    """Tension (Tsai-Wu, the TAI criterion): Xt / E11 and Yt / E22."""
    cfg, mesh, elements = _pristine_setup(ply_angle)
    eps = _solve_failure_strain_analytic(cfg, mesh, elements, strain_sign=+1, criterion="tsai_wu")
    expected = getattr(_MATERIAL, strength) / getattr(_MATERIAL, modulus)
    assert eps == pytest.approx(expected, rel=_REL)
