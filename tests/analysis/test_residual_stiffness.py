"""Residual in-plane stiffness alongside the strength knockdown (issue #114).

Every tier reports the intact laminate's CLT membrane moduli
(``pristine_E11_MPa`` / ``E22`` / ``G12``). Only fe3d computes a residual:
the E11 knockdown is its damaged over undamaged panel's secant axial
modulus, the loaded-edge reaction force at the same end strain in the two
first-ply failure solves. The closed-form tiers have no damage-zone
stiffness to homogenise and leave the residual fields ``None``.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, DamageState, DelaminationEllipse, MeshParams
from bvidfe.analysis import fe_tier as ft
from bvidfe.analysis.results import AnalysisResults
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.mapping import ImpactEvent

_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
_LAM = Laminate(MATERIAL_LIBRARY["IM7/8552"], _QI8, 0.152)
_UNDAMAGED = DamageState(delaminations=[], dent_depth_mm=0.0)
_DELAMINATED = DamageState(
    delaminations=[
        DelaminationEllipse(2, (20.0, 15.0), 12.0, 8.0, 0.0),
        DelaminationEllipse(3, (20.0, 15.0), 16.0, 12.0, 30.0),
    ],
    dent_depth_mm=0.1,
)
_STIFFNESS_KEYS = [
    f"{prefix}{name}{suffix}"
    for name in ("E11", "E22", "G12")
    for prefix, suffix in (("pristine_", "_MPa"), ("residual_", "_MPa"), ("", "_knockdown"))
]


def _cfg(tier: str, loading: str, damage: DamageState) -> AnalysisConfig:
    """40 x 30 mm panel at a 4 mm mesh: an fe3d run takes about a second."""
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=_QI8,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(40.0, 30.0),
        loading=loading,
        tier=tier,
        damage=damage,
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=4.0),
    )


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("tier", ["empirical", "semi_analytical", "fe3d"])
def test_pristine_moduli_are_the_clt_constants_in_every_tier(tier, loading):
    r = BvidAnalysis(_cfg(tier, loading, _DELAMINATED)).run()
    Ex, Ey, Gxy, nuxy = _LAM.effective_engineering_constants()
    assert (r.pristine_E11_MPa, r.pristine_E22_MPa, r.pristine_G12_MPa) == (Ex, Ey, Gxy)
    # Units: MPa. An IM7/8552 quasi-isotropic laminate is ~60 GPa and
    # isotropic in-plane, G = E / (2 (1 + nu)).
    assert 55e3 < r.pristine_E11_MPa < 70e3
    assert r.pristine_E22_MPa == pytest.approx(r.pristine_E11_MPa, rel=1e-9)
    assert r.pristine_G12_MPa == pytest.approx(r.pristine_E11_MPa / (2 * (1 + nuxy)), rel=1e-9)
    # Nothing models a residual E22 or G12.
    assert r.residual_E22_MPa is r.E22_knockdown is r.residual_G12_MPa is r.G12_knockdown is None


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("tier", ["empirical", "semi_analytical"])
def test_closed_form_tiers_leave_the_residual_stiffness_unmodelled(tier, loading):
    """The closed-form tiers model the damage as a strength equivalent (a
    hole, a notch, buckled sublaminates), with no stiffness to homogenise,
    and no residual-stiffness data exists to fit a correlation to."""
    for damage in (_UNDAMAGED, _DELAMINATED):
        r = BvidAnalysis(_cfg(tier, loading, damage)).run()
        assert r.residual_E11_MPa is None and r.E11_knockdown is None
        assert "residual not modelled by this tier" in r.summary()


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_fe3d_undamaged_panel_keeps_its_stiffness(loading):
    r = BvidAnalysis(_cfg("fe3d", loading, _UNDAMAGED)).run()
    assert r.E11_knockdown == 1.0
    assert r.residual_E11_MPa == r.pristine_E11_MPa


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_fe3d_residual_modulus_is_below_pristine(loading):
    r = BvidAnalysis(_cfg("fe3d", loading, _DELAMINATED)).run()
    assert r.E11_knockdown is not None and 0.0 < r.E11_knockdown < 1.0
    assert r.residual_E11_MPa == pytest.approx(r.pristine_E11_MPa * r.E11_knockdown, rel=1e-12)
    assert r.residual_E11_MPa < r.pristine_E11_MPa
    assert any(n.startswith("fe3d E11 knockdown") and loading in n for n in r.notes), r.notes


def test_fe3d_tension_keeps_more_stiffness_than_compression():
    """In tension the damage zone keeps its fibre-direction stiffness (its
    matrix-dominated terms alone are softened); in compression its whole
    in-plane stiffness is."""
    tension = BvidAnalysis(_cfg("fe3d", "tension", _DELAMINATED)).run().E11_knockdown
    compression = BvidAnalysis(_cfg("fe3d", "compression", _DELAMINATED)).run().E11_knockdown
    assert compression < tension < 1.0


def test_fe3d_E11_knockdown_falls_with_impact_energy():
    knockdowns = []
    for energy_J in (5.0, 15.0, 30.0):
        cfg = AnalysisConfig(
            material="IM7/8552",
            layup_deg=_QI8,
            ply_thickness_mm=0.152,
            panel=PanelGeometry(150.0, 100.0),
            impact=ImpactEvent(energy_J=energy_J),
            tier="fe3d",
        )
        knockdowns.append(BvidAnalysis(cfg).run().E11_knockdown)
    assert all(b < a for a, b in zip(knockdowns, knockdowns[1:])), knockdowns
    assert 0.0 < knockdowns[-1] < knockdowns[0] < 1.0


def _secant_by_energy(cfg, mesh, elements, u, strain):
    """2 U / (V strain^2), with the strain energy summed over the elements."""
    energy = 0.5 * sum(
        float(u[dofs] @ (elem.stiffness_matrix() @ u[dofs]))
        for elem, dofs in zip(elements, mesh.element_dof_maps)
    )
    volume = cfg.panel.Lx_mm * cfg.panel.Ly_mm * _LAM.thickness_mm
    return 2.0 * energy / (volume * strain**2)


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_reaction_force_ratio_agrees_with_direct_computation(loading):
    """The reported knockdown is the ratio of the two reaction-force moduli
    (``axial_secant_modulus``). Each equals the strain-energy modulus
    2U / (V eps^2) of the same solve (the reaction force does all the work
    of the prescribed end strain); the undamaged one is the CLT ``Ex``; the
    damaged one lies below the volume average of the damaged elements'
    stiffness under the undamaged field (the Voigt-type upper bound the
    issue proposed)."""
    tension = loading == "tension"
    sign = 1 if tension else -1
    strain = 0.05
    cfg = _cfg("fe3d", loading, _DELAMINATED)
    moduli = {}
    fields = {}
    for key, damage in (("damaged", _DELAMINATED), ("undamaged", _UNDAMAGED)):
        mesh, elements, _ = ft._fe3d_preflight(cfg, damage, _LAM, tension=tension)
        u = ft._reference_displacement(cfg, mesh, elements, sign, strain)
        moduli[key] = ft.axial_secant_modulus(cfg, mesh, elements, _LAM, sign)
        assert _secant_by_energy(cfg, mesh, elements, u, strain) == pytest.approx(
            moduli[key], rel=1e-6
        )
        fields[key] = (mesh, elements, u)
    assert moduli["undamaged"] == pytest.approx(_LAM.effective_engineering_constants()[0], rel=5e-3)

    r = BvidAnalysis(cfg).run()
    assert r.E11_knockdown == pytest.approx(moduli["damaged"] / moduli["undamaged"], rel=1e-9)

    mesh, elements, _ = fields["damaged"]
    u_undamaged = fields["undamaged"][2]
    voigt = _secant_by_energy(cfg, mesh, elements, u_undamaged, strain)
    assert moduli["damaged"] < voigt < moduli["undamaged"]


def test_stiffness_fields_reach_to_dict_and_json():
    r = BvidAnalysis(_cfg("fe3d", "compression", _DELAMINATED)).run()
    d = json.loads(json.dumps(r.to_dict()))
    for key in _STIFFNESS_KEYS:
        assert key in d
        assert d[key] == getattr(r, key)
    assert d["residual_E22_MPa"] is None
    assert "E11_MPa" in r.summary() and "knockdown" in r.summary()


def test_results_without_stiffness_still_construct_and_summarise():
    """The stiffness fields default to ``None``, so an ``AnalysisResults``
    built without them (as before issue #114) still works."""
    r = AnalysisResults(
        residual_strength_MPa=100.0,
        pristine_strength_MPa=200.0,
        knockdown=0.5,
        damage=_UNDAMAGED,
        dpa_mm2=0.0,
        tier_used="empirical",
        config_snapshot={},
    )
    assert all(getattr(r, key) is None for key in _STIFFNESS_KEYS)
    assert "E11" not in r.summary()
    assert all(r.to_dict()[key] is None for key in _STIFFNESS_KEYS)


def test_fe3d_strength_is_unchanged_by_the_stiffness_output():
    """The modulus comes out of the existing solves: the strength wrappers
    still return what the ``*_and_modulus`` functions do."""
    cfg = _cfg("fe3d", "tension", _DELAMINATED)
    sigma, modulus = ft.fe3d_tai_and_modulus(cfg, _DELAMINATED, _LAM, math.inf)
    assert ft.fe3d_tai(cfg, _DELAMINATED, _LAM, math.inf) == sigma
    assert np.isfinite(modulus) and modulus > 0
    cfg = _cfg("fe3d", "compression", _DELAMINATED)
    growth, fpf, modulus = ft.fe3d_cai_channels_and_modulus(cfg, _DELAMINATED, _LAM, math.inf)
    assert ft.fe3d_cai_channels(cfg, _DELAMINATED, _LAM, math.inf) == (growth, fpf)
