"""Numerical equivalence between the vectorised and scalar-reference
implementations of ``_solve_failure_strain_analytic``.

The vectorised production routine in ``bvidfe.analysis.fe_tier`` replaces
a per-Gauss-point Python loop with batched ``larc05_index_batch`` /
``tsai_wu_index_batch`` calls. ``_solve_failure_strain_analytic_scalar_ref``
is the unchanged pre-vectorisation implementation, kept solely so this
test can prove the two paths return identical strain-at-failure values
on representative inputs. Both check every Gauss point on its own
(``averaging_radius_mm=None``); the production default averages the index
over a disk, which the reference does not. Any future regression in either path will
trip this test.

Three cases are exercised:

  * LaRC05 path  — the pure-quadratic branch: ``c = 1/sqrt(idx_ref)``.
  * Puck path    — the degree-1 homogeneous branch: ``c = 1/idx_ref``.
  * Tsai-Wu path — the affine-quadratic branch including the linear-
    fallback (``|b| < 1e-14``) and disc < 0 mask edges.

Both run on a small fe3d mesh (~100 elements) so the test stays under a
second.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from bvidfe.analysis import AnalysisConfig, MeshParams
from bvidfe.analysis.fe_mesh import build_fe_mesh
from bvidfe.analysis.fe_tier import (
    _build_elements,
    _failure_stress_material,
    _solve_failure_strain_analytic,
    _solve_failure_strain_analytic_scalar_ref,
)
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.puck import puck_index_batch
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.solver.boundary import uniaxial_x_bcs
from bvidfe.solver.static import solve_linear_static


def _build_setup(damage: DamageState, layup_deg=(0.0, 90.0, 0.0, 90.0)):
    cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=list(layup_deg),
        ply_thickness_mm=0.2,
        panel=PanelGeometry(20.0, 10.0),
        loading="compression",
        tier="fe3d",
        impact=ImpactEvent(5.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=MeshParams(elements_per_ply=1, in_plane_size_mm=2.5),
    )
    lam = Laminate(
        material=MATERIAL_LIBRARY[cfg.material],
        layup_deg=cfg.layup_deg,
        ply_thickness_mm=cfg.ply_thickness_mm,
    )
    mesh = build_fe_mesh(cfg, damage)
    elements = _build_elements(mesh, lam)
    return cfg, mesh, elements


@pytest.mark.parametrize("criterion", ["larc05", "tsai_wu", "puck"])
@pytest.mark.parametrize(
    "damage",
    [
        DamageState(delaminations=[], dent_depth_mm=0.0),
        DamageState(
            delaminations=[DelaminationEllipse(1, (10, 5), 6, 3, 0)],
            dent_depth_mm=0.2,
        ),
    ],
    ids=["pristine", "delaminated"],
)
def test_vectorised_matches_scalar_reference(criterion, damage):
    """Compression FPF: vectorised result must equal scalar-ref result."""
    cfg, mesh, elements = _build_setup(damage)
    eps_vec = _solve_failure_strain_analytic(
        cfg, mesh, elements, strain_sign=-1, criterion=criterion, averaging_radius_mm=None
    )
    eps_ref = _solve_failure_strain_analytic_scalar_ref(
        cfg, mesh, elements, strain_sign=-1, criterion=criterion
    )
    assert np.isfinite(eps_vec)
    assert np.isfinite(eps_ref)
    # Both paths share the same FE solve and the same algebra; the only
    # numerical-evaluation difference is BLAS-routed batched dot products
    # vs Python multiplications. Drift should be at floating-point noise
    # level.
    assert eps_vec == pytest.approx(eps_ref, rel=1e-10, abs=1e-12)


def test_vectorised_tension_path_also_matches_scalar_reference():
    """The tension path uses uniaxial_x_bcs (positive strain) and Tsai-Wu by
    default; verify the vectorised form is still equivalent under the
    +strain_sign branch."""
    damage = DamageState(
        delaminations=[DelaminationEllipse(1, (10, 5), 6, 3, 0)],
        dent_depth_mm=0.0,
    )
    cfg, mesh, elements = _build_setup(damage)
    eps_vec = _solve_failure_strain_analytic(
        cfg, mesh, elements, strain_sign=+1, criterion="tsai_wu", averaging_radius_mm=None
    )
    eps_ref = _solve_failure_strain_analytic_scalar_ref(
        cfg, mesh, elements, strain_sign=+1, criterion="tsai_wu"
    )
    assert np.isfinite(eps_vec) and np.isfinite(eps_ref)
    assert eps_vec == pytest.approx(eps_ref, rel=1e-10, abs=1e-12)


def test_fpf_strain_does_not_depend_on_panel_boundary():
    """The membrane solve holds the mid-surface at u_z = 0 whatever
    ``panel.boundary`` is (see ``uniaxial_x_bcs``): out-of-plane edge
    supports carry no load in a linear in-plane solve of a flat panel, and
    ``boundary`` acts through the closed-form buckling and delamination
    growth instead.

    Issue #32 had the builder ignore ``boundary`` silently; it then held
    ``u_z = 0`` on whole edge faces, which gave clamped < simply_supported
    < free here (0.952 / 0.956 / 0.958%) only because more restrained faces
    blocked more of the plies' through-thickness Poisson expansion. Even
    this unsymmetric layup, which would bend if it could, gives one value.
    """
    damage = DamageState(delaminations=[], dent_depth_mm=0.0)
    cfg, mesh, elements = _build_setup(damage)

    def eps_for(boundary: str) -> float:
        panel = dataclasses.replace(cfg.panel, boundary=boundary)
        cfg_b = dataclasses.replace(cfg, panel=panel)
        return _solve_failure_strain_analytic(
            cfg_b, mesh, elements, strain_sign=-1, criterion="larc05"
        )

    eps_ss = eps_for("simply_supported")
    assert eps_for("clamped") == eps_ss
    assert eps_for("free") == eps_ss


# ---------------------------------------------------------------------------
# Criterion dispatch. Before v0.2.1 any criterion other than "larc05" fell
# through to the Tsai-Wu branch, so "puck" (and typos) silently returned
# Tsai-Wu results.
# ---------------------------------------------------------------------------

# Angle-ply layup where inter-fiber failure contributes, so the three
# criteria give distinct failure strains (on cross-ply, fiber compression
# governs and Puck/LaRC05 coincide at |s1|/Xc).
_ANGLE_PLY = (45.0, -45.0, -45.0, 45.0)
_DELAM = DamageState(
    delaminations=[DelaminationEllipse(1, (10, 5), 6, 3, 0)],
    dent_depth_mm=0.2,
)


@pytest.mark.parametrize(
    "solve", [_solve_failure_strain_analytic, _solve_failure_strain_analytic_scalar_ref]
)
def test_unknown_criterion_raises(solve):
    cfg, mesh, elements = _build_setup(_DELAM)
    with pytest.raises(ValueError, match="unknown failure criterion 'tsaiwu'"):
        solve(cfg, mesh, elements, strain_sign=-1, criterion="tsaiwu")


@pytest.mark.parametrize("strain_sign", [-1, +1])
def test_puck_is_not_evaluated_as_another_criterion(strain_sign):
    cfg, mesh, elements = _build_setup(_DELAM, layup_deg=_ANGLE_PLY)
    eps = {
        c: _solve_failure_strain_analytic(cfg, mesh, elements, strain_sign=strain_sign, criterion=c)
        for c in ("tsai_wu", "larc05", "puck")
    }
    assert eps["puck"] < 0.05, "expected failure below the strain cap"
    assert eps["puck"] != pytest.approx(eps["tsai_wu"], rel=1e-6)
    assert eps["puck"] != pytest.approx(eps["larc05"], rel=1e-6)


@pytest.mark.parametrize("strain_sign", [-1, +1])
def test_puck_index_reaches_one_at_returned_strain(strain_sign):
    """Re-solving at the returned strain puts the max Puck index, on the
    stress the criterion sees (effective in the damage zone), at exactly 1."""
    cfg, mesh, elements = _build_setup(_DELAM, layup_deg=_ANGLE_PLY)
    eps = _solve_failure_strain_analytic(
        cfg, mesh, elements, strain_sign=strain_sign, criterion="puck", averaging_radius_mm=None
    )
    bcs = uniaxial_x_bcs(mesh.node_coords, strain_sign * eps, boundary=cfg.panel.boundary)
    u = solve_linear_static(elements, mesh.element_dof_maps, mesh.n_dof, bcs)
    material = MATERIAL_LIBRARY[cfg.material]
    max_idx = max(
        float(puck_index_batch(material, _failure_stress_material(elem, u[dofs], f_ip)).max())
        for elem, dofs, f_ip in zip(elements, mesh.element_dof_maps, mesh.in_plane_damage_factors)
    )
    assert max_idx == pytest.approx(1.0, rel=1e-8)
