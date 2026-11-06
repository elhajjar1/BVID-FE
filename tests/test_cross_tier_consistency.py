"""Cross-tier knockdown consistency checks.

Issue #11: existing tests cover each tier in isolation but never assert the
documented cross-tier relationships from docs/physics_models.md "Knockdown definition and
cross-tier comparability":

  * All three tiers share the same pristine baseline (the unnotched
    laminate strength), so ``pristine_strength_MPa`` must be identical.
  * For TAI, ``empirical`` and ``semi_analytical`` delegate to the same
    fibre-damage open-hole model (``fibre_damage_tai``), so their
    knockdowns must match exactly.
  * For CAI, ``semi_analytical`` takes the lower of the empirical Soutis
    term and the delamination growth stress, so its residual is always <=
    the empirical residual on the same input — and therefore
    ``semi_analytical.knockdown <= empirical.knockdown``.
  * On a pristine input (no impact damage, zero DPA) all three tiers must
    return ``knockdown == 1.0`` exactly (empirical / semi_analytical) or
    very close to 1.0 (fe3d, capped at sigma_pristine_MPa).
"""

from __future__ import annotations

import math

import pytest

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.analysis.semi_analytical import weakest_sublaminate_growth
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.damage.state import DamageState
from bvidfe.impact.mapping import ImpactEvent


def _base_kwargs(loading: str = "compression"):
    return dict(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading=loading,
        impact=ImpactEvent(20.0, ImpactorGeometry(), mass_kg=5.5),
    )


def _run(tier: str, loading: str = "compression"):
    return BvidAnalysis(AnalysisConfig(tier=tier, **_base_kwargs(loading))).run()


def test_pristine_strength_identical_across_tiers():
    """All three tiers share _pristine_strength()."""
    e = _run("empirical")
    s = _run("semi_analytical")
    assert e.pristine_strength_MPa == s.pristine_strength_MPa


def test_tai_empirical_equals_semi_analytical():
    """semi_analytical_tai delegates unchanged to fibre_damage_tai."""
    e = _run("empirical", loading="tension")
    s = _run("semi_analytical", loading="tension")
    assert e.knockdown == pytest.approx(s.knockdown, rel=1e-12)


def test_cai_semi_analytical_le_empirical():
    """semi_analytical_cai = min(empirical_soutis, delamination_growth)."""
    e = _run("empirical", loading="compression")
    s = _run("semi_analytical", loading="compression")
    assert s.knockdown <= e.knockdown + 1e-9


def test_pristine_input_yields_unity_knockdown():
    """With explicit DamageState() containing no delaminations and zero
    dent depth, both empirical and semi_analytical must return knockdown
    exactly 1.0 (no early-exit short-circuits, no floating-point drift)."""
    pristine = DamageState(delaminations=[], dent_depth_mm=0.0)
    common = dict(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.2,
        panel=PanelGeometry(100, 50),
        loading="compression",
        damage=pristine,
    )
    for tier in ("empirical", "semi_analytical"):
        cfg = AnalysisConfig(tier=tier, **common)
        r = BvidAnalysis(cfg).run()
        assert r.knockdown == pytest.approx(
            1.0, rel=1e-12
        ), f"tier={tier} knockdown={r.knockdown!r} should be 1.0 on pristine input"


def test_knockdown_is_finite_and_bounded():
    """No tier should ever return NaN, Inf, or knockdown outside [0, 1+ε]."""
    for tier in ("empirical", "semi_analytical"):
        r = _run(tier)
        assert math.isfinite(r.knockdown)
        assert 0.0 <= r.knockdown <= 1.0 + 1e-9


def test_residual_never_exceeds_pristine():
    """Per the per-tier sigma cap (semi_analytical & fe3d explicitly,
    empirical implicitly via Soutis kd <= 1) residual <= pristine."""
    for tier in ("empirical", "semi_analytical"):
        r = _run(tier)
        assert r.residual_strength_MPa <= r.pristine_strength_MPa + 1e-6


# --- biaxial load (issue #115) ---------------------------------------------


def _run_biaxial(tier: str, loading: str, beta: float, layup=None):
    kw = _base_kwargs(loading)
    if layup is not None:
        kw["layup_deg"] = layup
    return BvidAnalysis(AnalysisConfig(tier=tier, transverse_stress_ratio=beta, **kw)).run()


@pytest.mark.parametrize("layup", [None, [0, 90, 0, 90, 90, 0, 90, 0]])
@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", [-1.0, -0.5, 0.5, 1.0])
def test_biaxial_semi_analytical_le_empirical(layup, loading, beta):
    """Along every load ray the two closed-form tiers share the pristine
    strength, and semi_analytical takes the lower of the empirical tier's
    result (homothetic Soutis in compression, the hole model with F(beta) in
    tension) and the biaxial growth stress, so it is never above it. Under
    biaxial tension (beta > 0) no film is compressed and the two are
    equal."""
    e = _run_biaxial("empirical", loading, beta, layup)
    s = _run_biaxial("semi_analytical", loading, beta, layup)
    assert e.pristine_strength_MPa == s.pristine_strength_MPa
    assert s.knockdown <= e.knockdown + 1e-12
    if loading == "tension" and beta > 0:
        assert s.knockdown == e.knockdown
    assert 0.0 < s.knockdown <= 1.0


@pytest.mark.parametrize(
    ("loading", "beta"), [("compression", -1.0), ("tension", -1.0), ("tension", 0.5)]
)
def test_biaxial_fe3d_shares_pristine_and_growth_with_semi_analytical(loading, beta):
    """fe3d at a nonzero ratio (10 mm mesh, to stay fast): the pristine
    envelope point is every tier's, and its growth channel is
    semi_analytical's biaxial growth stress, on the same interface, so its
    residual is never above that growth stress."""
    from bvidfe.analysis import MeshParams

    kw = _base_kwargs(loading)
    s = BvidAnalysis(
        AnalysisConfig(tier="semi_analytical", transverse_stress_ratio=beta, **kw)
    ).run()
    f = BvidAnalysis(
        AnalysisConfig(
            tier="fe3d",
            transverse_stress_ratio=beta,
            mesh=MeshParams(in_plane_size_mm=10.0),
            **kw,
        )
    ).run()
    assert f.pristine_strength_MPa == s.pristine_strength_MPa
    assert f.critical_sublaminate == s.critical_sublaminate
    growth = weakest_sublaminate_growth(
        Laminate(MATERIAL_LIBRARY["IM7/8552"], kw["layup_deg"], 0.152),
        f.damage,
        loading=loading,
        transverse_stress_ratio=beta,
    )
    assert growth is not None
    if loading == "tension" and beta > 0:
        assert f.critical_sublaminate is None and math.isinf(growth.growth_stress_MPa)
    else:
        assert f.critical_sublaminate == growth.interface_index
        assert f.residual_strength_MPa <= growth.growth_stress_MPa
        assert f"grows at {growth.growth_stress_MPa:.1f} MPa" in " ".join(f.notes)
    assert 0.0 < f.knockdown <= 1.0
