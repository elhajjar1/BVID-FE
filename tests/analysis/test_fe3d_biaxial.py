"""Biaxial loading in the fe3d tier (issue #115, PR 3): two basis solves
superposed to the load ratio, experimental and not validated.

The main-suite tests run on small meshes (a 60 x 40 mm panel at 5 mm, or
150 x 100 mm at 10 mm); the 5 mm table of the design probe and a 2.5 mm
h-refinement check are ``slow``.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, MeshParams, PanelGeometry
from bvidfe.analysis import fe_tier
from bvidfe.analysis.semi_analytical import _semi_analytical_tai_result, semi_analytical_cai
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.failure.laminate_strength import _ply_stresses_per_unit_direction, unnotched_strength
from bvidfe.failure.larc05 import larc05_index
from bvidfe.failure.lekhnitskii import biaxial_point_stress_factor
from bvidfe.failure.tsai_wu import tsai_wu_index

_IM7 = MATERIAL_LIBRARY["IM7/8552"]
QI8 = [0, 45, -45, 90, 90, -45, 45, 0]
CP8 = [0, 90, 0, 90, 90, 0, 90, 0]
UNDAMAGED = DamageState(delaminations=[], dent_depth_mm=0.0)
# Two delaminations under a 60 x 40 mm panel's centre, every element of
# their footprint at 5 mm in the soft zone.
DAMAGE = DamageState(
    [
        DelaminationEllipse(2, (30.0, 20.0), 10.0, 6.0, 0.0),
        DelaminationEllipse(5, (30.0, 20.0), 9.0, 5.0, 90.0),
    ]
)
RATIOS = [-1.0, -0.5, 0.5, 1.0]


def _cfg(
    loading: str = "compression",
    beta: float = 0.0,
    layup=QI8,
    damage: DamageState | None = DAMAGE,
    panel=(60.0, 40.0),
    mesh: float = 5.0,
    **kw,
) -> AnalysisConfig:
    if "impact" in kw:
        damage = None
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=list(layup),
        ply_thickness_mm=0.152,
        panel=PanelGeometry(*panel),
        loading=loading,
        tier="fe3d",
        damage=damage,
        mesh=MeshParams(in_plane_size_mm=mesh),
        transverse_stress_ratio=beta,
        **kw,
    )


def _lam(cfg: AnalysisConfig) -> Laminate:
    return Laminate(_IM7, cfg.layup_deg, cfg.ply_thickness_mm)


def _fpf(cfg: AnalysisConfig, damage: DamageState) -> float:
    """fe3d's own |sigma_x| at first-ply failure along the load ray."""
    lam = _lam(cfg)
    if cfg.loading == "compression":
        return fe_tier.fe3d_cai_channels_and_modulus(cfg, damage, lam, math.inf)[1]
    return fe_tier.fe3d_tai_and_modulus(cfg, damage, lam, math.inf)[0]


@pytest.fixture
def solve_count(monkeypatch):
    """Count the linear solves of the fe3d tier."""
    calls = {"n": 0}
    real = fe_tier.solve_linear_static

    def counting(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    monkeypatch.setattr(fe_tier, "solve_linear_static", counting)
    return calls


# --- superposition -----------------------------------------------------------


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", RATIOS)
def test_the_ratio_is_achieved(loading, beta):
    """sigma_y / sigma_x of the superposed state, each the reaction force on
    its face over the gross section, recomputed from the combined
    displacement, equals beta on the damaged mesh."""
    cfg = _cfg(loading, beta)
    lam = _lam(cfg)
    sign = -1 if loading == "compression" else 1
    basis, _ = fe_tier._fe3d_basis(
        cfg, DAMAGE, lam, label="test", tension=loading == "tension", strain_sign=sign
    )
    assert basis.u_y is not None
    u, force_x, sigma_y = fe_tier._ratio_displacement(cfg, basis, lam, beta)
    h = lam.thickness_mm
    sigma_x = fe_tier._edge_force(basis.mesh, basis.elements, u, 0, 0) / (cfg.panel.Ly_mm * h)
    sigma_y_u = fe_tier._edge_force(basis.mesh, basis.elements, u, 1, 1) / (cfg.panel.Lx_mm * h)
    assert sigma_x == pytest.approx(force_x / (cfg.panel.Ly_mm * h), rel=1e-9)
    assert sigma_y_u == pytest.approx(sigma_y, rel=1e-9)
    assert sigma_y_u / sigma_x == pytest.approx(beta, abs=1e-8)
    assert np.sign(sigma_x) == sign
    # basis X's free y face carries nothing (up to the solver tolerance)
    assert abs(basis.force_xy) < 1e-6 * abs(basis.force_x)


@pytest.mark.parametrize("layup", [[0, 0, 0, 0], [90, 90, 90, 90], [0, 90, 90, 0]])
@pytest.mark.parametrize(
    ("loading", "criterion", "sign"), [("compression", "larc05", -1), ("tension", "tsai_wu", 1)]
)
@pytest.mark.parametrize("beta", [-1.0, -0.5, 0.0, 0.5, 1.0])
def test_undamaged_first_ply_failure_is_clt(layup, loading, criterion, sign, beta):
    """An undamaged panel whose plies carry no lateral stress at the free
    edges (unidirectional, or cross-ply with its small Poisson mismatch) is
    in the uniform CLT state in both basis solves, so its first-ply failure
    along the ray is the plane-stress CLT one with the same criterion:
    to 1e-6 for the unidirectional panels, 5e-4 for cross-ply."""
    cfg = _cfg(loading, beta, layup=layup, damage=UNDAMAGED)
    lam = _lam(cfg)
    stresses = _ply_stresses_per_unit_direction(
        lam, [False] * lam.n_plies, np.array([sign, sign * beta, 0.0])
    )
    clt = min(_clt_multiplier(s, criterion) for s in stresses)
    rel = 1e-6 if len(set(layup)) == 1 else 5e-4
    assert _fpf(cfg, UNDAMAGED) == pytest.approx(clt, rel=rel)


def _clt_multiplier(ply_stress: np.ndarray, criterion: str) -> float:
    """Smallest multiplier on the plane-stress ply stress at which the
    criterion's index reaches 1."""
    s1, s2, t12 = ply_stress
    v = np.array([s1, s2, 0.0, 0.0, 0.0, t12])
    if criterion == "larc05":
        i = larc05_index(_IM7, v)
        return 1.0 / math.sqrt(i) if i > 0 else math.inf
    i1, i2 = tsai_wu_index(_IM7, v), tsai_wu_index(_IM7, 2.0 * v)
    b = (i2 - 2.0 * i1) / 2.0
    a = i1 - b
    c = (-a + math.sqrt(a * a + 4.0 * b)) / (2.0 * b) if abs(b) > 1e-14 else 1.0 / a
    return c if c > 0 else math.inf


@pytest.mark.parametrize(
    ("loading", "criterion", "sign", "rel"),
    [("compression", "larc05", -1, 0.05), ("tension", "tsai_wu", 1, 0.10)],
)
def test_undamaged_quasi_isotropic_panel_is_near_clt(loading, criterion, sign, rel):
    """On a quasi-isotropic panel the free-edge stresses of basis X lower
    first-ply failure below CLT: by at most 5% in compression (LaRC05) and
    10% in tension (Tsai-Wu; worst at beta = 1, where Tsai-Wu sees the
    interlaminar stresses), as in the design probe. Never above it."""
    lam = Laminate(_IM7, QI8, 0.152)
    for beta in (-1.0, -0.5, 0.0, 0.5, 1.0):
        stresses = _ply_stresses_per_unit_direction(
            lam, [False] * lam.n_plies, np.array([sign, sign * beta, 0.0])
        )
        clt = min(_clt_multiplier(s, criterion) for s in stresses)
        fpf = _fpf(_cfg(loading, beta, damage=UNDAMAGED), UNDAMAGED)
        assert clt * (1 - rel) < fpf < clt * (1 + 1e-6), beta


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", RATIOS)
def test_undamaged_knockdown_is_one(loading, beta):
    r = BvidAnalysis(_cfg(loading, beta, damage=UNDAMAGED)).run()
    assert r.knockdown == 1.0
    assert r.pristine_strength_MPa == unnotched_strength(_lam(_cfg()), loading, beta)
    # the residual is the pristine envelope, so the factor is its ratio
    assert r.biaxial_interaction_factor == r.residual_strength_MPa / unnotched_strength(
        _lam(_cfg()), loading
    )


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_rotation_symmetry_on_a_square_panel(loading):
    """Equibiaxial load on a square panel: the layup turned by 90 deg has
    the same strength. Undamaged to 0.2%; with a centred circular
    delamination the knockdown to 2% (the two basis solves' edges are not
    each other's rotation: x platens with free y edges, y platens with the
    x faces held)."""
    layup = [0, 90, 90, 0]
    rotated = [90, 0, 0, 90]
    damage = DamageState([DelaminationEllipse(1, (30.0, 30.0), 10.0, 10.0, 0.0)])
    kw = dict(panel=(60.0, 60.0))
    u = _fpf(_cfg(loading, 1.0, layup, UNDAMAGED, **kw), UNDAMAGED)
    u_rot = _fpf(_cfg(loading, 1.0, rotated, UNDAMAGED, **kw), UNDAMAGED)
    assert u_rot == pytest.approx(u, rel=2e-3)
    kd = BvidAnalysis(_cfg(loading, 1.0, layup, damage, **kw)).run().knockdown
    kd_rot = BvidAnalysis(_cfg(loading, 1.0, rotated, damage, **kw)).run().knockdown
    assert kd_rot == pytest.approx(kd, rel=0.02)


# --- beta = 0 and continuity ---------------------------------------------------


# fe3d at beta = 0 before biaxial loading (the code of PR 2 of #115), on
# _cfg()'s 60 x 40 mm panel.
_PINNED_BETA0 = {
    "compression": 239.68139734212997,
    "tension": 662.741008758128,
}


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_zero_ratio_solves_once_per_mesh(loading, solve_count):
    """beta = 0 runs the uniaxial code: one solve per mesh (damaged and
    undamaged), no transverse basis, no cache, and the uniaxial result."""
    r = BvidAnalysis(_cfg(loading, 0.0)).run()
    assert solve_count["n"] == 2
    assert r.biaxial_interaction_factor is None
    assert r.residual_transverse_stress_MPa is None
    assert not any("Biaxial" in n or "biaxial" in n for n in r.notes)
    assert fe_tier._ACTIVE_BASIS_CACHE.get() is None
    assert r.residual_strength_MPa == _PINNED_BETA0[loading]


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_nonzero_ratio_solves_twice_per_mesh(loading, solve_count):
    """beta != 0: the x and y basis solves of each mesh (4); the uniaxial
    residual of biaxial_interaction_factor reuses the x solves."""
    r = BvidAnalysis(_cfg(loading, -0.5)).run()
    assert solve_count["n"] == 4
    assert r.biaxial_interaction_factor is not None
    assert fe_tier._ACTIVE_BASIS_CACHE.get() is None


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", [-1.0, 0.5])
def test_interaction_factor_is_over_the_uniaxial_run(loading, beta):
    """biaxial_interaction_factor divides by exactly the residual of the
    same config at beta = 0, from the same run's x solves; the signed
    sigma_y fields are s beta times the strengths."""
    bi = BvidAnalysis(_cfg(loading, beta)).run()
    uni = BvidAnalysis(_cfg(loading, 0.0)).run()
    assert bi.biaxial_interaction_factor == bi.residual_strength_MPa / uni.residual_strength_MPa
    s = -1.0 if loading == "compression" else 1.0
    assert bi.residual_transverse_stress_MPa == s * beta * bi.residual_strength_MPa
    assert bi.pristine_transverse_stress_MPa == s * beta * bi.pristine_strength_MPa
    assert bi.E11_knockdown == uni.E11_knockdown


@pytest.mark.parametrize("loading", ["compression", "tension"])
@pytest.mark.parametrize("beta", [-1e-9, 1e-9])
def test_continuity_at_zero(loading, beta):
    """The superposed path at beta -> 0 meets the uniaxial one."""
    uni = BvidAnalysis(_cfg(loading, 0.0)).run()
    bi = BvidAnalysis(_cfg(loading, beta)).run()
    assert bi.residual_strength_MPa == pytest.approx(uni.residual_strength_MPa, rel=1e-6)
    assert bi.knockdown == pytest.approx(uni.knockdown, rel=1e-6)
    assert bi.biaxial_interaction_factor == pytest.approx(1.0, rel=1e-6)


# --- TAI notch and growth --------------------------------------------------------


@pytest.mark.parametrize("beta", RATIOS)
def test_tai_notch_is_kh_times_F(beta):
    """The fibre-damage notch takes the empirical tier's Lekhnitskii factor:
    the hole's knockdown is k_h F(beta), and the elements under it see their
    stress times 1 / (k_h F)."""
    cfg0 = _cfg("tension", 0.0)
    cfg = _cfg("tension", beta)
    lam = _lam(cfg)
    R, kd0 = fe_tier.tai_fibre_damage_hole(cfg0, DAMAGE, lam)
    R_b, kd = fe_tier.tai_fibre_damage_hole(cfg, DAMAGE, lam)
    F = biaxial_point_stress_factor(lam, R, _IM7.wn_d0_mm, beta)
    assert R_b == R and F != 1.0
    assert kd == pytest.approx(kd0 * F, rel=1e-12)
    mesh = fe_tier.build_fe_mesh(cfg, DAMAGE)
    notch = fe_tier._tai_fibre_damage_notch(cfg, DAMAGE, lam, mesh)
    notch0 = fe_tier._tai_fibre_damage_notch(cfg0, DAMAGE, lam, mesh)
    assert notch is not None and notch0 is not None
    hole = notch0 > 1.0
    assert hole.any()
    np.testing.assert_allclose(notch[hole], 1.0 / (kd0 * F), rtol=1e-12)
    assert (notch[~hole] == 1.0).all()


def _impact_kw(layup, loading, beta):
    return dict(
        material="IM7/8552",
        layup_deg=layup,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading=loading,
        impact=ImpactEvent(energy_J=10.0),
        transverse_stress_ratio=beta,
    )


def test_tension_growth_with_lateral_compression_is_semi_analyticals():
    """Axial tension with lateral compression (beta < 0): a laterally
    compressed film grows its delamination, at semi_analytical's growth
    stress, and here it governs, so fe3d and semi_analytical agree."""
    kw = _impact_kw(CP8, "tension", -1.0)
    f = BvidAnalysis(
        AnalysisConfig(tier="fe3d", mesh=MeshParams(in_plane_size_mm=10.0), **kw)
    ).run()
    s = BvidAnalysis(AnalysisConfig(tier="semi_analytical", **kw)).run()
    lam = Laminate(_IM7, CP8, 0.152)
    growth = _semi_analytical_tai_result(
        lam,
        f.damage,
        f.pristine_strength_MPa,
        100.0,
        boundary="simply_supported",
        transverse_stress_ratio=-1.0,
    ).growth
    assert growth is not None and growth.film_stress_axis == "y"
    assert f.residual_strength_MPa == growth.growth_stress_MPa == s.residual_strength_MPa
    assert f.critical_sublaminate == s.critical_sublaminate == growth.interface_index
    assert any(n.startswith("fe3d TAI: delamination growth governs") for n in f.notes)


@pytest.mark.parametrize("beta", [0.5, 1.0])
def test_no_growth_in_biaxial_tension(beta):
    """Under biaxial tension every film is in tension: first-ply failure
    alone, as in the uniaxial TAI."""
    r = BvidAnalysis(_cfg("tension", beta)).run()
    assert r.critical_sublaminate is None
    assert any("no delamination grows" in n for n in r.notes)


@pytest.mark.parametrize("beta", RATIOS)
def test_compression_growth_is_semi_analyticals(beta):
    """The growth channel in compression is semi_analytical's biaxial
    growth stress, the same interface."""
    cfg = _cfg("compression", beta)
    lam = _lam(cfg)
    growth, _, _ = fe_tier.fe3d_cai_channels_and_modulus(cfg, DAMAGE, lam, math.inf)
    sa = semi_analytical_cai(lam, DAMAGE, 1.0, transverse_stress_ratio=beta)
    assert sa.growth is not None and growth is not None
    assert growth == sa.growth
    r = BvidAnalysis(cfg).run()
    assert r.critical_sublaminate == growth.interface_index
    assert r.residual_strength_MPa <= growth.growth_stress_MPa


# --- notes -------------------------------------------------------------------


@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_biaxial_notes(loading):
    r = BvidAnalysis(_cfg(loading, -0.5)).run()
    first = r.notes[0]
    assert first.startswith("Biaxial sigma_y/sigma_x = -0.50: experimental, not validated")
    assert "fe3d superposes two linear solves of each mesh" in first
    assert "softened for the sense of sigma_x alone" in first
    assert "E11_knockdown is the uniaxial (x solve) one" in first
    assert f"biaxial_interaction_factor {r.biaxial_interaction_factor:.3f}" in first
    if loading == "tension":
        assert "1 / (knockdown x F)" in first
        assert any(n.startswith("TAI: equivalent hole") and "F = " in n for n in r.notes)
    else:
        assert any("along the load ray" in n and "whole-panel" in n for n in r.notes)
    assert r.tier_used == "fe3d"


def test_buckling_notes_use_the_biaxial_kernel():
    """Buckling onset (notes and buckling_eigenvalues only) is along the
    load ray: lateral tension raises the whole-panel buckling load."""
    from bvidfe.analysis.semi_analytical import panel_buckling_load

    lam = _lam(_cfg())
    for beta in (-0.5, 0.5):
        cfg = _cfg("compression", beta)
        _, _, notes = fe_tier.fe3d_cai_buckling(cfg, DAMAGE, lam, math.inf)
        N = panel_buckling_load(
            lam, 60.0, 40.0, loading="compression", transverse_stress_ratio=beta
        )
        assert f"fe3d whole-panel buckling at {N / lam.thickness_mm:.1f} MPa" in notes[0]
    assert panel_buckling_load(lam, 60.0, 40.0, transverse_stress_ratio=-0.5) > panel_buckling_load(
        lam, 60.0, 40.0
    )


# --- the basis cache ---------------------------------------------------------


def test_basis_cache_reuses_solves_across_ratios(solve_count):
    """Within fe3d_basis_cache every ratio of the same config and damage
    shares the damaged and undamaged basis solves, and each result equals
    an uncached run."""
    uncached = {b: BvidAnalysis(_cfg("compression", b)).run() for b in (0.0, -1.0, 0.5)}
    solve_count["n"] = 0
    with fe_tier.fe3d_basis_cache() as cache:
        cached = {b: BvidAnalysis(_cfg("compression", b)).run() for b in (0.0, -1.0, 0.5)}
        assert len(cache) == 2
    assert solve_count["n"] == 4
    for b in cached:
        assert cached[b].to_dict() == uncached[b].to_dict()


def test_basis_cache_evicts_the_oldest():
    cache = fe_tier.Fe3dBasisCache(max_entries=1)
    cache.put(("a",), "A")  # type: ignore[arg-type]
    cache.put(("b",), "B")  # type: ignore[arg-type]
    assert len(cache) == 1 and cache.get(("a",)) is None and cache.get(("b",)) == "B"
    with pytest.raises(ValueError, match="max_entries"):
        fe_tier.Fe3dBasisCache(max_entries=0)


def test_basis_cache_keys_on_loading_and_damage(solve_count):
    """A different damage state or loading sense is not a cache hit."""
    with fe_tier.fe3d_basis_cache(max_entries=8):
        BvidAnalysis(_cfg("compression", 0.5)).run()
        BvidAnalysis(_cfg("tension", 0.5)).run()
        other = DamageState([replace(DAMAGE.delaminations[0], major_mm=12.0)])
        BvidAnalysis(_cfg("compression", 0.5, damage=other)).run()
    # damaged and undamaged x and y solves for compression and tension, and
    # the other damage's damaged x and y (its undamaged mesh is a hit)
    assert solve_count["n"] == 4 + 4 + 2


# --- slow: the design probe's 5 mm table and an h-refinement check ----------------


# fe3d knockdowns of the 8-ply QI IM7/8552 150 x 100 mm panel at 10 J,
# 5 mm mesh (the design probe's panel), along the ratio.
_PINNED_5MM = {
    "compression": [
        0.3309263769793906,
        0.3159781827615896,
        0.2772585806120491,
        0.25785151961317104,
        0.1723540980642513,
    ],
    "tension": [
        0.7682895796409515,
        0.7620540533539746,
        0.7627693684743642,
        0.7796753623056926,
        0.8308226976895805,
    ],
}
# The same panel at 2.5 mm, beta = -1 and +1.
_PINNED_2_5MM = {
    "compression": [0.35450126234639356, 0.19348909855204677],
    "tension": [0.7663953874665895, 0.919461319882038],
}


@pytest.mark.slow
@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_pinned_knockdowns_at_5mm(loading, solve_count):
    """The design probe's table (beta = 0 is today's uniaxial fe3d), from
    one sweep: four solves for the five ratios."""
    from bvidfe.sweep.parametric_sweep import sweep_transverse_ratios

    cfg = AnalysisConfig(
        tier="fe3d", mesh=MeshParams(in_plane_size_mm=5.0), **_impact_kw(QI8, loading, 0.0)
    )
    df = sweep_transverse_ratios(cfg, [-1.0, -0.5, 0.0, 0.5, 1.0])
    assert solve_count["n"] == 4
    assert df["knockdown"].tolist() == pytest.approx(_PINNED_5MM[loading], rel=1e-6)


@pytest.mark.slow
@pytest.mark.parametrize("loading", ["compression", "tension"])
def test_h_refinement_at_2_5mm(loading):
    """Halving the mesh moves the knockdown at beta = -1 and +1 by under
    13% (compression +7% / +12%, tension -0.2% / +11%; the uniaxial one
    moves +3.7% / +0.7%): the 10 J damage zone is only a few elements
    across at 5 mm (see FPF_AVERAGING_RADIUS_MM)."""
    from bvidfe.sweep.parametric_sweep import sweep_transverse_ratios

    cfg = AnalysisConfig(
        tier="fe3d", mesh=MeshParams(in_plane_size_mm=2.5), **_impact_kw(QI8, loading, 0.0)
    )
    kd = sweep_transverse_ratios(cfg, [-1.0, 1.0])["knockdown"].tolist()
    assert kd == pytest.approx(_PINNED_2_5MM[loading], rel=1e-6)
    kd5 = [_PINNED_5MM[loading][0], _PINNED_5MM[loading][-1]]
    for fine, coarse in zip(kd, kd5):
        assert abs(fine / coarse - 1) < 0.13
