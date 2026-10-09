"""Unnotched laminate strength: ply-discount last-ply failure.

The pristine strength used to be a thickness-weighted average of the ply
strengths along the load, sum t_i (X cos^2 + Y sin^2) / sum t_i. Plies share
strain, not stress, so that credits the off-axis plies with strength they
never reach before the 0 deg fibres break: it overpredicted the NCAMP
unnotched strengths by 17-185% (see tests/validation/test_unnotched_strength_ncamp.py).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from bvidfe import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import PanelGeometry
from bvidfe.core.laminate import Laminate, _reduced_stiffness
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState
from bvidfe.failure.laminate_strength import PLY_DISCOUNT_FACTOR, unnotched_strength

_M = MATERIAL_LIBRARY["AS4/8552"]
_QI = [45, 0, -45, 90, 90, -45, 0, 45]


@pytest.mark.parametrize("name", sorted(MATERIAL_LIBRARY))
def test_unidirectional_0_fails_at_the_fibre_strengths(name):
    m = MATERIAL_LIBRARY[name]
    lam = Laminate(m, [0, 0, 0, 0], 0.2)
    assert unnotched_strength(lam, "tension") == pytest.approx(m.Xt, rel=1e-9)
    assert unnotched_strength(lam, "compression") == pytest.approx(m.Xc, rel=1e-9)


@pytest.mark.parametrize("name", sorted(MATERIAL_LIBRARY))
def test_unidirectional_90_fails_when_its_matrix_fails(name):
    """Every ply cracked and nothing left to carry load: the laminate fails."""
    m = MATERIAL_LIBRARY[name]
    lam = Laminate(m, [90, 90, 90, 90], 0.2)
    assert unnotched_strength(lam, "tension") == pytest.approx(m.Yt, rel=1e-9)
    assert unnotched_strength(lam, "compression") == pytest.approx(m.Yc, rel=1e-9)


def test_cross_ply_tension_is_carried_by_the_0_plies_after_the_90s_crack():
    lam = Laminate(_M, [0, 90, 90, 0], 0.2)
    assert unnotched_strength(lam, "tension") == pytest.approx(_M.Xt / 2, rel=0.01)


def test_ply_thickness_weights_the_load_share():
    """0 deg ply twice as thick as the 90: it carries 2/3 of the section."""
    lam = Laminate(_M, [0, 90], [0.20, 0.10])
    assert unnotched_strength(lam, "tension") == pytest.approx(_M.Xt * 2 / 3, rel=0.01)


def test_quasi_isotropic_compression_is_0_deg_fibre_failure():
    """No ply cracks first, so the laminate fails when the 0 deg ply reaches
    Xc: sigma_1 = eps_x (Q11 - Q12 nu_xy) and sigma_x = Ex eps_x."""
    lam = Laminate(_M, _QI, 0.188)
    Ex, _, _, nu_xy = lam.effective_engineering_constants()
    Q = _reduced_stiffness(_M)
    expected = _M.Xc * Ex / (Q[0, 0] - Q[0, 1] * nu_xy)
    assert unnotched_strength(lam, "compression") == pytest.approx(expected, rel=1e-9)


def test_quasi_isotropic_tension_lies_between_first_ply_and_fibre_failure():
    """The 90s and +-45s crack before the 0s break: the cracked plies shed
    load onto the 0s, so the laminate fails below the intact-laminate fibre
    failure stress but well above first-ply (matrix) failure."""
    lam = Laminate(_M, _QI, 0.188)
    Ex, _, _, nu_xy = lam.effective_engineering_constants()
    Q = _reduced_stiffness(_M)
    intact_fibre = _M.Xt * Ex / (Q[0, 0] - Q[0, 1] * nu_xy)
    sigma = unnotched_strength(lam, "tension")
    assert 0.8 * intact_fibre < sigma < intact_fibre


def test_old_ply_average_overpredicts_the_quasi_isotropic_laminate():
    lam = Laminate(_M, _QI, 0.188)
    n = len(_QI)
    ply_average = (
        sum(
            _M.Xc * math.cos(math.radians(t)) ** 2 + _M.Yc * math.sin(math.radians(t)) ** 2
            for t in _QI
        )
        / n
    )
    assert unnotched_strength(lam, "compression") < 0.7 * ply_average


def test_strength_is_insensitive_to_the_ply_discount_factor(monkeypatch):
    import bvidfe.failure.laminate_strength as ls

    lam = Laminate(_M, _QI, 0.188)
    base = unnotched_strength(lam, "tension")
    monkeypatch.setattr(ls, "PLY_DISCOUNT_FACTOR", PLY_DISCOUNT_FACTOR / 10)
    assert unnotched_strength(lam, "tension") == pytest.approx(base, rel=0.02)


def test_strength_is_positive_and_finite_for_every_preset_and_loading():
    for m in MATERIAL_LIBRARY.values():
        for layup in ([45, -45, -45, 45], _QI, [0, 45, 90, -45, -45, 90, 45, 0]):
            lam = Laminate(m, layup, 0.15)
            for loading in ("tension", "compression"):
                s = unnotched_strength(lam, loading)
                assert np.isfinite(s) and s > 0, (m.name, layup, loading)


def test_bvid_analysis_uses_the_unnotched_strength_as_pristine():
    cfg = AnalysisConfig(
        material="AS4/8552",
        layup_deg=_QI,
        ply_thickness_mm=0.188,
        panel=PanelGeometry(150.0, 100.0),
        loading="compression",
        tier="empirical",
        damage=DamageState([], dent_depth_mm=0.0),
    )
    lam = Laminate(_M, _QI, 0.188)
    r = BvidAnalysis(cfg).run()
    assert r.pristine_strength_MPa == pytest.approx(unnotched_strength(lam, "compression"))
