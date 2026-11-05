"""``transverse_stress_ratio = 0`` changes no result (issue #115, PR 1).

Every numeric ``AnalysisResults`` field, ``notes`` and ``warnings`` are
equal (``==``) between a config without the field and one with ``-0.0``
(which ``AnalysisConfig`` stores as ``0.0``, so this also covers ``0.0``), for the example configs, the synthetic dataset and the first case
of each validation dataset, in all three tiers (fe3d on a 10 mm mesh). The
uniaxial path is also checked not to touch the biaxial code: the
Lekhnitskii factor is never evaluated, and the empirical tier's strengths
are the pre-change call chain on the pre-change pristine strength
(``_reference_uniaxial``).
"""

from __future__ import annotations

import importlib.util
import json
import sys
import warnings
from dataclasses import replace
from pathlib import Path

import pytest

from bvidfe import AnalysisConfig, BvidAnalysis, ImpactEvent, MeshParams, PanelGeometry
from bvidfe.analysis import bvid
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.io import load_cscan_json
from bvidfe.failure import soutis_openhole
from bvidfe.failure.soutis_openhole import (
    SOUTIS_REFERENCE_AREA_MM2,
    lekhnitskii_kt_infinity,
    soutis_cai,
)
from tests.failure.test_laminate_strength_biaxial import _reference_uniaxial

REPO = Path(__file__).resolve().parent.parent
QI8 = [0, 45, -45, 90, 90, -45, 45, 0]


def _load_validator():
    spec = importlib.util.spec_from_file_location(
        "validate_bvid_public", REPO / "validation" / "validate_bvid_public.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["validate_bvid_public"] = mod
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop("validate_bvid_public", None)
    return mod


def _dataset_configs() -> dict:
    """The config ``run_case`` builds for the first case of every dataset."""
    v = _load_validator()
    out = {}
    for path in sorted(v.DATASET_DIR.glob("*.json")):
        _, cases, _ = v.load_dataset(path)
        seen = {}

        class Spy:
            def __init__(self, config):
                seen["config"] = config
                raise _Stop

        v.BvidAnalysis = Spy
        try:
            v.run_case(cases[0], "empirical")
        except _Stop:
            pass
        out[path.stem] = seen["config"]
    return out


class _Stop(Exception):
    pass


def _example_configs() -> dict:
    base = dict(
        material="IM7/8552",
        layup_deg=list(QI8),
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150.0, 100.0),
    )
    return {
        "quickstart_cai_10J": AnalysisConfig(**base, impact=ImpactEvent(energy_J=10.0)),
        "tai_25J": AnalysisConfig(**base, loading="tension", impact=ImpactEvent(energy_J=25.0)),
        "sample_cscan": AnalysisConfig(
            **base, damage=load_cscan_json(REPO / "examples" / "sample_cscan.json")
        ),
        "coupon_on_window": AnalysisConfig(
            **{**base, "panel": PanelGeometry(80.0, 50.0)},
            impact=ImpactEvent(energy_J=8.0),
            impact_panel=PanelGeometry(60.0, 60.0, boundary="clamped"),
        ),
    }


_CONFIGS = {**_example_configs(), **_dataset_configs()}
_TIERS = ["empirical", "semi_analytical", "fe3d"]
_LOADINGS = ["compression", "tension"]


def _variant(cfg: AnalysisConfig, tier: str, loading: str, ratio) -> AnalysisConfig:
    kw = dict(tier=tier, loading=loading)
    if tier == "fe3d":
        kw["mesh"] = MeshParams(in_plane_size_mm=10.0)
    out = replace(cfg, **kw)
    if ratio is not None:
        out = replace(out, transverse_stress_ratio=ratio)
    return out


def _comparable(result) -> str:
    d = result.to_dict()
    return json.dumps(d, sort_keys=True, default=repr)


@pytest.fixture
def no_biaxial_code(monkeypatch):
    """The uniaxial path must not evaluate the Lekhnitskii factor."""

    def boom(*a, **k):
        raise AssertionError("biaxial_point_stress_factor called at ratio 0")

    monkeypatch.setattr(bvid, "biaxial_point_stress_factor", boom)
    monkeypatch.setattr(soutis_openhole, "biaxial_point_stress_factor", boom)


def _run(cfg):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return BvidAnalysis(cfg).run()


@pytest.mark.parametrize("name", sorted(_CONFIGS))
@pytest.mark.parametrize("tier", _TIERS)
@pytest.mark.parametrize("loading", _LOADINGS)
def test_zero_ratio_is_bit_identical_to_no_ratio(name, tier, loading, no_biaxial_code):
    cfg = _CONFIGS[name]
    if tier == "fe3d" and len(cfg.layup_deg) > 16:
        pytest.skip("thick layup: covered by the closed-form tiers")
    results = [_run(_variant(cfg, tier, loading, r)) for r in (None, -0.0)]
    first = _comparable(results[0])
    for r in results[1:]:
        assert _comparable(r) == first
        assert r.notes == results[0].notes and r.warnings == results[0].warnings
    base = results[0]
    assert base.config_snapshot["transverse_stress_ratio"] == 0.0
    for field in (
        "pristine_transverse_stress_MPa",
        "residual_transverse_stress_MPa",
        "biaxial_interaction_factor",
    ):
        assert getattr(base, field) is None


@pytest.mark.parametrize("name", sorted(_CONFIGS))
@pytest.mark.parametrize("loading", _LOADINGS)
def test_empirical_strengths_are_the_pre_change_call_chain(name, loading, no_biaxial_code):
    cfg = _variant(_CONFIGS[name], "empirical", loading, 0.0)
    r = _run(cfg)
    mat = MATERIAL_LIBRARY[cfg.material] if isinstance(cfg.material, str) else cfg.material
    lam = Laminate(mat, cfg.layup_deg, cfg.ply_thickness_mm)
    sigma_0 = _reference_uniaxial(lam, loading)
    if loading == "compression":
        expected = soutis_cai(mat, r.dpa_mm2, SOUTIS_REFERENCE_AREA_MM2, sigma_0)
    else:
        expected = soutis_openhole.fibre_damage_tai(
            mat, r.damage, sigma_0, lekhnitskii_kt_infinity(lam), cfg.panel.Ly_mm
        )
    assert r.pristine_strength_MPa == sigma_0
    assert r.residual_strength_MPa == expected
    assert r.knockdown == expected / sigma_0
