#!/usr/bin/env python3
"""Visual verification of GUI tab outputs across a matrix of inputs.

Builds each GUI tab (Summary, Damage Map, Knockdown Curve, Damage Severity)
off-screen and saves a rendered PNG for every config in a small matrix.
The intent is purely visual QA — if two configurations that should look
different actually produce identical images, that is a silent rendering
regression and this script surfaces it.

Output: ``/tmp/bvidfe_visual/<group>/<case>_<tab>.png``. Groups correspond
to the variable being swept (boundary / shape / energy / diameter / mass).

Usage:
    QT_QPA_PLATFORM=offscreen ./.venv/bin/python scripts/visual_matrix.py
"""
from __future__ import annotations

import hashlib
import os
import warnings
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

warnings.simplefilter("ignore")

from PyQt6.QtWidgets import QApplication

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.gui.tabs.damage_map_tab import DamageMapTab
from bvidfe.gui.tabs.knockdown_tab import KnockdownTab
from bvidfe.gui.tabs.summary_tab import SummaryTab
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.sweep.parametric_sweep import sweep_energies


OUT = Path("/tmp/bvidfe_visual")
OUT.mkdir(parents=True, exist_ok=True)

# Default config used as the baseline — each sweep varies one parameter
# around this.
MAT = "IM7/8552"
LAYUP = [0, 45, -45, 90, 90, -45, 45, 0]
PLY_T = 0.152
# 300x200 panel so DPA doesn't saturate at the 80% cap and we actually
# see input-sensitivity in the downstream tabs.
DEFAULT_PANEL = PanelGeometry(300, 200)
DEFAULT_IMPACT = ImpactEvent(10.0, ImpactorGeometry(), mass_kg=5.5)
DEFAULT_TIER = "empirical"


def _base_kw(**over):
    kw = dict(
        material=MAT,
        layup_deg=LAYUP,
        ply_thickness_mm=PLY_T,
        panel=DEFAULT_PANEL,
        loading="compression",
        tier=DEFAULT_TIER,
        impact=DEFAULT_IMPACT,
    )
    kw.update(over)
    return kw


def _image_hash(widget) -> str:
    """Grab the rendered bytes of a QWidget, hash them. Used to verify two
    configs that *should* differ are actually producing different pixels."""
    pixmap = widget.grab()
    img = pixmap.toImage()
    n = img.width() * img.height() * 4  # RGBA
    ptr = img.bits()
    ptr.setsize(n)
    return hashlib.sha256(bytes(ptr)).hexdigest()[:12]


def render_case(group: str, case: str, cfg_kw: dict) -> dict[str, str]:
    """Run one analysis, render Summary + Damage Map + Knockdown Curve,
    save each tab as a PNG under OUT/group/. Returns a dict of image hashes."""
    cfg = AnalysisConfig(**cfg_kw)
    r = BvidAnalysis(cfg).run()

    group_dir = OUT / group
    group_dir.mkdir(exist_ok=True, parents=True)

    # --- Summary ---
    st = SummaryTab()
    st.resize(600, 400)
    st.update(r)
    st.show()
    QApplication.processEvents()
    summary_path = group_dir / f"{case}_summary.png"
    st.grab().save(str(summary_path))

    # --- Damage Map ---
    dm = DamageMapTab()
    dm.resize(600, 500)
    dm.update(r, panel=cfg.panel)
    dm.show()
    QApplication.processEvents()
    dmap_path = group_dir / f"{case}_damage_map.png"
    dm.grab().save(str(dmap_path))

    # --- Knockdown curve (use a small energy sweep) ---
    sw = sweep_energies(cfg, energies_J=[3, 5, 8, 12, 20])
    kt = KnockdownTab()
    kt.resize(700, 500)
    kt.update_series(sw["energy_J"].tolist(), sw["knockdown"].tolist(),
                     tier_label=cfg.tier)
    kt.show()
    QApplication.processEvents()
    kd_path = group_dir / f"{case}_knockdown.png"
    kt.grab().save(str(kd_path))

    return {
        "residual_MPa": f"{r.residual_strength_MPa:.1f}",
        "KD": f"{r.knockdown:.4f}",
        "DPA": f"{r.damage.projected_damage_area_mm2:.0f}",
        "summary_hash": _image_hash(st),
        "damage_map_hash": _image_hash(dm),
        "knockdown_hash": _image_hash(kt),
    }


def main():
    app = QApplication.instance() or QApplication([])

    groups = {
        "boundary": [
            ("simply_supported", _base_kw(panel=PanelGeometry(300, 200, "simply_supported"))),
            ("clamped", _base_kw(panel=PanelGeometry(300, 200, "clamped"))),
            ("free", _base_kw(panel=PanelGeometry(300, 200, "free"))),
        ],
        "shape": [
            ("hemispherical", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(16.0, "hemispherical"), 5.5))),
            ("flat", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(16.0, "flat"), 5.5))),
            ("conical", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(16.0, "conical"), 5.5))),
        ],
        "diameter": [
            ("d08mm", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(8.0), 5.5))),
            ("d16mm", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(16.0), 5.5))),
            ("d40mm", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(40.0), 5.5))),
        ],
        "mass": [
            ("m1kg", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(), 1.0))),
            ("m5p5kg", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(), 5.5))),
            ("m20kg", _base_kw(impact=ImpactEvent(
                10.0, ImpactorGeometry(), 20.0))),
        ],
        "energy": [
            ("E_03J", _base_kw(impact=ImpactEvent(3.0, ImpactorGeometry(), 5.5))),
            ("E_10J", _base_kw(impact=ImpactEvent(10.0, ImpactorGeometry(), 5.5))),
            ("E_30J", _base_kw(impact=ImpactEvent(30.0, ImpactorGeometry(), 5.5))),
        ],
        "loading": [
            ("compression", _base_kw(loading="compression")),
            ("tension", _base_kw(loading="tension")),
        ],
        "tier": [
            ("empirical", _base_kw(tier="empirical")),
            ("semi_analytical", _base_kw(tier="semi_analytical")),
        ],
    }

    all_reports: list[tuple[str, str, dict]] = []
    for group, cases in groups.items():
        print(f"\n=== group: {group} ===")
        hashes = {}
        for case, cfg_kw in cases:
            print(f"  [{case}]", flush=True)
            rep = render_case(group, case, cfg_kw)
            hashes[case] = rep
            all_reports.append((group, case, rep))
            print(f"    residual={rep['residual_MPa']} MPa   "
                  f"KD={rep['KD']}   DPA={rep['DPA']}   "
                  f"dmap_hash={rep['damage_map_hash']}   "
                  f"kd_hash={rep['knockdown_hash']}")
        # Alert if all cases in this group produced identical rendered damage maps
        # (would indicate the tab didn't respond to the input variation).
        unique_dmap = {r["damage_map_hash"] for r in hashes.values()}
        unique_kd = {r["knockdown_hash"] for r in hashes.values()}
        if len(unique_dmap) < len(hashes):
            print(f"  ⚠ {group}: {len(hashes)} configs produced only "
                  f"{len(unique_dmap)} distinct damage-map images")
        if len(unique_kd) < len(hashes):
            print(f"  ⚠ {group}: {len(hashes)} configs produced only "
                  f"{len(unique_kd)} distinct knockdown-curve images")

    print(f"\nRendered {len(all_reports)} cases into {OUT}/")


if __name__ == "__main__":
    main()
