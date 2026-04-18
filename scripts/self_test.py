#!/usr/bin/env python3
"""BVID-FE installation self-test.

Runs one analysis per tier on a standard configuration and prints a
summary. Exits 0 on success, 1 on any failure. Useful after a fresh
install to confirm everything is wired correctly.

Usage:
    python scripts/self_test.py
    python scripts/self_test.py --include-fe3d   # slow; adds ~10 s

Designed to be runnable from the repo root with `./.venv/bin/python`.
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings

from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


BASE_KW = dict(
    material="IM7/8552",
    layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
    ply_thickness_mm=0.152,
    panel=PanelGeometry(150, 100),
    loading="compression",
    impact=ImpactEvent(20.0, ImpactorGeometry(), mass_kg=5.5),
)


def run_tier(tier: str) -> tuple[bool, float, float, float]:
    """Run one analysis on the given tier. Returns (ok, KD, residual, seconds)."""
    cfg_kw = dict(BASE_KW)
    cfg_kw["tier"] = tier
    if tier == "fe3d":
        cfg_kw["mesh"] = MeshParams()  # default conservative mesh
    cfg = AnalysisConfig(**cfg_kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0 = time.time()
        try:
            result = BvidAnalysis(cfg).run()
            return True, result.knockdown, result.residual_strength_MPa, time.time() - t0
        except Exception as exc:
            print(f"  ERROR in {tier}: {type(exc).__name__}: {exc}")
            return False, float("nan"), float("nan"), time.time() - t0


def main() -> int:
    parser = argparse.ArgumentParser(description="BVID-FE installation self-test.")
    parser.add_argument(
        "--include-fe3d",
        action="store_true",
        help="Also run the fe3d tier (adds ~10 s).",
    )
    args = parser.parse_args()

    print("BVID-FE self-test")
    print("=================")
    print("Config: IM7/8552 quasi-iso [0/45/-45/90]_2s, 150x100 mm, 20 J impact")
    print()

    tiers = ["empirical", "semi_analytical"]
    if args.include_fe3d:
        tiers.append("fe3d")

    all_ok = True
    results = []
    for tier in tiers:
        print(f"  Running {tier}...", flush=True)
        ok, kd, residual, secs = run_tier(tier)
        results.append((tier, ok, kd, residual, secs))
        all_ok = all_ok and ok

    print()
    print(f"{'Tier':<18} {'OK?':<4} {'Knockdown':<10} {'Residual [MPa]':<16} {'Time [s]':<8}")
    print(f"{'':-<18} {'':-<4} {'':-<10} {'':-<16} {'':-<8}")
    for tier, ok, kd, res, secs in results:
        mark = "OK" if ok else "FAIL"
        print(f"{tier:<18} {mark:<4} {kd:<10.3f} {res:<16.1f} {secs:<8.2f}")

    print()
    if all_ok:
        print("All tiers completed successfully. BVID-FE install is working.")
        return 0
    else:
        print("One or more tiers failed. See errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
