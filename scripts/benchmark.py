#!/usr/bin/env python3
"""BVID-FE performance benchmark.

Times each tier across a small set of representative problems and prints
a table. Useful for:
- Verifying performance characteristics after a code change
- Picking realistic mesh sizes for your workflow
- Documenting runtime expectations in papers / reports

Usage:
    python scripts/benchmark.py
    python scripts/benchmark.py --fe3d-sizes 1,2,3   # more fe3d mesh points
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings

from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


def _base_cfg(tier: str, mesh_params: MeshParams | None = None) -> AnalysisConfig:
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier=tier,
        impact=ImpactEvent(20.0, ImpactorGeometry(), mass_kg=5.5),
        mesh=mesh_params,
    )


def time_run(cfg: AnalysisConfig) -> tuple[float, float]:
    """Return (elapsed_s, knockdown)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0 = time.time()
        result = BvidAnalysis(cfg).run()
        return time.time() - t0, result.knockdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fe3d-sizes",
        type=str,
        default="1",
        help="Comma-separated list of elements_per_ply values for fe3d benchmark rows "
        "(default '1').",
    )
    parser.add_argument(
        "--skip-fe3d",
        action="store_true",
        help="Skip all fe3d benchmark rows (keeps the bench under 2 s total).",
    )
    args = parser.parse_args()

    print("BVID-FE performance benchmark")
    print("=============================")
    print("Panel: 150x100 mm, IM7/8552 quasi-iso [0/45/-45/90]_2s, 20 J impact")
    print()

    rows: list[tuple[str, str, float, float]] = []

    # --- empirical ---
    elapsed, kd = time_run(_base_cfg("empirical"))
    rows.append(("empirical", "n/a", elapsed, kd))

    # --- semi_analytical ---
    elapsed, kd = time_run(_base_cfg("semi_analytical"))
    rows.append(("semi_analytical", "n/a", elapsed, kd))

    # --- fe3d across mesh sizes ---
    if not args.skip_fe3d:
        fe3d_epl = [int(s) for s in args.fe3d_sizes.split(",") if s.strip()]
        for epl in fe3d_epl:
            mp = MeshParams(elements_per_ply=epl, in_plane_size_mm=5.0)
            elapsed, kd = time_run(_base_cfg("fe3d", mp))
            rows.append(("fe3d", f"{epl} epl, 5 mm", elapsed, kd))

    print(f"{'Tier':<18} {'Mesh':<18} {'Time [s]':<10} {'Knockdown':<10}")
    print(f"{'':-<18} {'':-<18} {'':-<10} {'':-<10}")
    for tier, mesh, elapsed, kd in rows:
        print(f"{tier:<18} {mesh:<18} {elapsed:<10.2f} {kd:<10.3f}")

    print()
    print("Notes:")
    print("  - Empirical and semi_analytical are effectively instant (closed-form).")
    print("  - fe3d scales roughly linearly with element count in the assembly step,")
    print("    plus a cubic factor in the sparse LU of the solve step.")
    print("  - Use tier='empirical' for energy sweeps. fe3d is intended for")
    print("    stress-field context, not for high-throughput parametric studies.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
