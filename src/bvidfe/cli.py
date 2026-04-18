"""BVID-FE command-line interface.

Runs a BvidAnalysis from command-line arguments and prints the result as JSON.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import List, Sequence

import bvidfe
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


def _parse_panel(spec: str) -> PanelGeometry:
    try:
        a, b = spec.lower().split("x")
        return PanelGeometry(Lx_mm=float(a), Ly_mm=float(b))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--panel must be '<Lx>x<Ly>' (got {spec!r})") from exc


def _parse_layup(spec: str) -> List[float]:
    try:
        return [float(x) for x in spec.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--layup must be a comma-separated list of ply angles in degrees"
        ) from exc


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bvidfe",
        description="BVID-FE: Barely Visible Impact Damage residual-strength analysis.",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"bvidfe {bvidfe.__version__}",
    )
    p.add_argument("--material", required=True, help="Material preset name (e.g. IM7/8552)")
    p.add_argument(
        "--layup",
        required=True,
        type=_parse_layup,
        help="Comma-separated ply angles in degrees, e.g. 0,45,-45,90",
    )
    p.add_argument("--thickness", required=True, type=float, help="Ply thickness in millimeters")
    p.add_argument(
        "--panel",
        required=True,
        type=_parse_panel,
        help="Panel dimensions as LxY in millimeters, e.g. 150x100",
    )
    p.add_argument("--loading", required=True, choices=["compression", "tension"])
    p.add_argument("--tier", default="empirical", choices=["empirical", "semi_analytical", "fe3d"])
    p.add_argument("--energy", type=float, required=True, help="Impact energy in Joules")
    p.add_argument(
        "--impactor-diameter",
        type=float,
        default=16.0,
        help="Impactor diameter in mm (default 16.0)",
    )
    p.add_argument("--mass", type=float, default=5.5, help="Impactor mass in kg (default 5.5)")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    cfg = AnalysisConfig(
        material=args.material,
        layup_deg=args.layup,
        ply_thickness_mm=args.thickness,
        panel=args.panel,
        loading=args.loading,
        tier=args.tier,
        impact=ImpactEvent(
            energy_J=args.energy,
            impactor=ImpactorGeometry(diameter_mm=args.impactor_diameter),
            mass_kg=args.mass,
        ),
    )
    result = BvidAnalysis(cfg).run()
    json.dump(result.to_dict(), sys.stdout, indent=2, default=str)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
