"""BVID-FE command-line interface.

Runs a BvidAnalysis from command-line arguments and prints the result as JSON.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import replace
from typing import List, Sequence

import bvidfe
from bvidfe._types import _LOADING_MODES, _TIER_NAMES
from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.analysis.fe_tier import FE3DSizeError, enable_stderr_logging
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent


def _parse_panel(spec: str) -> PanelGeometry:
    try:
        a, b = spec.lower().split("x")
        Lx, Ly = float(a), float(b)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--panel must be '<Lx>x<Ly>' (got {spec!r})") from exc
    if not (math.isfinite(Lx) and math.isfinite(Ly)) or Lx <= 0 or Ly <= 0:
        raise argparse.ArgumentTypeError(
            f"--panel dimensions must be positive and finite (got {Lx}x{Ly})"
        )
    return PanelGeometry(Lx_mm=Lx, Ly_mm=Ly)


def _parse_layup(spec: str) -> List[float]:
    """Comma-separated ply angles, each finite and within [-90, 90] degrees
    (the same range the Streamlit app accepts)."""
    try:
        angles = [float(x) for x in spec.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--layup must be a comma-separated list of ply angles in degrees"
        ) from exc
    for i, a in enumerate(angles):
        if not math.isfinite(a):
            raise argparse.ArgumentTypeError(f"ply angle {i} must be finite (got {a})")
        if not -90.0 <= a <= 90.0:
            raise argparse.ArgumentTypeError(
                f"ply angle {i} must be within [-90, 90] degrees (got {a:g})"
            )
    return angles


def _positive_float(spec: str) -> float:
    """argparse type that rejects non-positive and non-finite numbers."""
    try:
        v = float(spec)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected a number (got {spec!r})") from exc
    if not math.isfinite(v):
        raise argparse.ArgumentTypeError(f"must be finite (got {v})")
    if v <= 0:
        raise argparse.ArgumentTypeError(f"must be > 0 (got {v})")
    return v


def _parse_thickness(spec: str):
    """argparse type for ``--thickness``.

    Accepts either a single positive number (uniform ply thickness) or a
    comma-separated list of positive numbers (per-ply thicknesses, length
    must match the layup at config-validation time).
    """
    if "," in spec:
        try:
            ts = [float(x) for x in spec.split(",")]
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                "--thickness must be a positive number or a comma-separated "
                "list of positive numbers (got " + repr(spec) + ")"
            ) from exc
        for i, t in enumerate(ts):
            if not math.isfinite(t) or t <= 0:
                raise argparse.ArgumentTypeError(
                    f"--thickness[{i}] must be > 0 and finite (got {t})"
                )
        return ts
    return _positive_float(spec)


def _existing_path(spec: str):
    """argparse type that requires a readable file at the given path."""
    from pathlib import Path

    p = Path(spec)
    if not p.exists():
        raise argparse.ArgumentTypeError(f"file not found: {spec}")
    if not p.is_file():
        raise argparse.ArgumentTypeError(f"not a file: {spec}")
    return p


_VALID_TIERS = tuple(sorted(_TIER_NAMES))


def _parse_tiers(spec: str) -> List[str]:
    """Parse ``--tier`` as one tier or a comma-separated list to compare.

    Validated against the known tiers, order-preserving, and de-duplicated
    so ``--tier empirical,empirical`` doesn't run the same solve twice.
    """
    tiers: List[str] = []
    for raw in spec.split(","):
        t = raw.strip()
        if not t:
            continue
        if t not in _VALID_TIERS:
            raise argparse.ArgumentTypeError(
                f"unknown tier {t!r}; choose from {', '.join(_VALID_TIERS)}"
            )
        if t not in tiers:
            tiers.append(t)
    if not tiers:
        raise argparse.ArgumentTypeError("--tier must name at least one tier")
    return tiers


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show ``(default: ...)`` only for options that have a real default,
    not ``None`` (required inputs) or ``False`` (flags)."""

    def _get_help_string(self, action: argparse.Action) -> str:
        if action.default is None or action.default is False:
            return action.help or ""
        return super()._get_help_string(action) or ""


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bvidfe",
        description="BVID-FE: Barely Visible Impact Damage residual-strength analysis.",
        formatter_class=_HelpFormatter,
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"bvidfe {bvidfe.__version__}",
    )
    # Defer the MATERIAL_LIBRARY import + choices population to parser-build
    # time so --material rejects typos early with the standard argparse
    # 'invalid choice' message listing the presets, instead of failing
    # downstream with a raw KeyError.
    from bvidfe.core.material import MATERIAL_LIBRARY

    p.add_argument(
        "--material",
        choices=sorted(MATERIAL_LIBRARY.keys()),
        help="Material preset name (e.g. IM7/8552)",
    )
    p.add_argument(
        "--layup",
        type=_parse_layup,
        help="Comma-separated ply angles in degrees, e.g. 0,45,-45,90",
    )
    p.add_argument(
        "--thickness",
        type=_parse_thickness,
        help="Ply thickness in millimeters. Either a single positive number "
        "(uniform laminate) or a comma-separated list of per-ply thicknesses "
        "with length equal to the number of plies in --layup, e.g. "
        "'0.10,0.10,0.20,0.20,0.20,0.20,0.10,0.10' for a hybrid stack.",
    )
    p.add_argument(
        "--panel",
        type=_parse_panel,
        help="Panel dimensions as 'Lx_mm x Ly_mm' in millimeters, "
        "e.g. 150x100. Lowercase 'x' separator, no spaces.",
    )
    p.add_argument(
        "--loading",
        choices=sorted(_LOADING_MODES),
        help="Load case: compression-after-impact or tension-after-impact",
    )
    p.add_argument(
        "--tier",
        type=_parse_tiers,
        # A string default is run through _parse_tiers, giving ["empirical"].
        default="empirical",
        help="Analysis fidelity tier(s). A single tier (e.g. semi_analytical) "
        "or a comma-separated list to compare in one invocation (e.g. "
        f"empirical,semi_analytical,fe3d). Choices: {', '.join(_VALID_TIERS)}.",
    )
    p.add_argument(
        "--energy",
        type=_positive_float,
        help="Impact energy in Joules (impact-driven workflow). Mutually exclusive with --cscan.",
    )
    p.add_argument(
        "--cscan",
        type=_existing_path,
        help="Path to a C-scan JSON file (inspection-driven workflow). "
        "Mutually exclusive with --energy. See docs/cscan_schema.md for the format.",
    )
    p.add_argument(
        "--impactor-diameter",
        type=_positive_float,
        default=16.0,
        help="Impactor diameter in mm",
    )
    p.add_argument("--mass", type=_positive_float, default=5.5, help="Impactor mass in kg")
    p.add_argument(
        "--quick",
        action="store_true",
        help="Print only the knockdown scalar (residual / pristine) to stdout instead of the full JSON. "
        "Useful for shell pipelines: e.g. `bvidfe ... --quick | xargs -I {} ...`.",
    )
    p.add_argument(
        "--quick-json",
        action="store_true",
        help="Like --quick but emits a single-line JSON object {knockdown, "
        "residual_strength_MPa, pristine_strength_MPa, tier_used} instead of "
        "a bare scalar. Recommended for scripted callers that want the "
        "value plus tier provenance without parsing the full --json output.",
    )
    p.add_argument(
        "--list-materials",
        action="store_true",
        help="List available material presets with key properties and exit.",
    )
    return p


def _list_materials() -> None:
    from bvidfe.core.material import MATERIAL_LIBRARY

    print(f"{'Name':<18} {'E11':>8} {'E22':>7} {'Xt':>7} {'Xc':>7} {'Yt':>5} {'Yc':>5}")
    print(
        f"{'':-<18} {'-' * 8:>8} {'-' * 7:>7} {'-' * 7:>7} {'-' * 7:>7} {'-' * 5:>5} {'-' * 5:>5}"
    )
    for name, m in MATERIAL_LIBRARY.items():
        print(
            f"{name:<18} {m.E11:>8.0f} {m.E22:>7.0f} {m.Xt:>7.0f} {m.Xc:>7.0f} {m.Yt:>5.0f} {m.Yc:>5.0f}"
        )
    print("\nUnits: MPa. Use --material <Name> to select.")


def _error(msg: str) -> int:
    """Print an argparse-style error for a failure after parsing; exit code 1."""
    sys.stderr.write(f"bvidfe: error: {msg}\n")
    return 1


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return _main(argv)
    except BrokenPipeError:
        # The reader went away (``bvidfe ... | head``). Point stdout at
        # devnull so the interpreter's final flush does not raise again; see
        # the "Note on SIGPIPE" in the Python ``signal`` docs.
        try:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
        except (OSError, ValueError):
            pass
        return 1


def _main(argv: Sequence[str] | None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.list_materials:
        _list_materials()
        return 0
    # Base required args (always needed)
    missing = [
        n
        for n in ("material", "layup", "thickness", "panel", "loading")
        if getattr(args, n) is None
    ]
    if missing:
        parser.error(f"missing required arguments: {', '.join('--' + m for m in missing)}")
    # Exactly one of --energy (impact-driven) or --cscan (inspection-driven)
    if args.energy is None and args.cscan is None:
        parser.error("must provide either --energy (impact-driven) or --cscan (inspection-driven)")
    if args.energy is not None and args.cscan is not None:
        parser.error("--energy and --cscan are mutually exclusive")

    # If --thickness was given as a per-ply list, its length must match the layup.
    if isinstance(args.thickness, list) and len(args.thickness) != len(args.layup):
        parser.error(
            f"--thickness has {len(args.thickness)} per-ply values but --layup "
            f"has {len(args.layup)} plies; they must match."
        )

    if args.quick and args.quick_json:
        parser.error("--quick and --quick-json are mutually exclusive")

    damage = None
    if args.cscan is not None:
        from bvidfe.damage.io import CScanSchemaError, load_cscan_json

        # args.cscan is already a validated Path (see _existing_path).
        try:
            damage = load_cscan_json(args.cscan)
        except CScanSchemaError as exc:
            parser.error(f"--cscan {args.cscan}: invalid C-scan file: {exc}")
        except (OSError, UnicodeDecodeError) as exc:
            parser.error(f"--cscan {args.cscan}: cannot read file: {exc}")
    try:
        cfg = AnalysisConfig(
            material=args.material,
            layup_deg=args.layup,
            ply_thickness_mm=args.thickness,
            panel=args.panel,
            loading=args.loading,
            tier=args.tier[0],
            impact=(
                None
                if args.energy is None
                else ImpactEvent(
                    energy_J=args.energy,
                    impactor=ImpactorGeometry(diameter_mm=args.impactor_diameter),
                    mass_kg=args.mass,
                )
            ),
            damage=damage,
        )
    except ValueError as exc:
        parser.error(str(exc))

    # fe3d stage logs go to stderr (the library itself only has a NullHandler).
    enable_stderr_logging()

    # One result per requested tier. Single-tier output is byte-identical to
    # the pre-multi-tier behaviour; multiple tiers switch to a comparison
    # shape (JSON array / NDJSON / TSV) keyed by tier_used.
    try:
        results = [BvidAnalysis(replace(cfg, tier=t)).run() for t in args.tier]
    except FE3DSizeError as exc:
        return _error(str(exc))
    multi = len(results) > 1

    if args.quick_json:
        for result in results:
            json.dump(
                {
                    "knockdown": result.knockdown,
                    "residual_strength_MPa": result.residual_strength_MPa,
                    "pristine_strength_MPa": result.pristine_strength_MPa,
                    "tier_used": result.tier_used,
                },
                sys.stdout,
            )
            sys.stdout.write("\n")
    elif args.quick:
        for result in results:
            if multi:
                print(f"{result.tier_used}\t{result.knockdown:.6f}")
            else:
                print(f"{result.knockdown:.6f}")
    else:
        if multi:
            json.dump([r.to_dict() for r in results], sys.stdout, indent=2, default=str)
        else:
            json.dump(results[0].to_dict(), sys.stdout, indent=2, default=str)
        sys.stdout.write("\n")
    # Flush here so a closed pipe raises inside main()'s BrokenPipeError guard.
    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
