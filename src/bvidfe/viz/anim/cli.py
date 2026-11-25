"""``bvidfe-animate``: render the BVID-FE animations from the command line.

Three subcommands, ``energy-sweep``, ``damage-stack`` and ``cai-ramp``,
take the analysis options of ``bvidfe`` (``--material``, ``--layup``,
``--thickness``, ``--panel``, ``--energy`` or ``--cscan``, ...; ``cai-ramp``
also the fe3d mesh options ``--in-plane-size`` and ``--elements-per-ply``)
and the output options (``--out``, ``--format``, ``--fps``, ``--duration``,
``--width`` / ``--height``, ``--preview``).

Examples::

    bvidfe-animate energy-sweep --material IM7/8552 --layup 0,45,-45,90,0,45,-45,90 \\
        --thickness 0.152 --panel 150x100 --energy-max 30 --out damage.mp4
    bvidfe-animate damage-stack --material IM7/8552 --layup 0,45,-45,90,90,-45,45,0 \\
        --thickness 0.152 --panel 150x100 --energy 20 --preview --out stack.png
    bvidfe-animate cai-ramp --material IM7/8552 --layup 0,45,-45,90,90,-45,45,0 \\
        --thickness 0.152 --panel 150x100 --energy 20 --in-plane-size 5 --out ramp.mp4
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import bvidfe
from bvidfe.cli import (
    _HelpFormatter,
    _add_config_arguments,
    _check_config_args,
    _config_from_args,
    _positive_float,
    _positive_int,
)
from bvidfe.viz.anim.spec import AnimationSpec

PROG = "bvidfe-animate"
_SUBCOMMANDS = ("energy-sweep", "damage-stack", "cai-ramp")


def _non_negative(spec: str) -> float:
    try:
        v = float(spec)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected a number (got {spec!r})") from exc
    if not math.isfinite(v) or v < 0:
        raise argparse.ArgumentTypeError(f"must be finite and >= 0 (got {v})")
    return v


def _add_output_arguments(p: argparse.ArgumentParser, default_out: str) -> None:
    g = p.add_argument_group("output")
    g.add_argument(
        "--out",
        default=default_out,
        metavar="PATH",
        help="Output file. A PNG sequence goes into a directory: PATH without "
        "its suffix plus '_frames' (or PATH itself when it has no suffix). With "
        "--preview, the stills are PATH's stem plus '_t0.00.png' ...",
    )
    g.add_argument(
        "--format",
        choices=("mp4", "gif", "png", "auto"),
        default="mp4",
        help="mp4 (H.264, needs ffmpeg), gif (a reduced-size preview), png "
        "(one file per frame) or auto (mp4 if ffmpeg is found, else png)",
    )
    g.add_argument("--fps", type=_positive_int, help="Frames per second (default: 30; 12 for gif)")
    g.add_argument(
        "--duration",
        type=_positive_float,
        default=10.0,
        metavar="S",
        help="Clip length in seconds",
    )
    g.add_argument(
        "--width",
        type=_positive_int,
        help="Frame width in px, 16:9 with --height (default: 1920; 640 for gif)",
    )
    g.add_argument(
        "--height",
        type=_positive_int,
        help="Frame height in px (default: 1080; 360 for gif)",
    )
    g.add_argument(
        "--preview",
        action="store_true",
        help="Render three stills (t = 0, 0.5 and 1) as PNG instead of the clip",
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=PROG,
        description="Render BVID-FE animations (1920 x 1080 MP4 by default). Every frame "
        "carries the citation of BVID-FE.",
        formatter_class=_HelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"{PROG} {bvidfe.__version__}")
    sub = p.add_subparsers(dest="scene", metavar="SCENE")
    sub.required = True
    es = sub.add_parser(
        "energy-sweep",
        help="Damage growing with impact energy, with area, dent and knockdown curves",
        description="The delaminations (plan view), damage area, dent depth and "
        "knockdown as the impact energy rises: one analysis per keyframe energy "
        "(empirical or semi_analytical tier).",
        formatter_class=_HelpFormatter,
    )
    _add_config_arguments(es, mesh=False)
    g = es.add_argument_group("energy sweep")
    g.add_argument(
        "--energy-min",
        type=_non_negative,
        metavar="J",
        help="First keyframe energy in J (default: 0.8 x the onset energy)",
    )
    g.add_argument(
        "--energy-max",
        type=_positive_float,
        metavar="J",
        help="Last keyframe energy in J (default: --energy)",
    )
    g.add_argument(
        "--n-energies",
        type=_positive_int,
        metavar="N",
        help="Number of keyframe energies (default: 60 empirical, 40 semi_analytical)",
    )
    _add_output_arguments(es, "energy_sweep.mp4")
    es.set_defaults(subparser=es)
    ds = sub.add_parser(
        "damage-stack",
        help="Exploded 3D stack of the delaminations, camera orbiting",
        description="The delaminations through the thickness as an exploded stack "
        "of plies, the camera orbiting 360 deg; damage from the impact model "
        "(--energy) or a C-scan (--cscan).",
        formatter_class=_HelpFormatter,
    )
    _add_config_arguments(ds, tier=False, mesh=False)
    ds.add_argument(
        "--no-intro",
        action="store_true",
        help="Leave out the plies dropping into place, for a loopable orbit",
    )
    _add_output_arguments(ds, "damage_stack.mp4")
    ds.set_defaults(subparser=ds)
    cr = sub.add_parser(
        "cai-ramp",
        help="fe3d load ramp to failure, elements coloured by their load ratio",
        description="The fe3d load ramp of the damaged panel to its residual "
        "strength: the panel in 3D, cut open through the first element to fail, "
        "each element coloured by its load ratio to local failure, beside the "
        "stress against end shortening. Two fe3d solves (four for a biaxial "
        "load), always fe3d whatever the analysis tier; linear elastic to "
        "first-ply failure.",
        formatter_class=_HelpFormatter,
    )
    _add_config_arguments(cr, tier=False)
    cr.add_argument(
        "--pointwise",
        action="store_true",
        help="Colour by each element's lowest Gauss point instead of the "
        "disk-averaged check that sets first-ply failure (it can exceed 1 first)",
    )
    _add_output_arguments(cr, "cai_ramp.mp4")
    cr.set_defaults(subparser=cr)
    return p


def _spec(parser: argparse.ArgumentParser, args: argparse.Namespace) -> AnimationSpec:
    try:
        return AnimationSpec(
            width=args.width,
            height=args.height,
            fps=args.fps,
            duration_s=args.duration,
            format=args.format,
        )
    except ValueError as exc:
        parser.error(str(exc))


def _progress(i: int, n: int) -> None:
    step = max(1, n // 10)
    if i == n or i % step == 0:
        sys.stderr.write(f"\r{PROG}: frame {i}/{n}")
        if i == n:
            sys.stderr.write("\n")
        sys.stderr.flush()


def energy_sweep_tiers() -> tuple:
    from bvidfe.viz.anim.scenes.energy_sweep import SUPPORTED_TIERS

    return SUPPORTED_TIERS


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    sub: argparse.ArgumentParser = args.subparser
    if args.scene == "energy-sweep":
        if args.cscan is not None:
            sub.error("--cscan has no impact energy to sweep; use damage-stack for a C-scan")
        if args.energy is None:
            if args.energy_max is None:
                sub.error("give --energy or --energy-max (the last energy of the sweep)")
            args.energy = args.energy_max
        if len(args.tier) > 1:
            sub.error("--tier takes one tier for an animation")
        if args.tier[0] not in energy_sweep_tiers():
            sub.error(
                f"--tier {args.tier[0]}: the energy sweep runs one analysis per keyframe "
                f"energy; use {' or '.join(energy_sweep_tiers())}"
            )
    if args.scene == "cai-ramp":
        args.tier = ["fe3d"]  # the ramp is fe3d's; the mesh options apply
    _check_config_args(sub, args)
    spec = _spec(sub, args)
    cfg = _config_from_args(sub, args)

    from bvidfe.viz.anim import encode_timeline, preview_frames
    from bvidfe.analysis.fe_tier import FE3DSizeError
    from bvidfe.viz.anim.scenes import cai_ramp, damage_stack, energy_sweep
    from bvidfe.viz.anim.timeline import Timeline

    tl: Timeline
    try:
        if args.scene == "energy-sweep":
            tl = energy_sweep.build_timeline(
                cfg,
                energy_min=args.energy_min,
                energy_max=args.energy_max,
                n_energies=args.n_energies,
            )
        elif args.scene == "damage-stack":
            tl = damage_stack.build_timeline(cfg, intro=not args.no_intro)
        else:
            tl = cai_ramp.build_timeline(cfg, ratio="pointwise" if args.pointwise else "averaged")
    except (ValueError, FE3DSizeError) as exc:
        sys.stderr.write(f"{PROG}: error: {exc}\n")
        return 1

    out = Path(args.out)
    if args.preview:
        times = (0.0, 0.5, 1.0)
        directory = out.parent
        preview_frames(tl, times, spec, out_dir=directory, stem=out.stem)
        paths: List[Path] = [directory / f"{out.stem}_t{t:.2f}.png" for t in times]
        for p in paths:
            print(p)
        return 0
    from bvidfe.viz.anim.writer import FFmpegNotFoundError

    try:
        result = encode_timeline(tl, out, spec, progress=_progress)
    except FFmpegNotFoundError as exc:
        sys.stderr.write(f"{PROG}: error: {exc}\n")
        return 1
    size_mb = result.size_bytes / 1e6
    sys.stderr.write(
        f"{PROG}: wrote {result.path} ({result.format}, {result.width} x {result.height}, "
        f"{result.n_frames} frames at {result.fps} fps, {size_mb:.1f} MB) in "
        f"{result.elapsed_s:.1f} s\n"
    )
    print(result.path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
