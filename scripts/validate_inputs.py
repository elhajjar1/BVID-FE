#!/usr/bin/env python3
"""Input-variable validation sweep for BVID-FE.

Holds layup / ply thickness / panel dimensions / material constant and
sweeps each of: boundary, energy, impactor diameter, impactor shape,
impactor mass. For each level, runs empirical / semi_analytical / fe3d
under both compression and tension. Also runs a small energy sweep per
(boundary, tier, loading) combo.

Output: prints compact tables. Flags any variable whose level changes
produce zero variation in the downstream knockdown.

Usage:
    ./.venv/bin/python scripts/validate_inputs.py
"""
from __future__ import annotations

import sys
import warnings
from dataclasses import asdict
from itertools import product

import numpy as np

# Make stdout line-buffered so progress is visible as the sweep runs.
try:
    sys.stdout.reconfigure(line_buffering=True)  # py3.7+
except AttributeError:
    pass

from bvidfe.analysis import AnalysisConfig, BvidAnalysis, MeshParams
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.mapping import ImpactEvent, impact_to_damage
from bvidfe.impact.olsson import onset_energy, threshold_load

warnings.simplefilter("ignore")

# ---- fixed inputs ----
MATERIAL = "IM7/8552"
LAYUP = [0, 45, -45, 90, 90, -45, 45, 0]
PLY_T = 0.152
PANEL_LX, PANEL_LY = 150.0, 100.0

# ---- nominal impact (for single-point sweeps) ----
# Nominal energy chosen so DPA stays *under* the 80% panel-area cap at the
# default 150x100 panel — otherwise the cap saturates every prediction and
# all per-variable sensitivity disappears.
NOMINAL_ENERGY_J = 10.0
NOMINAL_DIAM_MM = 16.0
NOMINAL_SHAPE = "hemispherical"
NOMINAL_MASS_KG = 5.5
NOMINAL_BOUNDARY = "simply_supported"

# ---- levels ----
BOUNDARIES = ["simply_supported", "clamped", "free"]
ENERGIES = [3.0, 5.0, 8.0, 12.0, 20.0]
DIAMETERS = [8.0, 12.0, 16.0, 25.0, 40.0]
SHAPES = ["hemispherical", "flat", "conical"]
MASSES = [1.0, 2.5, 5.5, 10.0, 20.0]

TIERS = ["empirical", "semi_analytical", "fe3d"]
LOADINGS = ["compression", "tension"]

# Coarse fe3d mesh for the validation script only — we care about trends,
# not absolute values. Halves the per-run cost vs. the default.
FE3D_MESH = MeshParams(elements_per_ply=1, in_plane_size_mm=10.0)


def _make_cfg(
    *,
    boundary=NOMINAL_BOUNDARY,
    energy=NOMINAL_ENERGY_J,
    diameter=NOMINAL_DIAM_MM,
    shape=NOMINAL_SHAPE,
    mass=NOMINAL_MASS_KG,
    tier="empirical",
    loading="compression",
):
    return AnalysisConfig(
        material=MATERIAL,
        layup_deg=LAYUP,
        ply_thickness_mm=PLY_T,
        panel=PanelGeometry(PANEL_LX, PANEL_LY, boundary=boundary),
        loading=loading,
        tier=tier,
        impact=ImpactEvent(
            energy_J=energy,
            impactor=ImpactorGeometry(diameter_mm=diameter, shape=shape),
            mass_kg=mass,
        ),
        mesh=FE3D_MESH if tier == "fe3d" else None,
    )


def _run(cfg):
    try:
        r = BvidAnalysis(cfg).run()
        n_ell = len(r.damage.delaminations)
        dpa = r.damage.projected_damage_area_mm2
        dent = r.damage.dent_depth_mm
        fbr = r.damage.fiber_break_radius_mm
        return dict(
            ok=True,
            kd=r.knockdown,
            resid=r.residual_strength_MPa,
            dpa=dpa,
            dent=dent,
            fbr=fbr,
            n_ell=n_ell,
        )
    except Exception as exc:  # pragma: no cover
        return dict(ok=False, err=f"{type(exc).__name__}: {exc}")


def _impact_inspect(cfg):
    """Pre-compute Olsson Pc, E_onset, and DamageState without running a tier."""
    lam = Laminate(MATERIAL_LIBRARY[MATERIAL], LAYUP, PLY_T)
    pc = threshold_load(lam, cfg.panel, cfg.impact.impactor)
    eo = onset_energy(
        lam,
        cfg.panel,
        cfg.impact.impactor,
        location_xy_mm=(cfg.panel.Lx_mm / 2, cfg.panel.Ly_mm / 2),
    )
    ds = impact_to_damage(cfg.impact, lam, cfg.panel)
    return pc, eo, ds


def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


def table_header(columns):
    print("  ".join(columns))
    print("  ".join("-" * len(c) for c in columns))


def fmt_row(values):
    parts = []
    for v in values:
        if isinstance(v, float):
            parts.append(f"{v:>10.3f}")
        elif isinstance(v, int):
            parts.append(f"{v:>10d}")
        else:
            parts.append(f"{str(v):>10s}")
    print("  ".join(parts))


def sweep_one_variable(var_name, levels, set_kw):
    """Sweep one input variable across all (tier, loading) combos.
    Returns a flags list for any (tier, loading) combo where knockdown did not vary."""
    section(f"SWEEP {var_name}  (levels={levels})")

    # Pre-stage: how the impact-mapping stage responds to just this variable
    print(f"\n[impact-mapping response to {var_name}]")
    print(f"{'level':>18}  {'Pc [N]':>12}  {'E_onset [J]':>12}  {'DPA [mm^2]':>12}  "
          f"{'dent [mm]':>10}  {'fbr [mm]':>8}  {'#ell':>5}")
    for lvl in levels:
        cfg = _make_cfg(**{set_kw: lvl})
        pc, eo, ds = _impact_inspect(cfg)
        print(
            f"{str(lvl):>18}  {pc:>12.1f}  {eo:>12.3f}  "
            f"{ds.projected_damage_area_mm2:>12.1f}  {ds.dent_depth_mm:>10.3f}  "
            f"{ds.fiber_break_radius_mm:>8.3f}  {len(ds.delaminations):>5d}"
        )

    # Tier x loading table
    flags = []
    for tier in TIERS:
        for loading in LOADINGS:
            print(f"\n[{tier} / {loading}]")
            print(f"{'level':>18}  {'knockdown':>10}  {'residual_MPa':>14}  "
                  f"{'DPA':>10}  {'dent':>8}  {'n_ell':>6}")
            kds = []
            for lvl in levels:
                cfg = _make_cfg(tier=tier, loading=loading, **{set_kw: lvl})
                res = _run(cfg)
                if res["ok"]:
                    kds.append(res["kd"])
                    print(
                        f"{str(lvl):>18}  {res['kd']:>10.4f}  {res['resid']:>14.1f}  "
                        f"{res['dpa']:>10.1f}  {res['dent']:>8.3f}  {res['n_ell']:>6d}"
                    )
                else:
                    print(f"{str(lvl):>18}  ERR: {res['err']}")
                    kds.append(None)
            # flag if all knockdowns are identical across varied levels
            finite = [k for k in kds if k is not None]
            if len(finite) >= 2 and max(finite) - min(finite) < 1e-6:
                flags.append((var_name, tier, loading, finite[0]))
    return flags


def energy_sweep_per_combo(boundary):
    """Run an energy sweep (5..80 J) for every (tier, loading) combination."""
    section(f"ENERGY SWEEP @ boundary={boundary}")
    print(f"{'tier':>16} {'loading':>12}  " + "  ".join(f"{e:>7.1f}J" for e in ENERGIES))
    print("-" * (30 + 10 * len(ENERGIES)))
    for tier in TIERS:
        for loading in LOADINGS:
            kds = []
            for e in ENERGIES:
                cfg = _make_cfg(boundary=boundary, energy=e, tier=tier, loading=loading)
                res = _run(cfg)
                kds.append(res["kd"] if res["ok"] else float("nan"))
            kds_str = "  ".join(f"{k:>8.4f}" for k in kds)
            print(f"{tier:>16} {loading:>12}  {kds_str}")


def main():
    section("FIXED INPUTS")
    print(f"material:      {MATERIAL}")
    print(f"layup_deg:     {LAYUP}")
    print(f"ply_thickness: {PLY_T} mm  (total h = {len(LAYUP)*PLY_T:.3f} mm)")
    print(f"panel:         {PANEL_LX} x {PANEL_LY} mm")

    flags_all = []
    flags_all += sweep_one_variable("boundary", BOUNDARIES, "boundary")
    flags_all += sweep_one_variable("energy_J", ENERGIES, "energy")
    flags_all += sweep_one_variable("diameter_mm", DIAMETERS, "diameter")
    flags_all += sweep_one_variable("shape", SHAPES, "shape")
    flags_all += sweep_one_variable("mass_kg", MASSES, "mass")

    # Energy sweeps per boundary (to cross-check tier/boundary coupling)
    for bnd in BOUNDARIES:
        energy_sweep_per_combo(bnd)

    # Flags summary
    section("FLAGS — variables with zero effect on knockdown (per tier/loading)")
    if not flags_all:
        print("(none — every swept variable produced at least some variation)")
    else:
        print("The following (variable, tier, loading) tuples produced identical "
              "knockdowns across all levels — meaning the variable has no effect "
              "in that pipeline:\n")
        for var, tier, loading, kd in flags_all:
            print(f"  - {var:<14} @ {tier:<16} / {loading:<11}  kd = {kd:.4f}")


if __name__ == "__main__":
    main()
