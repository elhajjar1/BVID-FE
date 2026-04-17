"""Parametric sweep utilities producing pandas DataFrames and CSV output."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import List, Optional, Sequence

import pandas as pd

from bvidfe.analysis import AnalysisConfig, BvidAnalysis


def _run_one(cfg: AnalysisConfig) -> dict:
    """Run a single analysis and return a dict of key fields."""
    result = BvidAnalysis(cfg).run()
    return {
        "knockdown": result.knockdown,
        "residual_MPa": result.residual_strength_MPa,
        "pristine_MPa": result.pristine_strength_MPa,
        "dpa_mm2": result.dpa_mm2,
        "dent_mm": result.damage.dent_depth_mm,
        "n_delaminations": len(result.damage.delaminations),
        "tier_used": result.tier_used,
    }


def _write_csv(df: pd.DataFrame, csv_path: Optional[Path]) -> None:
    if csv_path is not None:
        df.to_csv(Path(csv_path), index=False)


def sweep_energies(
    base_cfg: AnalysisConfig,
    energies_J: Sequence[float],
    csv_path: Optional[Path | str] = None,
) -> pd.DataFrame:
    """Sweep impact energies; base_cfg must have `impact` set."""
    if base_cfg.impact is None:
        raise ValueError("sweep_energies requires base_cfg.impact to be set")
    rows: List[dict] = []
    for E in energies_J:
        new_impact = replace(base_cfg.impact, energy_J=float(E))
        cfg = replace(base_cfg, impact=new_impact)
        row = _run_one(cfg)
        row["energy_J"] = float(E)
        rows.append(row)
    df = pd.DataFrame(rows)
    _write_csv(df, Path(csv_path) if csv_path else None)
    return df


def sweep_layups(
    base_cfg: AnalysisConfig,
    layups: Sequence[Sequence[float]],
    csv_path: Optional[Path | str] = None,
) -> pd.DataFrame:
    """Sweep layup sequences."""
    rows: List[dict] = []
    for layup in layups:
        cfg = replace(base_cfg, layup_deg=list(layup))
        row = _run_one(cfg)
        row["layup"] = "/".join(f"{a:g}" for a in layup)
        rows.append(row)
    df = pd.DataFrame(rows)
    _write_csv(df, Path(csv_path) if csv_path else None)
    return df


def sweep_thicknesses(
    base_cfg: AnalysisConfig,
    ply_thicknesses_mm: Sequence[float],
    csv_path: Optional[Path | str] = None,
) -> pd.DataFrame:
    """Sweep ply thickness values."""
    rows: List[dict] = []
    for t in ply_thicknesses_mm:
        cfg = replace(base_cfg, ply_thickness_mm=float(t))
        row = _run_one(cfg)
        row["ply_thickness_mm"] = float(t)
        rows.append(row)
    df = pd.DataFrame(rows)
    _write_csv(df, Path(csv_path) if csv_path else None)
    return df
