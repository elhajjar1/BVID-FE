"""QThread workers that run BvidAnalysis off the UI thread."""

from __future__ import annotations

import traceback
from dataclasses import replace
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
from PyQt6.QtCore import QThread, pyqtSignal

from bvidfe.analysis import AnalysisConfig, BvidAnalysis


class AnalysisWorker(QThread):
    """Runs `BvidAnalysis(config).run()` in a background thread."""

    resultReady = pyqtSignal(object)  # AnalysisResults
    error = pyqtSignal(str)
    progress = pyqtSignal(int)

    def __init__(self, config: AnalysisConfig, parent=None) -> None:
        super().__init__(parent)
        self.config = config

    def run(self) -> None:  # type: ignore[override]
        try:
            self.progress.emit(10)
            result = BvidAnalysis(self.config).run()
            self.progress.emit(100)
            self.resultReady.emit(result)
        except Exception:
            self.error.emit(traceback.format_exc())


class SweepWorker(QThread):
    """Runs a parametric energy sweep in a background thread."""

    resultReady = pyqtSignal(object)  # pandas.DataFrame
    error = pyqtSignal(str)
    progress = pyqtSignal(int)

    def __init__(
        self,
        base_config: AnalysisConfig,
        energies_J: Sequence[float],
        csv_path: Optional[str | Path] = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.base_config = base_config
        self.energies_J = list(energies_J)
        self.csv_path = csv_path

    def run(self) -> None:  # type: ignore[override]
        """Run one BvidAnalysis per energy and report progress after each point.

        This in-lines the sweep_energies logic so we can emit per-energy
        progress signals instead of only 10% at the start and 100% at the end.
        The functional behavior is identical to sweep_energies (same dataframe
        columns and CSV format).
        """
        try:
            n = len(self.energies_J)
            rows: list[dict] = []
            if self.base_config.impact is None:
                raise ValueError("sweep requires base_config.impact to be set")
            self.progress.emit(5)
            for i, E in enumerate(self.energies_J):
                new_impact = replace(self.base_config.impact, energy_J=float(E))
                cfg = replace(self.base_config, impact=new_impact)
                result = BvidAnalysis(cfg).run()
                rows.append(
                    {
                        "energy_J": float(E),
                        "knockdown": result.knockdown,
                        "residual_MPa": result.residual_strength_MPa,
                        "pristine_MPa": result.pristine_strength_MPa,
                        "dpa_mm2": result.dpa_mm2,
                        "dent_mm": result.damage.dent_depth_mm,
                        "n_delaminations": len(result.damage.delaminations),
                        "tier_used": result.tier_used,
                    }
                )
                # Linear progress from 5 to 95 percent across the N energies
                self.progress.emit(5 + int(90 * (i + 1) / n))
            df = pd.DataFrame(rows)
            if self.csv_path is not None:
                df.to_csv(Path(self.csv_path), index=False)
            self.progress.emit(100)
            self.resultReady.emit(df)
        except Exception:
            self.error.emit(traceback.format_exc())
