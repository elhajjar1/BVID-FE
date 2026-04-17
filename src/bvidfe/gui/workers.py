"""QThread workers that run BvidAnalysis off the UI thread."""

from __future__ import annotations

import traceback
from pathlib import Path
from typing import Optional, Sequence

from PyQt6.QtCore import QThread, pyqtSignal

from bvidfe.analysis import AnalysisConfig, BvidAnalysis
from bvidfe.sweep.parametric_sweep import sweep_energies


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
            self.progress.emit(25)
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
        try:
            self.progress.emit(10)
            df = sweep_energies(self.base_config, self.energies_J, csv_path=self.csv_path)
            self.progress.emit(100)
            self.resultReady.emit(df)
        except Exception:
            self.error.emit(traceback.format_exc())
