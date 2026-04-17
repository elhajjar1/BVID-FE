"""QThread workers that run BvidAnalysis off the UI thread."""

from __future__ import annotations

import logging
import threading
import time
import traceback
from dataclasses import replace
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
from PyQt6.QtCore import QThread, pyqtSignal

from bvidfe.analysis import AnalysisConfig, BvidAnalysis

_log = logging.getLogger("bvidfe.gui")


class AnalysisWorker(QThread):
    """Runs `BvidAnalysis(config).run()` in a background thread.

    BvidAnalysis is synchronous and can take tens of seconds on the fe3d
    tier. To keep the status bar from appearing frozen, we run the
    analysis in a daemon worker thread and emit heartbeat progress from
    the QThread itself every few seconds until the work thread finishes.
    """

    resultReady = pyqtSignal(object)  # AnalysisResults
    error = pyqtSignal(str)
    progress = pyqtSignal(int)

    HEARTBEAT_INTERVAL_S: float = 2.0

    def __init__(self, config: AnalysisConfig, parent=None) -> None:
        super().__init__(parent)
        self.config = config

    def run(self) -> None:  # type: ignore[override]
        result_box: list[object] = [None]
        error_box: list[str] = []
        t_start = time.time()
        _log.info(
            "AnalysisWorker started: tier=%s loading=%s", self.config.tier, self.config.loading
        )

        def _do_work() -> None:
            try:
                result_box[0] = BvidAnalysis(self.config).run()
            except Exception:  # noqa: BLE001
                error_box.append(traceback.format_exc())

        t = threading.Thread(target=_do_work, daemon=True)
        t.start()

        # Tick heartbeat 10% -> 90% while waiting for the work thread.
        # Each tick bumps the progress until we hit the cap; the actual
        # 100% comes from the completion branch below.
        pct = 10
        self.progress.emit(pct)
        while t.is_alive():
            t.join(timeout=self.HEARTBEAT_INTERVAL_S)
            if t.is_alive() and pct < 90:
                pct = min(90, pct + 5)
                self.progress.emit(pct)
                _log.info("AnalysisWorker heartbeat: %d%% (%.1fs)", pct, time.time() - t_start)

        if error_box:
            _log.warning(
                "AnalysisWorker error after %.1fs: %s",
                time.time() - t_start,
                error_box[0].splitlines()[-1],
            )
            self.error.emit(error_box[0])
            return
        _log.info("AnalysisWorker done (%.1fs)", time.time() - t_start)
        self.progress.emit(100)
        self.resultReady.emit(result_box[0])


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

                # Run this energy in a daemon thread so we can heartbeat
                # progress between the coarse per-energy steps.
                result_box: list[object] = [None]
                err_box: list[str] = []

                def _do_one() -> None:
                    try:
                        result_box[0] = BvidAnalysis(cfg).run()
                    except Exception:  # noqa: BLE001
                        err_box.append(traceback.format_exc())

                th = threading.Thread(target=_do_one, daemon=True)
                th.start()
                base_pct = 5 + int(90 * i / n)
                next_pct = 5 + int(90 * (i + 1) / n)
                tick = base_pct
                while th.is_alive():
                    th.join(timeout=2.0)
                    if th.is_alive() and tick + 1 < next_pct:
                        tick += 1
                        self.progress.emit(tick)
                if err_box:
                    raise RuntimeError(err_box[0])
                result = result_box[0]
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
                self.progress.emit(next_pct)
            df = pd.DataFrame(rows)
            if self.csv_path is not None:
                df.to_csv(Path(self.csv_path), index=False)
            self.progress.emit(100)
            self.resultReady.emit(df)
        except Exception:
            self.error.emit(traceback.format_exc())
