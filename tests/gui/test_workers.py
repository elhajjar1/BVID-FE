import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from bvidfe.analysis import AnalysisConfig
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.gui.workers import AnalysisWorker, SweepWorker
from bvidfe.impact.mapping import ImpactEvent


@pytest.fixture
def cfg():
    return AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90, 90, -45, 45, 0],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(30.0, ImpactorGeometry(), mass_kg=5.5),
    )


def test_analysis_worker_emits_result_ready(qtbot, cfg):
    worker = AnalysisWorker(cfg)
    with qtbot.waitSignal(worker.resultReady, timeout=10_000) as blocker:
        worker.start()
    (result,) = blocker.args
    assert 0 < result.knockdown <= 1.0


def test_analysis_worker_emits_error_on_bad_config(qtbot):
    # Build a config that will raise at run time (tier not recognized)
    bad_cfg = AnalysisConfig(
        material="IM7/8552",
        layup_deg=[0, 90, 0, 90],
        ply_thickness_mm=0.152,
        panel=PanelGeometry(100, 50),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(10.0, ImpactorGeometry(), mass_kg=5.5),
    )
    # Corrupt it to force a failure inside run()
    bad_cfg.tier = "bogus"  # type: ignore[assignment]

    worker = AnalysisWorker(bad_cfg)
    with qtbot.waitSignal(worker.error, timeout=10_000):
        worker.start()


def test_sweep_worker_emits_result_ready(qtbot, cfg):
    worker = SweepWorker(cfg, energies_J=[5, 10], csv_path=None)
    with qtbot.waitSignal(worker.resultReady, timeout=20_000) as blocker:
        worker.start()
    (df,) = blocker.args
    assert len(df) == 2
    assert "knockdown" in df.columns
