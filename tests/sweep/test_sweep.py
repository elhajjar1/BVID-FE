import pandas as pd

from bvidfe.analysis import AnalysisConfig
from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.impact.mapping import ImpactEvent
from bvidfe.sweep.parametric_sweep import (
    sweep_energies,
    sweep_layups,
    sweep_thicknesses,
)


def _base_impact_cfg(**overrides):
    kw = dict(
        material="IM7/8552",
        layup_deg=[0, 45, -45, 90] * 4,
        ply_thickness_mm=0.152,
        panel=PanelGeometry(150, 100),
        loading="compression",
        tier="empirical",
        impact=ImpactEvent(10.0, ImpactorGeometry(), mass_kg=5.5),
    )
    kw.update(overrides)
    return AnalysisConfig(**kw)


def test_sweep_energies_returns_dataframe_with_expected_columns():
    cfg = _base_impact_cfg()
    df = sweep_energies(cfg, energies_J=[5, 10, 20, 30])
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    for col in ["energy_J", "knockdown", "residual_MPa", "dpa_mm2", "dent_mm"]:
        assert col in df.columns


def test_sweep_energies_writes_csv(tmp_path):
    cfg = _base_impact_cfg()
    csv_path = tmp_path / "sweep.csv"
    sweep_energies(cfg, energies_J=[5, 10], csv_path=csv_path)
    assert csv_path.exists()
    loaded = pd.read_csv(csv_path)
    assert len(loaded) == 2


def test_sweep_layups():
    cfg = _base_impact_cfg()
    layups = [
        [0, 45, -45, 90] * 4,
        [0, 90] * 8,
        [0, 60, -60] * 5 + [0],
    ]
    df = sweep_layups(cfg, layups=layups)
    assert len(df) == 3
    assert "layup" in df.columns
    assert "knockdown" in df.columns


def test_sweep_thicknesses():
    cfg = _base_impact_cfg()
    df = sweep_thicknesses(cfg, ply_thicknesses_mm=[0.125, 0.152, 0.2])
    assert len(df) == 3
    assert "ply_thickness_mm" in df.columns
