import matplotlib

matplotlib.use("Agg")

from bvidfe.core.geometry import PanelGeometry
from bvidfe.damage.state import DamageState, DelaminationEllipse
from bvidfe.viz.plots_2d import (
    plot_damage_map,
    plot_knockdown_curve,
    plot_tier_comparison,
)


def _sample_damage():
    return DamageState(
        [
            DelaminationEllipse(3, (75, 50), 28, 18, 45),
            DelaminationEllipse(4, (78, 52), 32, 20, 50),
            DelaminationEllipse(5, (80, 54), 36, 22, 55),
        ],
        dent_depth_mm=0.45,
    )


def test_plot_damage_map_saves_png(tmp_path):
    panel = PanelGeometry(150, 100)
    damage = _sample_damage()
    fig = plot_damage_map(damage, panel)
    out = tmp_path / "damage.png"
    fig.savefig(out, dpi=100)
    assert out.exists() and out.stat().st_size > 5000


def test_plot_knockdown_curve_saves_png(tmp_path):
    energies = [5, 10, 15, 20, 25, 30, 35]
    kd = [0.95, 0.85, 0.70, 0.55, 0.42, 0.32, 0.26]
    fig = plot_knockdown_curve(energies, kd, tier_label="empirical")
    out = tmp_path / "kd.png"
    fig.savefig(out, dpi=100)
    assert out.exists() and out.stat().st_size > 5000


def test_plot_tier_comparison_saves_png(tmp_path):
    energies = [5, 10, 15, 20, 25, 30, 35]
    results = {
        "empirical": [0.95, 0.85, 0.70, 0.55, 0.42, 0.32, 0.26],
        "semi_analytical": [0.96, 0.80, 0.60, 0.45, 0.35, 0.28, 0.24],
    }
    fig = plot_tier_comparison(energies, results)
    out = tmp_path / "cmp.png"
    fig.savefig(out, dpi=100)
    assert out.exists() and out.stat().st_size > 5000


def test_plot_empty_damage_handles_gracefully():
    panel = PanelGeometry(150, 100)
    damage = DamageState([], dent_depth_mm=0.0)
    fig = plot_damage_map(damage, panel)
    assert fig is not None
