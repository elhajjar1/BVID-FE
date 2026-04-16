import numpy as np
import pytest

from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.evaluator import FailureEvaluator, LaminateFailureReport


def test_evaluator_tsai_wu_single_stress_at_Xt():
    m = MATERIAL_LIBRARY["IM7/8552"]
    ev = FailureEvaluator(m, criterion="tsai_wu")
    stress_field = np.array([[[m.Xt, 0, 0, 0, 0, 0]]])  # 1 elem, 1 gp, 6 components
    rpt = ev.evaluate(stress_field)
    assert isinstance(rpt, LaminateFailureReport)
    assert abs(rpt.max_index - 1.0) < 0.05
    assert rpt.critical_element == 0
    assert rpt.critical_gauss_point == 0


def test_evaluator_larc05_matches_expected_mode():
    m = MATERIAL_LIBRARY["IM7/8552"]
    ev = FailureEvaluator(m, criterion="larc05")
    stress_field = np.array([[[0, m.Yt, 0, 0, 0, 0]]])
    rpt = ev.evaluate(stress_field)
    assert abs(rpt.max_index - 1.0) < 1e-6


def test_evaluator_critical_element_picked():
    m = MATERIAL_LIBRARY["IM7/8552"]
    ev = FailureEvaluator(m, criterion="tsai_wu")
    stress_field = np.zeros((3, 2, 6))
    stress_field[1, 1, :] = [m.Xt, 0, 0, 0, 0, 0]  # element 1, gp 1 is critical
    rpt = ev.evaluate(stress_field)
    assert rpt.critical_element == 1
    assert rpt.critical_gauss_point == 1


def test_evaluator_unknown_criterion_raises():
    m = MATERIAL_LIBRARY["IM7/8552"]
    with pytest.raises(ValueError):
        FailureEvaluator(m, criterion="bogus")
