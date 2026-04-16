"""Failure-criterion evaluator across an element x Gauss-point stress field."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from bvidfe.core.material import OrthotropicMaterial
from bvidfe.failure.larc05 import larc05_index
from bvidfe.failure.tsai_wu import tsai_wu_index

CriterionName = Literal["tsai_wu", "larc05"]


@dataclass
class LaminateFailureReport:
    """Summary of a failure evaluation across a laminate stress field."""

    max_index: float
    critical_element: int
    critical_gauss_point: int
    criterion: str


class FailureEvaluator:
    """Applies a failure criterion across every (element, gauss-point) stress."""

    def __init__(self, material: OrthotropicMaterial, criterion: CriterionName = "tsai_wu"):
        if criterion not in ("tsai_wu", "larc05"):
            raise ValueError(f"unknown criterion {criterion!r}")
        self.material = material
        self.criterion = criterion

    def _index(self, stress):
        if self.criterion == "tsai_wu":
            return tsai_wu_index(self.material, stress)
        return larc05_index(self.material, stress)

    def evaluate(self, stress_field: np.ndarray) -> LaminateFailureReport:
        """stress_field shape (n_elem, n_gp, 6) in the material frame."""
        assert stress_field.ndim == 3 and stress_field.shape[2] == 6
        max_idx = -1.0
        crit_e = 0
        crit_g = 0
        for e in range(stress_field.shape[0]):
            for g in range(stress_field.shape[1]):
                idx = self._index(stress_field[e, g])
                if idx > max_idx:
                    max_idx = idx
                    crit_e = e
                    crit_g = g
        return LaminateFailureReport(
            max_index=max_idx,
            critical_element=crit_e,
            critical_gauss_point=crit_g,
            criterion=self.criterion,
        )
