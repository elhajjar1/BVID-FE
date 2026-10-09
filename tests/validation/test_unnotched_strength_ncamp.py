"""Unnotched laminate strength against NCAMP measurements.

NCAMP qualification reports, RTD means normalized to the nominal cured ply
thickness, for the three standard layups:

* 25/50/25: [45/0/-45/90]2S (16 plies)
* 10/80/10: [45/-45/0/45/-45/90/45/-45/45/-45]S (20 plies)
* 50/40/10: [0/45/0/90/0/-45/0/45/0/-45]S (20 plies)

AS4/8552: CAM-RP-2010-002 Rev A (2011), Table 2-2, CPT 0.0074 in.
IM7/8552: CAM-RP-2009-015 Rev B (2019), sections 2.3.9-2.3.14, CPT 0.0072 in.

The IM7/8552 preset is not built from the NCAMP lamina data, so its
agreement is looser than AS4/8552's.

The old thickness-weighted ply-strength average was 1.17-2.85x these values
(mean 1.75x).
"""

from __future__ import annotations

import numpy as np
import pytest

from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.laminate_strength import unnotched_strength

_KSI = 6.894757  # MPa


def _sym(half):
    return half + half[::-1]


_LAYUPS = {
    "25/50/25": _sym([45, 0, -45, 90] * 2),
    "10/80/10": _sym([45, -45, 0, 45, -45, 90, 45, -45, 45, -45]),
    "50/40/10": _sym([0, 45, 0, 90, 0, -45, 0, 45, 0, -45]),
}
_CPT_MM = {"AS4/8552": 0.0074 * 25.4, "IM7/8552": 0.0072 * 25.4}
# (material, layup) -> (UNT, UNC) RTD normalized mean strength, ksi
_NCAMP_KSI = {
    ("AS4/8552", "25/50/25"): (88.61, 81.32),
    ("AS4/8552", "10/80/10"): (63.62, 62.42),
    ("AS4/8552", "50/40/10"): (152.32, 131.05),
    ("IM7/8552", "25/50/25"): (104.69, 87.05),
    ("IM7/8552", "10/80/10"): (67.01, 66.44),
    ("IM7/8552", "50/40/10"): (175.63, 120.84),
}
_CASES = [
    (mat, layup, loading, ksi * _KSI)
    for (mat, layup), pair in _NCAMP_KSI.items()
    for loading, ksi in zip(("tension", "compression"), pair)
]


def _predict(mat, layup, loading):
    lam = Laminate(MATERIAL_LIBRARY[mat], _LAYUPS[layup], _CPT_MM[mat])
    return unnotched_strength(lam, loading)


@pytest.mark.parametrize(
    ("mat", "layup", "loading", "measured_MPa"),
    _CASES,
    ids=[f"{m}-{lay}-{ld}" for m, lay, ld, _ in _CASES],
)
def test_unnotched_strength_within_20_percent_of_ncamp(mat, layup, loading, measured_MPa):
    ratio = _predict(mat, layup, loading) / measured_MPa
    assert 0.8 <= ratio <= 1.2, ratio


def test_unnotched_strength_is_unbiased_against_ncamp():
    ratios = np.array([_predict(m, lay, ld) / meas for m, lay, ld, meas in _CASES])
    assert 0.95 <= ratios.mean() <= 1.05, ratios
    assert np.exp(np.abs(np.log(ratios)).mean()) - 1 < 0.15, ratios
