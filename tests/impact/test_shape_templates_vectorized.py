"""Regression tests for the vectorized ``_build`` and interpolant-cached
Brent root-find in ``bvidfe.impact.shape_templates`` (issue #87).

Two guarantees:
1. The vectorized ``_template_arrays`` reproduces the scalar reference
   implementation exactly (``ar``, ``rel``, ``orient``) for several
   layups.
2. ``distribute_damage`` lands on a known-good total DPA pinned from
   the pre-refactor implementation, to within 1e-6 relative.
"""

from __future__ import annotations

import math

import numpy as np

from bvidfe.damage.state import DamageState
from bvidfe.impact.shape_templates import (
    _aspect_ratio,
    _orientation_deg,
    _relative_size,
    _template_arrays,
    distribute_damage,
)


def _reference_arrays(layup, n_interfaces):
    """Pure-Python reference for (ar, rel, orient) per interface."""
    ar = [_aspect_ratio(layup[i + 1] - layup[i]) for i in range(n_interfaces)]
    orient = [_orientation_deg(layup[i], layup[i + 1]) for i in range(n_interfaces)]
    rel = [_relative_size(i, n_interfaces) for i in range(n_interfaces)]
    return ar, rel, orient


def test_template_arrays_matches_scalar_reference():
    layups = [
        [0, 45, -45, 90, 0, 90, -45, 45, 0],  # 9 plies
        [0, 90, 0],  # 3 plies
        [0, 30, -30, 60, -60, 90, 0],  # 7 plies
        [0, 0, 0, 0, 0, 0],  # 6 plies, all aligned
        [45, -45, 45, -45, 45, -45, 45],  # 7 plies, alternating
    ]
    for layup in layups:
        n_interfaces = len(layup) - 1
        ar_arr, rel_arr, orient_arr = _template_arrays(layup, n_interfaces)
        ar_ref, rel_ref, orient_ref = _reference_arrays(layup, n_interfaces)
        np.testing.assert_allclose(ar_arr, ar_ref, rtol=0, atol=0)
        np.testing.assert_allclose(rel_arr, rel_ref, rtol=0, atol=0)
        np.testing.assert_allclose(orient_arr, orient_ref, rtol=0, atol=0)


def test_vectorized_build_semi_axes_match_scalar_recipe():
    """For several scalars, vectorized major/minor axes equal
    scalar * AR * rel and scalar * rel respectively."""
    layup = [0, 45, -45, 90, 90, -45, 45, 0]
    n_interfaces = len(layup) - 1
    ar_arr, rel_arr, orient_arr = _template_arrays(layup, n_interfaces)
    for scalar in [0.5, 1.7, 5.0, 12.3]:
        majors = scalar * ar_arr * rel_arr
        minors = scalar * rel_arr
        # Compare against the Python recipe (pre-refactor formula).
        for i in range(n_interfaces):
            ar_i = _aspect_ratio(layup[i + 1] - layup[i])
            rel_i = _relative_size(i, n_interfaces)
            assert math.isclose(majors[i], scalar * ar_i * rel_i, rel_tol=0, abs_tol=0)
            assert math.isclose(minors[i], scalar * rel_i, rel_tol=0, abs_tol=0)


def test_distribute_damage_pinned_dpa_known_good():
    """Pin the total projected DPA. The Brent solver has ``xtol=1e-3`` on
    the scalar, so the DPA itself lands within ~1e-4 relative of the
    target. A drift larger than 1e-3 rel from the pinned DPA (or larger
    than 1% from the target) signals that the optimization compromised
    numerical equivalence."""
    layup = [0, 45, -45, 90, 90, -45, 45, 0]
    target = 800.0
    ellipses = distribute_damage(layup, target, 0.3, 0.0)
    ds = DamageState(ellipses, 0.3, 0.0)
    # Captured with the axial (mod-180) template, no 90/90 delamination.
    known_good = 800.0177263133414
    assert (
        abs(ds.projected_damage_area_mm2 - known_good) / known_good < 1e-3
    ), f"DPA drift: got {ds.projected_damage_area_mm2!r}, expected ~{known_good!r}"
    # And of course still meet the 1% contract against the target.
    assert abs(ds.projected_damage_area_mm2 - target) / target < 0.01


def test_distribute_damage_pinned_ellipse_axes():
    """Pin the per-interface semi-axes and orientations.

    Tolerance is ``rel_tol=1e-3`` because the Brent root-find sets
    ``xtol=1e-3`` on the scalar multiplier. The 90/90 mid-plane interface
    (index 3) carries no delamination; each orientation is the fibre angle
    of the back-face-side ply; -45/90 is a 45 deg mismatch (AR 2.125).
    """
    layup = [0, 45, -45, 90, 90, -45, 45, 0]
    ellipses = distribute_damage(layup, 800.0, 0.3, 0.0)
    ellipses.sort(key=lambda e: e.interface_index)
    expected = [
        (0, 7.0790573754156325, 3.3313211178426503, 45.0),
        (1, 13.533492041235768, 4.164151397303313, -45.0),
        (2, 10.618586063123448, 4.996981676763975, 90.0),
        (4, 14.158114750831267, 6.6626422356853015, -45.0),
        (5, 24.360285674224382, 7.495472515145964, 45.0),
        (6, 17.69764343853908, 8.328302794606627, 0.0),
    ]
    assert len(ellipses) == len(expected)
    for e, (idx, maj, minr, orient) in zip(ellipses, expected):
        assert e.interface_index == idx
        assert math.isclose(e.major_mm, maj, rel_tol=1e-3)
        assert math.isclose(e.minor_mm, minr, rel_tol=1e-3)
        assert math.isclose(e.orientation_deg, orient, rel_tol=0, abs_tol=0)
