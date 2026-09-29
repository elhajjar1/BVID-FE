"""Ply-angle convention and material-frame stress recovery for Hex8Element.

Failure criteria take stress in the ply material frame
``[s1, s2, s3, t23, t13, t12]`` (1 = fiber). Two things must therefore agree:

* the rotated stiffness must be that of a ply whose fibers sit at
  +ply_angle_deg from x, the same convention CLT uses
  (``laminate._transform_reduced_stiffness``); and
* stress recovered for the failure check must be rotated back into that
  ply's material frame.

Before v0.2.1 the element built ``T C T^T`` with ``T`` the global-to-ply
stress rotation, which is the stiffness of a -theta ply (the xx-xy coupling
had the opposite sign to CLT), and fe3d handed global-frame stress to the
criteria, so a 90 deg ply's transverse stress was checked against Xc.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from bvidfe.core.laminate import _transform_reduced_stiffness
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.elements.hex8 import Hex8Element

_IN_PLANE = [0, 1, 5]  # Voigt xx, yy, xy
_ANGLES = [-60.0, -45.0, -30.0, 0.0, 15.0, 30.0, 45.0, 90.0]


def _unit_cube_nodes():
    return np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )


def _plane_stress_reduced(C: np.ndarray) -> np.ndarray:
    """3x3 plane-stress reduced stiffness (xx, yy, xy) from a 6x6 stiffness."""
    S = np.linalg.inv(C)
    return np.linalg.inv(S[np.ix_(_IN_PLANE, _IN_PLANE)])


def _uniform_strain_displacements(nodes: np.ndarray, eps_voigt: np.ndarray) -> np.ndarray:
    """Nodal displacements (24,) for a uniform engineering strain field."""
    exx, eyy, ezz, gyz, gxz, gxy = eps_voigt
    E = np.array(
        [
            [exx, gxy / 2, gxz / 2],
            [gxy / 2, eyy, gyz / 2],
            [gxz / 2, gyz / 2, ezz],
        ]
    )
    return (nodes @ E.T).reshape(-1)


@pytest.mark.parametrize("angle", _ANGLES)
def test_rotated_stiffness_matches_clt_convention(angle):
    """The element's in-plane plane-stress stiffness equals CLT's Q-bar for
    the same ply angle, including the sign of the Q16 / Q26 couplings."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    elem = Hex8Element(_unit_cube_nodes(), m, ply_angle_deg=angle)
    Q = _plane_stress_reduced(m.get_stiffness_matrix())
    expected = _transform_reduced_stiffness(Q, math.radians(angle))
    np.testing.assert_allclose(
        _plane_stress_reduced(elem._C_global),
        expected,
        rtol=1e-9,
        atol=1e-9 * np.abs(expected).max(),
    )


@pytest.mark.parametrize("angle", _ANGLES)
def test_uniaxial_stress_along_fiber_reads_as_pure_sigma1(angle):
    """A uniform uniaxial stress along the fiber direction (cos t, sin t)
    must come back from the material-frame recovery as pure s1."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    nodes = _unit_cube_nodes()
    elem = Hex8Element(nodes, m, ply_angle_deg=angle)
    c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    sigma_global = 100.0 * np.array([c * c, s * s, 0.0, 0.0, 0.0, c * s])
    eps_global = np.linalg.solve(elem._C_global, sigma_global)
    u = _uniform_strain_displacements(nodes, eps_global)

    np.testing.assert_allclose(
        elem.stress_at_gauss_points(u), np.tile(sigma_global, (8, 1)), atol=1e-8
    )
    expected_material = np.tile([100.0, 0.0, 0.0, 0.0, 0.0, 0.0], (8, 1))
    np.testing.assert_allclose(
        elem.stress_at_gauss_points_material(u), expected_material, atol=1e-8
    )


def test_ninety_degree_ply_maps_global_xx_to_transverse_sigma2():
    m = MATERIAL_LIBRARY["IM7/8552"]
    nodes = _unit_cube_nodes()
    elem = Hex8Element(nodes, m, ply_angle_deg=90.0)
    sigma_global = np.array([-50.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    u = _uniform_strain_displacements(nodes, np.linalg.solve(elem._C_global, sigma_global))
    sigma_material = elem.stress_at_gauss_points_material(u)
    np.testing.assert_allclose(sigma_material[:, 1], -50.0, atol=1e-8)  # s2
    np.testing.assert_allclose(sigma_material[:, 0], 0.0, atol=1e-8)  # s1


def test_zero_degree_material_stress_equals_global_stress():
    m = MATERIAL_LIBRARY["IM7/8552"]
    nodes = _unit_cube_nodes()
    elem = Hex8Element(nodes, m, ply_angle_deg=0.0)
    rng = np.random.default_rng(3)
    u = rng.standard_normal(24) * 1e-3
    np.testing.assert_array_equal(
        elem.stress_at_gauss_points_material(u), elem.stress_at_gauss_points(u)
    )
