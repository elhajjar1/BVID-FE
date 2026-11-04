import math

import numpy as np
import pytest
import scipy.linalg

from bvidfe.analysis.semi_analytical import (
    DegenerateThinSublaminateWarning,
    SemiAnalyticalResult,
    _clamped_ellipse_load_factor,
    _reduced_bending_stiffness,
    find_critical_interface,
    semi_analytical_cai,
    semi_analytical_tai,
    sublaminate_buckling_load,
)
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.damage.state import DamageState, DelaminationEllipse


def test_buckling_load_positive_for_realistic_case():
    m = MATERIAL_LIBRARY["IM7/8552"]
    layup = [0, 45, -45, 90] * 4  # 16 plies
    lam = Laminate(m, layup, 0.152)
    # Delamination at interface 3 (between ply 3 and ply 4)
    ell = DelaminationEllipse(3, (75, 50), major_mm=20, minor_mm=10, orientation_deg=0)
    N = sublaminate_buckling_load(lam, ell)
    assert N > 0


def test_buckling_load_decreases_with_ellipse_size():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ell_small = DelaminationEllipse(3, (0, 0), 10, 6, 0)
    ell_large = DelaminationEllipse(3, (0, 0), 40, 20, 0)
    N_small = sublaminate_buckling_load(lam, ell_small)
    N_large = sublaminate_buckling_load(lam, ell_large)
    assert N_small > N_large  # smaller ellipse => harder to buckle


def test_find_critical_interface_picks_largest_ellipse_near_surface():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState(
        [
            DelaminationEllipse(1, (0, 0), 10, 5, 0),  # near impact face
            DelaminationEllipse(7, (0, 0), 20, 10, 0),  # mid-plane, larger
            DelaminationEllipse(14, (0, 0), 8, 4, 0),  # near back face
        ],
        dent_depth_mm=0.3,
    )
    idx = find_critical_interface(ds, lam)
    # Interface 7 has largest area, but let's check scoring:
    # interface 7 z-upper = 7 * 0.152 = 1.064 (distance from top)
    # interface 7 z-lower = (16-8) * 0.152 = 1.216 (distance from bottom)
    # max |z| = 1.216
    # score_7 = pi*20*10 * 1.216 ≈ 764
    # score_1 = pi*10*5 * (16-2)*0.152 = pi*50 * 2.128 ≈ 334 (distance from bottom larger)
    # score_14 = pi*8*4 * (16-15)*0.152 = small
    # Expect 7 to win
    assert idx == 7


def test_find_critical_interface_returns_none_for_empty():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState([], dent_depth_mm=0.0)
    idx = find_critical_interface(ds, lam)
    assert idx is None


def test_sublaminate_above_interface_uses_correct_plies():
    """Critical sublaminate includes plies above the interface (from top surface to the delamination)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    layup = [0, 45, -45, 90, 90, -45, 45, 0]  # 8 plies
    lam = Laminate(m, layup, 0.152)
    # Delamination at interface 1 (between ply 1 and ply 2)
    # Upper sublaminate = plies 0, 1 (2 plies); thinner sublaminate buckles first.
    ell = DelaminationEllipse(1, (0, 0), 20, 10, 0)
    N = sublaminate_buckling_load(lam, ell)
    assert N > 0


def test_semi_analytical_cai_pristine_returns_pristine():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState([], dent_depth_mm=0.0)
    result = semi_analytical_cai(lam, ds, 500.0, 15000.0)
    assert isinstance(result, SemiAnalyticalResult)
    assert result.residual_strength_MPa == 500.0
    assert result.critical_interface_index is None
    assert result.critical_buckling_load_N is None


def test_semi_analytical_cai_returns_less_than_pristine():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState(
        [DelaminationEllipse(3, (0, 0), 30, 20, 0)],
        dent_depth_mm=0.5,
    )
    result = semi_analytical_cai(lam, ds, 500.0, 15000.0)
    assert isinstance(result, SemiAnalyticalResult)
    assert 0 < result.residual_strength_MPa < 500.0
    assert result.critical_interface_index == 3
    assert result.critical_buckling_load_N is not None
    assert result.critical_buckling_load_N > 0


def test_semi_analytical_cai_positive_for_large_delamination():
    """A large delamination still leaves a positive residual (min of Soutis
    and delamination growth; see test_sublaminate_growth.py)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds_large = DamageState(
        [DelaminationEllipse(3, (0, 0), 60, 40, 0)],
        dent_depth_mm=0.5,
    )
    result_large = semi_analytical_cai(lam, ds_large, 500.0, 15000.0)
    assert result_large.residual_strength_MPa > 0


def test_semi_analytical_result_is_frozen_dataclass():
    """``SemiAnalyticalResult`` is an immutable, named-attribute container —
    callers must not be able to mutate its fields after construction, and
    the structured access must replace positional-tuple unpacking."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState(
        [DelaminationEllipse(3, (0, 0), 30, 20, 0)],
        dent_depth_mm=0.5,
    )
    result = semi_analytical_cai(lam, ds, 500.0, 15000.0)
    # Frozen dataclass: assignment to a field raises.
    with pytest.raises((AttributeError, Exception)):
        result.residual_strength_MPa = 0.0  # type: ignore[misc]
    # Named-attribute access (the whole point of this dataclass).
    assert hasattr(result, "residual_strength_MPa")
    assert hasattr(result, "critical_interface_index")
    assert hasattr(result, "critical_buckling_load_N")


def test_semi_analytical_tai_pristine_returns_pristine():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds = DamageState([], dent_depth_mm=0.0)
    sigma_tai = semi_analytical_tai(lam, ds, 800.0)
    assert sigma_tai == 800.0


def test_semi_analytical_tai_monotonic_in_dpa():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    ds_small = DamageState([DelaminationEllipse(3, (0, 0), 10, 5, 0)], 0.3)
    ds_large = DamageState([DelaminationEllipse(3, (0, 0), 30, 20, 0)], 0.3)
    t1 = semi_analytical_tai(lam, ds_small, 800.0)
    t2 = semi_analytical_tai(lam, ds_large, 800.0)
    assert t1 > t2


def _reference_clamped_ellipse_N_cr(D, a, b, theta_deg, degree=14):
    """Independent Rayleigh-Ritz reference for a clamped ellipse, worked in
    the PANEL frame: semi-axes ``a`` (at ``theta_deg`` from x) and ``b``,
    the panel-frame bending stiffness ``D`` and uniaxial compression along
    x. Monomials ``x^p y^q`` (``p + q <= degree``) times the clamp factor
    ``g^2``, ``g = 1 - r^T Q r``; the kernel under test instead rotates the
    plies into the ellipse's frame and uses a Legendre basis there."""
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    rot = np.array([[c, -s], [s, c]])
    Q = rot @ np.diag([1.0 / a**2, 1.0 / b**2]) @ rot.T
    rg, rw = np.polynomial.legendre.leggauss(degree + 6)
    r = 0.5 * (rg + 1.0)
    n_t = 2 * degree + 12
    th = np.linspace(0.0, 2.0 * np.pi, n_t, endpoint=False)
    R, T = np.meshgrid(r, th, indexing="ij")
    xl, yl = a * R * np.cos(T), b * R * np.sin(T)
    x, y = (c * xl - s * yl).ravel(), (s * xl + c * yl).ravel()
    w = (np.outer(0.5 * rw, np.full(n_t, 2.0 * np.pi / n_t)) * R).ravel() * a * b
    g = 1.0 - (Q[0, 0] * x * x + 2.0 * Q[0, 1] * x * y + Q[1, 1] * y * y)
    gx, gy = -2.0 * (Q[0, 0] * x + Q[0, 1] * y), -2.0 * (Q[0, 1] * x + Q[1, 1] * y)
    G0, Gx, Gy = g * g, 2.0 * g * gx, 2.0 * g * gy
    Gxx = 2.0 * gx * gx - 4.0 * g * Q[0, 0]
    Gyy = 2.0 * gy * gy - 4.0 * g * Q[1, 1]
    Gxy = 2.0 * gx * gy - 4.0 * g * Q[0, 1]
    L = max(a, b)
    X, Y = x / L, y / L
    F = {k: [] for k in ("x", "y", "xx", "yy", "xy")}
    for p in range(degree + 1):
        for q in range(degree + 1 - p):
            m = X**p * Y**q
            mx = p * X ** max(p - 1, 0) * Y**q / L
            my = q * X**p * Y ** max(q - 1, 0) / L
            mxx = p * (p - 1) * X ** max(p - 2, 0) * Y**q / L**2
            myy = q * (q - 1) * X**p * Y ** max(q - 2, 0) / L**2
            mxy = p * q * X ** max(p - 1, 0) * Y ** max(q - 1, 0) / L**2
            F["x"].append(Gx * m + G0 * mx)
            F["y"].append(Gy * m + G0 * my)
            F["xx"].append(Gxx * m + 2.0 * Gx * mx + G0 * mxx)
            F["yy"].append(Gyy * m + 2.0 * Gy * my + G0 * myy)
            F["xy"].append(Gxy * m + Gx * my + Gy * mx + G0 * mxy)
    F = {k: np.array(v) for k, v in F.items()}

    def integral(f, h):
        return (F[f] * w) @ F[h].T

    K = (
        D[0, 0] * integral("xx", "xx")
        + D[1, 1] * integral("yy", "yy")
        + D[0, 1] * (integral("xx", "yy") + integral("yy", "xx"))
        + 4.0 * D[2, 2] * integral("xy", "xy")
        + 2.0 * D[0, 2] * (integral("xx", "xy") + integral("xy", "xx"))
        + 2.0 * D[1, 2] * (integral("yy", "xy") + integral("xy", "yy"))
    )
    return 1.0 / scipy.linalg.eigh(integral("x", "x"), K, eigvals_only=True)[-1]


# --- #132 PR0: the sublaminate is a clamped ellipse with D* = D - B A^-1 B ---


def test_buckling_load_is_the_clamped_ellipse_with_semi_axes():
    """Successor of the #29 regression tests (full rectangle sides, not
    semi-axes, in the old simply supported closed form): the sublaminate is
    now the clamped ellipse itself, with ``major_mm`` / ``minor_mm`` as its
    semi-axes. It matches an independent panel-frame Ritz on the coupled
    stiffness ``D*`` and scales as 1 / size^2."""
    mat = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(mat, [0.0, 90.0, 90.0, 0.0], 0.152)
    # Interface 1: [0, 90] over [90, 0]; the tie goes to the upper stack.
    D_star = _reduced_bending_stiffness(mat, [0.0, 90.0], 0.152)
    N = sublaminate_buckling_load(lam, DelaminationEllipse(1, (0.0, 0.0), 20.0, 20.0, 0.0))
    assert N == pytest.approx(_reference_clamped_ellipse_N_cr(D_star, 20.0, 20.0, 0.0), rel=1e-6)
    half = sublaminate_buckling_load(lam, DelaminationEllipse(1, (0.0, 0.0), 10.0, 10.0, 0.0))
    assert half == pytest.approx(4.0 * N, rel=1e-10)
    # The old simply supported 40 x 40 mm rectangle on the plain D was
    # 1.56x lower: the clamped edge raises the load about 2.8x and the
    # [0/90] stack's bending-extension coupling takes much of it back.
    D = Laminate(mat, [0.0, 90.0], 0.152).abd_matrices()[2]
    navier_40 = (math.pi / 40.0) ** 2 * (D[0, 0] + 2.0 * (D[0, 1] + 2.0 * D[2, 2]) + D[1, 1])
    assert N / navier_40 == pytest.approx(1.56, abs=0.01)


@pytest.mark.parametrize(
    "layup, iface, major, minor, film_MPa",
    [
        # Clamped-ellipse Ritz on D* from the #132 design probes (an
        # independent monomial-basis script; film stress N_cr / h, MPa).
        ([0, 90, 90, 0, 0, 90, 90, 0], 0, 12.0, 12.0, 27.0),
        ([0, 45, -45, 90, 90, -45, 45, 0], 0, 12.0, 12.0, 27.0),
        ([0, 45, -45, 90, 90, -45, 45, 0], 1, 12.0, 12.0, 49.4),
        ([0, 45, -45, 90, 90, -45, 45, 0], 2, 12.0, 12.0, 111.2),
        ([0, 45, -45, 90, 90, -45, 45, 0], 1, 12.0, 6.0, 137.1),
    ],
)
def test_buckling_load_matches_the_design_probe(layup, iface, major, minor, film_MPa):
    lam = Laminate(MATERIAL_LIBRARY["IM7/8552"], layup, 0.152)
    N = sublaminate_buckling_load(lam, DelaminationEllipse(iface, (0, 0), major, minor, 0.0))
    assert N / ((iface + 1) * 0.152) == pytest.approx(film_MPa, abs=0.06)


def test_bending_extension_coupling_lowers_the_sublaminate_buckling_load():
    """The [0/45] sublaminate over interface 1 of the 8-ply QI laminate is
    unsymmetric: ignoring B (the plain D about its mid-plane, as before
    #132) the clamped R = 12 mm circle buckles at 103.5 MPa, with
    ``D* = D - B A^-1 B`` at 49.4 MPa (design probe). D - D* is positive
    semi-definite, and D* = D for a symmetric stack."""
    mat = MATERIAL_LIBRARY["IM7/8552"]
    D = Laminate(mat, [0.0, 45.0], 0.152).abd_matrices()[2]
    D_star = _reduced_bending_stiffness(mat, [0.0, 45.0], 0.152)
    assert np.linalg.eigvalsh(D - D_star).min() > -1e-9 * D[0, 0]
    assert np.linalg.eigvalsh(D - D_star).max() > 0.1 * D[0, 0]
    plain = _clamped_ellipse_load_factor(D, 12.0, 12.0, 1.0, 0.0, 0.0) / 0.304
    coupled = _clamped_ellipse_load_factor(D_star, 12.0, 12.0, 1.0, 0.0, 0.0) / 0.304
    assert plain == pytest.approx(103.5, abs=0.1)
    assert coupled == pytest.approx(49.4, abs=0.1)
    lam = Laminate(mat, [0, 45, -45, 90, 90, -45, 45, 0], 0.152)
    N = sublaminate_buckling_load(lam, DelaminationEllipse(1, (0, 0), 12.0, 12.0, 0.0))
    assert N / 0.304 == pytest.approx(coupled, rel=1e-12)
    symmetric = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]
    D_sym = Laminate(mat, symmetric, 0.152).abd_matrices()[2]
    assert _reduced_bending_stiffness(mat, symmetric, 0.152) == pytest.approx(
        D_sym, rel=1e-9, abs=1e-9 * D_sym[0, 0]
    )


def test_reduced_bending_stiffness_does_not_depend_on_the_reference_plane():
    """Moving the reference plane by e turns (A, B, D) into
    (A, B + e A, D + 2 e B + e^2 A); D - B A^-1 B is unchanged."""
    mat = MATERIAL_LIBRARY["IM7/8552"]
    A, B, D = Laminate(mat, [0.0, 45.0, -45.0], 0.152).abd_matrices()
    e = 0.137
    B_e, D_e = B + e * A, D + 2.0 * e * B + e * e * A
    D_star = _reduced_bending_stiffness(mat, [0.0, 45.0, -45.0], 0.152)
    assert D_e - B_e @ np.linalg.solve(A, B_e) == pytest.approx(D_star, rel=1e-9, abs=1e-9)


def test_panel_boundary_does_not_scale_the_sublaminate():
    """The sublaminate's edges are the delamination front (#132): the
    panel's boundary factor (clamped 1.9, free 0.5) no longer applies."""
    lam = _lam8()
    e = DelaminationEllipse(2, (75, 50), 20.0, 8.0, 30.0)
    ss = sublaminate_buckling_load(lam, e, boundary="simply_supported")
    assert sublaminate_buckling_load(lam, e, boundary="clamped") == ss
    assert sublaminate_buckling_load(lam, e, boundary="free") == ss


def test_slender_ellipse_is_clipped_to_the_trusted_aspect_ratio():
    """Above 10:1 the shorter semi-axis is widened to a tenth of the longer
    one, with a DegenerateThinSublaminateWarning."""
    lam = _lam8()
    clipped = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 2.0, 30.0))
    with pytest.warns(DegenerateThinSublaminateWarning, match="aspect ratio"):
        slender = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 0.5, 30.0))
    assert slender == pytest.approx(clipped, rel=1e-12)
    with pytest.warns(DegenerateThinSublaminateWarning):
        across = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 0.5, 20.0, 30.0))
    assert across == pytest.approx(
        sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 2.0, 20.0, 30.0)), rel=1e-12
    )


# --- #165: the sublaminate follows the ellipse's orientation ---

_QI8 = [0, 45, -45, 90, 90, -45, 45, 0]


def _lam8():
    return Laminate(MATERIAL_LIBRARY["IM7/8552"], _QI8, 0.152)


def test_buckling_load_rotated_90_equals_swapped_axes():
    """An ellipse elongated across the load (orientation 90 deg) must buckle
    like the unrotated ellipse with its semi-axes swapped; before #165
    ``orientation_deg`` was ignored and both gave the along-load value."""
    lam = _lam8()
    across = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 90.0))
    swapped = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 8.0, 20.0, 0.0))
    along = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 0.0))
    assert across == pytest.approx(swapped, rel=1e-12)
    assert across != pytest.approx(along, rel=1e-3)
    # Fibre directions are axial: 180 deg is 0 deg and -90 deg is 90 deg.
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 180.0)
    ) == pytest.approx(along, rel=1e-12)
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, -90.0)
    ) == pytest.approx(across, rel=1e-12)


def test_buckling_load_follows_the_ellipse_not_its_bounding_box():
    """A tilted ellipse buckles as itself (#132), not as its axis-aligned
    bounding box (the #165 enclosing rectangle). Clamped plates obey domain
    monotonicity, so the 20 x 8 mm ellipse at any angle lies between the
    circles it contains (R = 8) and is contained in (R = 20)."""
    lam = _lam8()
    outer = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 20.0, 20.0, 0.0))
    inner = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 8.0, 8.0, 0.0))
    for theta in (0.0, 30.0, 45.0, -45.0, 60.0, 90.0):
        e = DelaminationEllipse(2, (75, 50), 20.0, 8.0, theta)
        tilted = sublaminate_buckling_load(lam, e)
        assert outer < tilted < inner, theta
        if theta % 90.0:
            ex, ey = e.half_extents_mm()
            box = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), ex, ey, 0.0))
            assert tilted != pytest.approx(box, rel=1e-2), theta


def test_buckling_load_of_a_circle_does_not_depend_on_its_orientation():
    """A circle's orientation is arbitrary; the plies are rotated into the
    ellipse's frame and the load out of it, so a sign error in either would
    show up here on an off-axis sublaminate."""
    lam = _lam8()
    ref = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), 12.0, 12.0, 0.0))
    for theta in (17.0, 45.0, 90.0, 133.0, -60.0):
        e = DelaminationEllipse(2, (75, 50), 12.0, 12.0, theta)
        assert sublaminate_buckling_load(lam, e) == pytest.approx(ref, rel=1e-10)


@pytest.mark.parametrize("major, minor, theta", [(20.0, 14.0, 30.0), (20.0, 8.0, -45.0)])
def test_rotated_buckling_load_matches_a_panel_frame_reference(major, minor, theta):
    """The ellipse-frame kernel equals an independent Ritz worked in the
    panel frame on the rotated ellipse, with the panel-frame D*."""
    lam = _lam8()
    D_star = _reduced_bending_stiffness(MATERIAL_LIBRARY["IM7/8552"], [0.0, 45.0, -45.0], 0.152)
    N = sublaminate_buckling_load(lam, DelaminationEllipse(2, (75, 50), major, minor, theta))
    assert N == pytest.approx(
        _reference_clamped_ellipse_N_cr(D_star, major, minor, theta), rel=1e-6
    )


def test_buckling_load_pins_issue_165_probe():
    """The 20 x 8 mm ellipse at interface 2 of the 8-ply QI laminate from
    #165: 77.7 N/mm along the load and 70.7 N/mm across it as a clamped
    ellipse on D* (#132). The #165 enclosing rectangle gave 78.4 and 39.5."""
    lam = _lam8()
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 0.0)
    ) == pytest.approx(77.71, abs=0.05)
    assert sublaminate_buckling_load(
        lam, DelaminationEllipse(2, (75, 50), 20.0, 8.0, 90.0)
    ) == pytest.approx(70.73, abs=0.05)


def test_tai_width_still_uses_the_extents():
    """The TAI damage width reads the rotated ellipse's half-extents
    (#165); the sublaminate buckling no longer does (#132)."""
    from bvidfe.failure.soutis_openhole import transverse_damage_width_mm

    for theta in (0.0, 30.0, 45.0, 90.0, -45.0):
        e = DelaminationEllipse(3, (75, 50), 18.0, 7.0, theta)
        _, ey = e.half_extents_mm()
        assert transverse_damage_width_mm(DamageState([e])) == pytest.approx(2.0 * ey, rel=1e-12)


def test_template_damage_growth_end_to_end_pins():
    """End-to-end regression pins on 150 x 100 mm IM7/8552 panels (16 mm,
    5.5 kg impactor), re-pinned for the clamped-ellipse sublaminate (#132
    PR0; the #165 values on the enclosing rectangle in brackets). 8-ply QI
    at 5 J: growth governs at 309.1 MPa, interface 5 (281.8, interface 4):
    knockdown 0.514 (0.469). At 10 J the film buckling stress at the
    governing interface 4 is 156.1 MPa (78.4) and Soutis governs (0.448
    both). On the [0/90]2s panel the growth stress rose above the Soutis
    bound at 5 and 10 J: 0.544 / 0.447 (0.424 / 0.389 with growth
    governing)."""
    from bvidfe import AnalysisConfig, BvidAnalysis
    from bvidfe.analysis.semi_analytical import weakest_sublaminate_growth
    from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
    from bvidfe.impact.mapping import ImpactEvent

    def run(layup, energy):
        cfg = AnalysisConfig(
            material="IM7/8552",
            layup_deg=layup,
            ply_thickness_mm=0.152,
            panel=PanelGeometry(150.0, 100.0),
            loading="compression",
            tier="semi_analytical",
            impact=ImpactEvent(energy, ImpactorGeometry(diameter_mm=16.0), mass_kg=5.5),
        )
        return BvidAnalysis(cfg).run()

    r5 = run(_QI8, 5.0)
    assert r5.knockdown == pytest.approx(0.514, abs=5e-4)
    assert r5.residual_strength_MPa == pytest.approx(309.1, abs=0.1)
    assert r5.critical_sublaminate == 5
    r10 = run(_QI8, 10.0)
    assert r10.knockdown == pytest.approx(0.448, abs=5e-4)
    g10 = weakest_sublaminate_growth(_lam8(), r10.damage)
    assert g10.interface_index == 4
    assert g10.film_buckling_stress_MPa == pytest.approx(156.1, abs=0.1)
    assert g10.growth_stress_MPa > r10.residual_strength_MPa  # Soutis governs
    xp = [0, 90, 0, 90, 90, 0, 90, 0]
    for energy, knockdown in ((5.0, 0.544), (10.0, 0.447)):
        r = run(xp, energy)
        assert r.knockdown == pytest.approx(knockdown, abs=5e-4)
        lam_xp = Laminate(MATERIAL_LIBRARY["IM7/8552"], xp, 0.152)
        growth = weakest_sublaminate_growth(lam_xp, r.damage)
        assert growth.growth_stress_MPa > r.residual_strength_MPa  # Soutis governs
