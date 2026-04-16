import numpy as np
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.core.laminate import Laminate


def test_symmetric_laminate_has_zero_B():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0, 45, -45, 90, 90, -45, 45, 0], ply_thickness_mm=0.152)
    A, B, D = lam.abd_matrices()
    assert np.allclose(B, 0.0, atol=1e-6)


def test_quasi_isotropic_A_is_isotropic():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0, 45, -45, 90] * 2, ply_thickness_mm=0.125)
    A, _, _ = lam.abd_matrices()
    # A11 == A22, A16 == A26 == 0 for quasi-iso
    assert abs(A[0, 0] - A[1, 1]) / A[0, 0] < 0.02
    assert abs(A[0, 2]) / A[0, 0] < 0.02
    assert abs(A[1, 2]) / A[0, 0] < 0.02


def test_effective_Ex_matches_ply_when_all_zero():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0] * 8, ply_thickness_mm=0.125)
    Ex, Ey, Gxy, nuxy = lam.effective_engineering_constants()
    assert abs(Ex - m.E11) / m.E11 < 0.01


def test_thickness_property():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0, 90, 0, 90], ply_thickness_mm=0.2)
    assert abs(lam.thickness_mm - 0.8) < 1e-9


def test_D_eff_positive():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(material=m, layup_deg=[0, 45, -45, 90] * 4, ply_thickness_mm=0.152)
    assert lam.flexural_rigidity_Deff() > 0
