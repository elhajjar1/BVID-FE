import dataclasses
import math

from bvidfe.core.geometry import ImpactorGeometry, PanelGeometry
from bvidfe.core.laminate import Laminate
from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.impact.olsson import NAVIER_N, onset_energy, threshold_load


def test_threshold_load_scales_with_sqrt_G_IIc():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    Pc1 = threshold_load(lam, pan, imp)
    m4 = dataclasses.replace(m, G_IIc=m.G_IIc * 4)
    lam2 = Laminate(m4, [0, 45, -45, 90] * 4, 0.152)
    Pc2 = threshold_load(lam2, pan, imp)
    # Pc ∝ sqrt(G_IIc), so 4x G_IIc => 2x Pc
    assert math.isclose(Pc2 / Pc1, 2.0, rel_tol=0.05)


def test_onset_energy_positive_and_monotonic_in_thickness():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam_thin = Laminate(m, [0, 90] * 4, 0.125)  # 8 plies
    lam_thick = Laminate(m, [0, 90] * 12, 0.125)  # 24 plies
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    E_thin = onset_energy(lam_thin, pan, imp)
    E_thick = onset_energy(lam_thick, pan, imp)
    assert E_thin > 0 and E_thick > 0
    assert E_thick > E_thin  # thicker plate harder to damage


def test_navier_n_is_11_by_default():
    assert NAVIER_N == 11


def test_onset_energy_scales_with_material_G_IIc():
    m = MATERIAL_LIBRARY["IM7/8552"]
    lam1 = Laminate(m, [0, 45, -45, 90] * 4, 0.152)
    m4 = dataclasses.replace(m, G_IIc=m.G_IIc * 4)
    lam2 = Laminate(m4, [0, 45, -45, 90] * 4, 0.152)
    pan = PanelGeometry(150, 100)
    imp = ImpactorGeometry()
    E1 = onset_energy(lam1, pan, imp)
    E2 = onset_energy(lam2, pan, imp)
    # onset energy ∝ Pc^2 ∝ G_IIc, so 4x => 4x (approximately)
    assert E2 > E1
