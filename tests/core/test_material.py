import pytest
from bvidfe.core.material import OrthotropicMaterial, MATERIAL_LIBRARY


def test_material_library_has_five_presets():
    for name in ["AS4/3501-6", "AS4/8552", "IM7/8552", "T700/2510", "T800/epoxy"]:
        assert name in MATERIAL_LIBRARY
        m = MATERIAL_LIBRARY[name]
        assert isinstance(m, OrthotropicMaterial)
        assert m.E11 > m.E22 > 0


def test_orthotropic_material_stiffness_matrix_is_symmetric():
    m = MATERIAL_LIBRARY["IM7/8552"]
    C = m.get_stiffness_matrix()  # 6x6
    assert C.shape == (6, 6)
    import numpy as np

    assert np.allclose(C, C.T, atol=1e-6)


def test_material_rejects_negative_modulus():
    with pytest.raises(ValueError):
        OrthotropicMaterial(
            name="bad",
            E11=-1.0,
            E22=10.0,
            nu12=0.3,
            G12=5.0,
            G13=5.0,
            G23=3.0,
            Xt=1.0,
            Xc=1.0,
            Yt=1.0,
            Yc=1.0,
            S12=1.0,
            S23=1.0,
            G_Ic=0.1,
            G_IIc=0.5,
            rho=1.6e-6,
        )


@pytest.mark.parametrize("name", sorted(MATERIAL_LIBRARY))
def test_preset_name_matches_key_and_compliance_inverts_stiffness(name):
    import numpy as np

    m = MATERIAL_LIBRARY[name]
    assert m.name == name
    S, C = m.get_compliance_matrix(), m.get_stiffness_matrix()
    assert S.shape == C.shape == (6, 6)
    np.testing.assert_allclose(S @ C, np.eye(6), atol=1e-9)


_KSI = 6.894757  # MPa
_MSI = 6894.757  # MPa


def test_as4_8552_lamina_properties_are_ncamp_rtd_values():
    """NCAMP CAM-RP-2010-002 Rev A, Table 2-1, RTD means normalized to
    CPT = 0.0074 in where NCAMP normalizes (fiber-dominated properties)."""
    m = MATERIAL_LIBRARY["AS4/8552"]
    assert m.E11 == pytest.approx(18.46 * _MSI, rel=1e-3)
    assert m.E22 == pytest.approx(1.34 * _MSI, rel=1e-3)
    assert m.nu12 == pytest.approx(0.302, abs=0.005)
    assert m.G12 == pytest.approx(0.70 * _MSI, rel=1e-3)
    assert m.Xt == pytest.approx(289.47 * _KSI, rel=1e-3)
    assert m.Xc == pytest.approx(202.70 * _KSI, rel=1e-3)
    assert m.Yt == pytest.approx(9.27 * _KSI, rel=1e-2)
    assert m.Yc == pytest.approx(38.85 * _KSI, rel=1e-3)
    assert m.S12 == pytest.approx(13.28 * _KSI, rel=1e-2)
    # Laminate density 1.58-1.60 g/cc across the NCAMP panels.
    assert 1.58e-6 <= m.rho <= 1.60e-6


def test_as4_8552_derived_constants():
    import math

    m = MATERIAL_LIBRARY["AS4/8552"]
    # Transverse isotropy with nu23 = 0.45 (Gonzalez et al., AS4/8552).
    assert m.G13 == m.G12
    assert m.G23 == pytest.approx(m.E22 / (2 * (1 + 0.45)), rel=1e-3)
    # LaRC03 transverse shear strength from Yc at a 53 deg fracture plane.
    a0 = math.radians(53.0)
    S_T = m.Yc * math.cos(a0) * (math.sin(a0) + math.cos(a0) / math.tan(2 * a0))
    assert m.S23 == pytest.approx(S_T, rel=1e-2)
    # Same 8552 matrix: the IM7/8552 toughnesses, as the Girona AS4/8552
    # CAI simulations use.
    assert (m.G_Ic, m.G_IIc) == (0.28, 0.79)
