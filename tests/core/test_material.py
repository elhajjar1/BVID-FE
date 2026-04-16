import pytest
from bvidfe.core.material import OrthotropicMaterial, MATERIAL_LIBRARY


def test_material_library_has_four_presets():
    for name in ["AS4/3501-6", "IM7/8552", "T700/2510", "T800/epoxy"]:
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
