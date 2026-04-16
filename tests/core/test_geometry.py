import pytest
from bvidfe.core.geometry import PanelGeometry, ImpactorGeometry


def test_panel_default_boundary_simply_supported():
    p = PanelGeometry(Lx_mm=150, Ly_mm=100)
    assert p.boundary == "simply_supported"


def test_panel_rejects_nonpositive_dims():
    with pytest.raises(ValueError):
        PanelGeometry(Lx_mm=0, Ly_mm=100)
    with pytest.raises(ValueError):
        PanelGeometry(Lx_mm=150, Ly_mm=-10)


def test_impactor_default_is_hemispherical_16mm():
    i = ImpactorGeometry()
    assert i.diameter_mm == 16.0
    assert i.shape == "hemispherical"


def test_impactor_rejects_bad_shape():
    with pytest.raises(ValueError):
        ImpactorGeometry(shape="triangular")  # type: ignore[arg-type]


def test_impactor_rejects_nonpositive_diameter():
    with pytest.raises(ValueError):
        ImpactorGeometry(diameter_mm=0)
