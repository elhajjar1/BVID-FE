import numpy as np

from bvidfe.elements.cohesive import CohesiveSurfaceElement


def test_cohesive_precracked_returns_zero_traction():
    elem = CohesiveSurfaceElement(
        sigma_n_max=60.0,
        tau_max=90.0,
        G_Ic=0.28,
        G_IIc=0.79,
        is_precracked=True,
    )
    T = elem.traction(np.array([0.1, 0.0, 0.0]))  # some separation
    assert np.allclose(T, 0.0)


def test_cohesive_pristine_elastic_regime():
    elem = CohesiveSurfaceElement(sigma_n_max=60.0, tau_max=90.0, G_Ic=0.28, G_IIc=0.79)
    # Very small separation => elastic, positive traction
    T = elem.traction(np.array([1e-6, 0.0, 0.0]))
    assert T[0] > 0
    # Purely compressive normal (negative delta_n) => elastic compressive traction (no damage)
    Tc = elem.traction(np.array([-1e-6, 0.0, 0.0]))
    assert Tc[0] < 0


def test_cohesive_mode_I_reaches_peak_and_softens():
    elem = CohesiveSurfaceElement(sigma_n_max=60.0, tau_max=90.0, G_Ic=0.28, G_IIc=0.79)
    # Sweep delta_n and find peak
    deltas = np.linspace(0, 0.05, 201)
    tractions = [elem.traction(np.array([d, 0, 0]))[0] for d in deltas]
    peak = max(tractions)
    peak_i = tractions.index(peak)
    assert peak > 0
    # After peak, traction decreases
    assert tractions[peak_i + 5] < peak


def test_cohesive_mode_I_energy_integrates_to_GIc():
    elem = CohesiveSurfaceElement(sigma_n_max=60.0, tau_max=90.0, G_Ic=0.28, G_IIc=0.79)
    # Integrate T dot d(delta) over opening from 0 to delta_f
    delta_f = 2 * 0.28 / 60.0
    deltas = np.linspace(0, delta_f + 1e-4, 2001)
    tractions = [elem.traction(np.array([d, 0, 0]))[0] for d in deltas]
    # Trapezoidal integration: integral T d(delta) = G_Ic
    area = np.trapezoid(tractions, deltas)
    assert abs(area - 0.28) / 0.28 < 0.02  # within 2%


def test_cohesive_fully_open_zero_traction():
    elem = CohesiveSurfaceElement(sigma_n_max=60.0, tau_max=90.0, G_Ic=0.28, G_IIc=0.79)
    delta_f = 2 * 0.28 / 60.0
    T = elem.traction(np.array([delta_f * 1.5, 0, 0]))
    assert np.allclose(T, 0.0, atol=1e-6)


def test_cohesive_stiffness_matrix_shape():
    elem = CohesiveSurfaceElement(sigma_n_max=60.0, tau_max=90.0, G_Ic=0.28, G_IIc=0.79)
    # Nodal displacements: (8 nodes, 3 DOFs each) = 24
    K = elem.stiffness_matrix_point()
    assert K.shape == (3, 3)  # per Gauss point constitutive tangent
