import pytest

from bvidfe.core.material import MATERIAL_LIBRARY
from bvidfe.failure.tsai_wu import tsai_wu_index


def test_uniaxial_Xt_gives_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    s = [m.Xt, 0, 0, 0, 0, 0]
    assert abs(tsai_wu_index(m, s) - 1.0) < 0.05


def test_uniaxial_negative_Xc_gives_index_one():
    m = MATERIAL_LIBRARY["IM7/8552"]
    s = [-m.Xc, 0, 0, 0, 0, 0]
    assert abs(tsai_wu_index(m, s) - 1.0) < 0.05


def test_zero_stress_gives_zero_index():
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert tsai_wu_index(m, [0] * 6) == 0.0


def test_zt_default_falls_back_to_yt():
    """When Zt is not set, direction-3 tensile strength must equal direction-2
    (the unidirectional / transverse-isotropy assumption)."""
    m = MATERIAL_LIBRARY["IM7/8552"]
    assert m.Zt is None and m.Zc is None and m.S13 is None
    assert tsai_wu_index(m, [0, 0, m.Yt, 0, 0, 0]) == tsai_wu_index(m, [0, m.Yt, 0, 0, 0, 0])
    assert tsai_wu_index(m, [0, 0, -m.Yc, 0, 0, 0]) == tsai_wu_index(m, [0, -m.Yc, 0, 0, 0, 0])


def test_zt_override_changes_through_thickness_strength():
    """When Zt is explicitly set lower than Yt (typical for void-prone CFRP),
    the direction-3 tensile strength predicted by Tsai-Wu should track Zt and
    differ from the in-plane direction-2 strength."""
    from dataclasses import replace

    base = MATERIAL_LIBRARY["IM7/8552"]
    weaker = replace(base, Zt=0.5 * base.Yt, Zc=0.5 * base.Yc)
    # Uniaxial Tsai-Wu reaches 1 exactly at the strength in that direction.
    assert tsai_wu_index(weaker, [0, 0, weaker.Zt, 0, 0, 0]) == pytest.approx(1.0, rel=1e-12)
    assert tsai_wu_index(weaker, [0, weaker.Yt, 0, 0, 0, 0]) == pytest.approx(1.0, rel=1e-12)
    # Lower Zt -> the same stress is nearer failure through the thickness.
    s = 0.6 * weaker.Yt
    assert tsai_wu_index(weaker, [0, 0, s, 0, 0, 0]) > tsai_wu_index(weaker, [0, s, 0, 0, 0, 0])


def test_uniaxial_xz_shear_uses_s13():
    """Pure tau_xz loading (Voigt index 4) should produce a Tsai-Wu index
    governed by F55 = 1/S13^2 — previously F55 used S12, which was a real
    formula bug. With the unidirectional default S13_resolved == S12 the
    legacy behaviour is preserved; once S13 is overridden the index responds."""
    from dataclasses import replace

    base = MATERIAL_LIBRARY["IM7/8552"]
    s = [0.0, 0.0, 0.0, 0.0, base.S12, 0.0]  # tau_xz at the in-plane shear strength
    # Default behaviour: S13 falls back to S12, so index ≈ 1.0
    assert abs(tsai_wu_index(base, s) - 1.0) < 0.05
    # Override S13 to be twice as strong — same applied stress, lower index
    stiffer = replace(base, S13=2.0 * base.S12)
    assert tsai_wu_index(stiffer, s) < 0.5


def test_tsai_wu_index_matches_explicit_nested_sums():
    """Locks the vectorised tsai_wu_index against the explicit nested-sum
    formula F_i s_i + F_ij s_i s_j. The vectorised form uses
    F.dot(s) + s.dot(Q.dot(s)), which is mathematically identical but
    different floating-point evaluation order — this test catches any
    accidental regression to a wrong indexing or a sign flip."""
    import numpy as np

    from bvidfe.failure.tsai_wu import _tsai_wu_coefficients

    m = MATERIAL_LIBRARY["IM7/8552"]
    F, Q = _tsai_wu_coefficients(m)
    rng = np.random.default_rng(42)
    for _ in range(50):
        s = rng.standard_normal(6) * 100.0
        explicit = sum(F[i] * s[i] for i in range(6)) + sum(
            Q[i][j] * s[i] * s[j] for i in range(6) for j in range(6)
        )
        assert tsai_wu_index(m, s.tolist()) == pytest.approx(explicit, rel=1e-12, abs=1e-9)
