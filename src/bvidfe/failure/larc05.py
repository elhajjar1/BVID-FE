"""LaRC05 (Hashin-3D reduction) composite failure criterion.

This is a minimal engineering implementation covering the four fundamental
LaRC05 modes (fiber tension, fiber compression, matrix tension, matrix
compression). For BVID CAI/TAI first-ply-failure prediction this is
sufficient; extended LaRC05 features (plane search for matrix cracking,
fiber kinking with non-linear shear) can be added in a future release.

Stress convention: Voigt 6-vector in the material frame,
    [sigma_11, sigma_22, sigma_33, tau_23, tau_13, tau_12]
with engineering shear strains.
"""

from __future__ import annotations

from typing import Sequence

from bvidfe.core.material import OrthotropicMaterial


def larc05_index(m: OrthotropicMaterial, stress: Sequence[float]) -> float:
    """Return max LaRC05 failure index across the four fundamental modes."""
    s1, s2, s3, t23, t13, t12 = stress

    modes: list[float] = []

    # Fiber tension (direction 1)
    if s1 >= 0:
        modes.append((s1 / m.Xt) ** 2)
    else:
        # Fiber compression (direction 1)
        modes.append((s1 / m.Xc) ** 2)

    # Matrix tension (direction 2)
    if s2 >= 0:
        modes.append((s2 / m.Yt) ** 2 + (t12 / m.S12) ** 2 + (t23 / m.S23) ** 2)
    else:
        # Matrix compression (direction 2)
        modes.append((s2 / m.Yc) ** 2 + (t12 / m.S12) ** 2 + (t23 / m.S23) ** 2)

    return max(modes)
