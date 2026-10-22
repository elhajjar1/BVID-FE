"""CLI numeric snapshot tests.

The BVID-FE CLI emits a JSON result document, but no existing test pins
the exact numeric signature of canonical runs. A dependency bump (numpy,
scipy) or an accidental formula tweak can shift the ``knockdown`` by a
few percent and pass CI unnoticed. This module locks down the numeric
output of three canonical CLI invocations.

Tolerance
---------
The empirical and semi_analytical tiers are purely deterministic
closed-form / direct-solve evaluations: rerunning the same command
three times produced bit-identical floats. We assert with
``rel=1e-6`` so a tweak that shifts ``knockdown`` by even 0.0001%
trips the test, while still leaving headroom for benign cross-platform
floating-point drift in the lower bits.

The 3D FE tier (``fe3d``) is intentionally excluded — a single
invocation already exceeds the 5-second budget for the snapshot suite,
and the iterative eigensolve introduces low-bit drift that would
require a much looser tolerance (~1e-4) to be useful. Cross-tier
agreement is covered by ``tests/test_cross_tier_consistency.py``.

Regenerating snapshots
----------------------
When a formula change is intentional, re-capture the snapshots by
running each CLI invocation below with ``--quick-json`` and pasting
the new numeric fields into ``CANONICAL_RUNS``. Example::

    python -m bvidfe.cli --tier empirical --material IM7/8552 \\
        --layup 0,45,-45,90,90,-45,45,0,0,45,-45,90,90,-45,45,0 \\
        --thickness 0.152 --panel 200x150 --loading compression \\
        --energy 20 --quick-json

Always rerun ``pytest tests/test_cli_reproducibility.py`` after
updating, and document the formula change in the PR description.
"""

from __future__ import annotations

import json
import subprocess
import sys
from typing import Any

import pytest

# Each entry: (id, CLI args, expected JSON snapshot fields).
# Snapshots captured 2026-05-20 against bvidfe 0.2.0.dev0 with
# numpy 2.x and scipy 1.14.x on linux x86_64.
CANONICAL_RUNS: list[tuple[str, list[str], dict[str, Any]]] = [
    (
        "empirical_quasi_iso_quad_20J",
        [
            "--tier",
            "empirical",
            "--material",
            "IM7/8552",
            "--layup",
            "0,45,-45,90,90,-45,45,0,0,45,-45,90,90,-45,45,0",
            "--thickness",
            "0.152",
            "--panel",
            "200x150",
            "--loading",
            "compression",
            "--energy",
            "20",
        ],
        {
            # Recaptured when soutis_k_s went from 2.5 to 3.5: knockdown
            # 0.439 -> 0.359; and when DPA became the energy balance
            # (7824 -> 1855 mm^2) with soutis_k_s = 4.7 and the fixed
            # reference area: 0.377. Corrected Olsson threshold load (onset
            # 1.2 -> 14.3 J on this panel; DPA 1855 -> 628 mm^2) with
            # olsson_alpha 0.0865 and soutis_k_s 4.8: 0.505.
            "knockdown": 0.5045652386609865,
            "residual_strength_MPa": 303.38835202806507,
            "pristine_strength_MPa": 601.2866697539372,
            "dpa_mm2": 627.6927028441303,
            "tier_used": "empirical",
        },
    ),
    (
        "semi_analytical_quasi_iso_quad_22J",
        [
            "--tier",
            "semi_analytical",
            "--material",
            "IM7/8552",
            "--layup",
            "0,45,-45,90,90,-45,45,0,0,45,-45,90,90,-45,45,0",
            "--thickness",
            "0.152",
            "--panel",
            "200x150",
            "--loading",
            "compression",
            "--energy",
            "22",
        ],
        {
            # Recaptured when the sublaminate channel changed from buckling
            # onset to buckling-driven delamination growth — knockdown
            # rose from 0.055 to 0.306 — and again when the pristine
            # strength became ply-discount last-ply failure (887.5 -> 601.3
            # MPa): the growth stress is absolute, so the knockdown rose to
            # 0.451. With soutis_k_s = 3.5 the Soutis cap now binds:
            # 0.450. Energy-balance DPA (3660 -> 868 mm^2), soutis_k_s =
            # 4.7 and the fixed reference area instead of the 200x150 panel:
            # 0.469. With the corrected Olsson threshold load the onset is
            # 14.3 J on this panel and 10 J leaves no damage, so the case
            # moved to 22 J (DPA 847 mm^2, close to the 868 mm^2 at 10 J
            # before): 0.467.
            "knockdown": 0.46720340120844095,
            "residual_strength_MPa": 280.92317721033606,
            "pristine_strength_MPa": 601.2866697539372,
            "dpa_mm2": 846.6807392998476,
            "tier_used": "semi_analytical",
        },
    ),
    (
        "empirical_cross_ply_tension_15J",
        [
            "--tier",
            "empirical",
            "--material",
            "IM7/8552",
            "--layup",
            "0,90,0,90,90,0,90,0",
            "--thickness",
            "0.2",
            "--panel",
            "150x100",
            "--loading",
            "tension",
            "--energy",
            "15",
        ],
        {
            # Recaptured when IM7/8552 Xt became NCAMP's cross-ply back-out
            # value (2560 -> 2238 MPa): pristine 1280.6 -> 1119.6 MPa. Then
            # energy-balance DPA (8902 -> 1389 mm^2, below the 80% cap):
            # knockdown 0.219 -> 0.260. Corrected Olsson threshold load
            # (onset 0.9 -> 8.6 J; DPA 1389 -> 706 mm^2) and olsson_alpha 0.0865:
            # 0.290.
            "knockdown": 0.2897702040883179,
            "residual_strength_MPa": 324.4170461742307,
            "pristine_strength_MPa": 1119.566613810138,
            "dpa_mm2": 705.9357964219786,
            "tier_used": "empirical",
        },
    ),
]


# Fields snapshot-asserted with pytest.approx; the categorical
# ``tier_used`` is compared exactly.
NUMERIC_FIELDS = (
    "knockdown",
    "residual_strength_MPa",
    "pristine_strength_MPa",
    "dpa_mm2",
)


@pytest.mark.parametrize(
    ("case_id", "cli_args", "expected"),
    CANONICAL_RUNS,
    ids=[c[0] for c in CANONICAL_RUNS],
)
def test_cli_numeric_snapshot(
    case_id: str,
    cli_args: list[str],
    expected: dict[str, Any],
) -> None:
    """Assert the CLI's JSON output matches the captured numeric snapshot.

    Invoked via ``python -m bvidfe.cli`` (not the installed ``bvidfe``
    console script) so the test works in environments where the package
    is importable but the console entry point is not on PATH.
    """
    result = subprocess.run(
        [sys.executable, "-m", "bvidfe.cli", *cli_args],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"CLI failed for {case_id}: stderr={result.stderr!r}"

    data = json.loads(result.stdout)

    assert data["tier_used"] == expected["tier_used"], (
        f"{case_id}: tier_used mismatch (got {data['tier_used']!r}, "
        f"expected {expected['tier_used']!r})"
    )
    for field in NUMERIC_FIELDS:
        assert field in data, f"{case_id}: missing field {field!r} in CLI output"
        assert data[field] == pytest.approx(expected[field], rel=1e-6), (
            f"{case_id}: {field} drifted "
            f"(got {data[field]!r}, expected {expected[field]!r}). "
            "If this drift is intentional, re-capture the snapshot — see "
            "the module docstring for the regeneration recipe."
        )
