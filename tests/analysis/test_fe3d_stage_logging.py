"""fe3d per-stage progress log (issue #110).

Every real fe3d stage logs one INFO line at its start and one at its end
(mesh build, element build, assembly, BCs, solve, first-ply-failure
evaluation), with DOF counts and wall time, for both the damaged and the
undamaged run that ``BvidAnalysis`` makes."""

from __future__ import annotations

import re

import pytest

from bvidfe import BvidAnalysis
from tests.analysis._fe3d_fixtures import impact_cfg

_STAGES = [
    "mesh build start",
    "mesh build done",
    "element build start",
    "element build done",
    "assembly start",
    "assembly done",
    "BCs start",
    "BCs done",
    "solve start",
    "solve done",
    "FPF evaluation start",
    "FPF evaluation done",
    "done",
]
_TIMED = re.compile(r"\(\d+\.\d\ds\)$|\(total \d+\.\d\ds\)$")


def _stage(msg: str) -> str:
    """The stage name of a log line, without its run label."""
    head = re.sub(r"^fe3d (FPF|TAI) \((damaged, \d+ delaminations?|undamaged)\):? ?", "", msg)
    return re.sub(r" \(\d+\.\d\ds\)$", "", head).split(":")[0]


@pytest.mark.parametrize("loading, run", [("compression", "FPF"), ("tension", "TAI")])
def test_fe3d_logs_every_stage_of_the_damaged_and_undamaged_runs(caplog, loading, run):
    with caplog.at_level("INFO", logger="bvidfe.fe3d"):
        BvidAnalysis(impact_cfg(loading, 30.0)).run()
    msgs = [
        r.getMessage()
        for r in caplog.records
        if r.name.startswith("bvidfe.fe3d") and "buckling" not in r.getMessage()
    ]
    # One stage sequence per run, damaged first. The solver's fallback notes
    # ("ILU-GMRES did not reach rtol ...; using the direct solver") are not
    # stages.
    msgs = [m for m in msgs if not m.startswith("ILU")]
    assert [_stage(m) for m in msgs] == _STAGES * 2, msgs
    damaged, undamaged = msgs[: len(_STAGES)], msgs[len(_STAGES) :]
    assert re.match(rf"fe3d {run} \(damaged, \d+ delaminations?\): mesh build start$", damaged[0])
    assert undamaged[0] == f"fe3d {run} (undamaged): mesh build start"
    for lines in (damaged, undamaged):
        # Every "done" line carries its wall time; the DOF counts are on the
        # mesh, assembly, BC and solve lines.
        for stage, line in zip(_STAGES, lines):
            if stage.endswith("done"):
                assert _TIMED.search(line), line
        assert re.search(r"mesh build done: \d+ elements, \d+ DOFs", lines[1])
        assert re.search(r"assembly start: \d+ elements, \d+ DOFs", lines[4])
        assert re.search(r"BCs done: \d+ constrained DOFs eliminated, \d+ free", lines[7])
        assert re.search(r"solve start: method=\w+, \d+ free DOFs", lines[8])
        assert re.search(r"solve done: (direct LU|ILU-GMRES, \d+ iterations)", lines[9])
        criterion = "larc05" if run == "FPF" else "tsai_wu"
        assert lines[10].startswith(f"FPF evaluation start: {criterion}, ")
        assert re.search(r"FPF evaluation done: failure stress \d+\.\d MPa", lines[11])
        assert re.search(r"done: first-ply failure at \d+\.\d MPa \(total", lines[12])


def test_stage_log_is_silent_above_info(caplog):
    """At WARNING the stage lines are not emitted."""
    with caplog.at_level("WARNING", logger="bvidfe.fe3d"):
        BvidAnalysis(impact_cfg("compression", 20.0)).run()
    assert not [r for r in caplog.records if r.name.startswith("bvidfe.fe3d")]
