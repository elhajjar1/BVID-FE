"""Run the Python code blocks of the user-facing docs.

Every ```python fence in the files below is executed, in order, in one
namespace per file (later blocks may use names an earlier block defined),
from a temporary working directory that holds a copy of ``examples/`` so
relative paths such as ``examples/sample_cscan.json`` resolve. A stale
example (a renamed argument, a variable an earlier block rebinds) then
fails here instead of in a reader's terminal.

To leave a block out, put this marker on the line just before its fence,
with a reason:

    <!-- docs-test: skip (fragment) -->

Use it for fragments that do not run on their own and for blocks that need
the fe3d solver or Streamlit, which are too slow or not installed here.
"""

from __future__ import annotations

import re
import shutil
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

REPO = Path(__file__).resolve().parents[1]
DOC_FILES = ["README.md", "docs/quickstart.md", "docs/python_api.md"]

SKIP_MARKER = re.compile(r"<!--\s*docs-test:\s*skip\b.*-->")
FENCE_OPEN = re.compile(r"^```python\s*$")
FENCE_CLOSE = re.compile(r"^```\s*$")


def python_blocks(path: Path) -> list[tuple[int, str, bool]]:
    """Return ``(first code line, code, skipped)`` for each ```python fence."""
    lines = path.read_text(encoding="utf-8").splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        if FENCE_OPEN.match(lines[i]):
            skipped = i > 0 and bool(SKIP_MARKER.search(lines[i - 1]))
            start = i + 1
            j = start
            while j < len(lines) and not FENCE_CLOSE.match(lines[j]):
                j += 1
            if j == len(lines):
                raise AssertionError(f"{path.name}:{i + 1}: unclosed ```python fence")
            blocks.append((start + 1, "\n".join(lines[start:j]), skipped))
            i = j
        i += 1
    return blocks


@pytest.mark.parametrize("doc", DOC_FILES)
def test_doc_python_blocks_run(doc, tmp_path, monkeypatch):
    path = REPO / doc
    blocks = python_blocks(path)
    assert any(not skipped for _, _, skipped in blocks), f"no runnable python blocks in {doc}"

    shutil.copytree(REPO / "examples", tmp_path / "examples")
    monkeypatch.chdir(tmp_path)

    namespace: dict = {"__name__": "__docs_test__"}
    try:
        for first_line, code, skipped in blocks:
            if skipped:
                continue
            # Pad with blank lines so tracebacks point at the line in the doc.
            source = "\n" * (first_line - 1) + code
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                exec(compile(source, str(path), "exec"), namespace)
    finally:
        plt.close("all")


def test_block_parser_honours_skip_marker(tmp_path):
    doc = tmp_path / "doc.md"
    doc.write_text(
        "text\n```python\nx = 1\n```\n<!-- docs-test: skip (fragment) -->\n"
        "```python\nundefined_name\n```\n```bash\nls\n```\n",
        encoding="utf-8",
    )
    assert python_blocks(doc) == [(3, "x = 1", False), (7, "undefined_name", True)]
