#!/usr/bin/env bash
# SessionStart hook for Claude Code.
#
# Reports any black formatting drift and ruff lint issues so the main Claude
# agent sees them in its initial context and can fix them before delegating
# to subagents (which historically commit without running formatters).
#
# This hook is informational only: it always exits 0 and never blocks the
# session.

set -u
# Intentionally not using -e: we want to keep running even if one of the
# checks reports failures, and we always want to exit 0.

REPO_ROOT="${CLAUDE_PROJECT_DIR:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
cd "$REPO_ROOT" || exit 0

# Only check directories that exist (keeps the hook safe on partial clones).
TARGETS=()
[ -d src ] && TARGETS+=("src")
[ -d tests ] && TARGETS+=("tests")

if [ "${#TARGETS[@]}" -eq 0 ]; then
    echo "[session-start] No src/ or tests/ directory found; skipping formatter checks."
    exit 0
fi

echo "[session-start] Running formatter/linter drift checks on: ${TARGETS[*]}"
echo

# ---- black --------------------------------------------------------------
if command -v black >/dev/null 2>&1; then
    black --check --quiet "${TARGETS[@]}" >/dev/null 2>&1
    BLACK_RC=$?
    if [ $BLACK_RC -eq 0 ]; then
        echo "[session-start] black: clean (no formatting drift)."
    else
        # Re-run without --quiet to capture the "would reformat" lines for the agent.
        DRIFT="$(black --check "${TARGETS[@]}" 2>&1 | grep -E '^would reformat' | sed 's/^would reformat //')"
        echo "[session-start] black: FORMATTING DRIFT DETECTED"
        if [ -n "$DRIFT" ]; then
            echo "  Files needing reformatting:"
            while IFS= read -r f; do
                [ -n "$f" ] && echo "    - $f"
            done <<< "$DRIFT"
        fi
        echo "  Action: run \`black src tests\` before committing."
        echo "  (CI runs \`black --check src tests\` and will fail until this is fixed.)"
    fi
else
    echo "[session-start] black: not installed (run \`pip install -e .[dev]\` to enable this check)."
fi

echo

# ---- ruff ---------------------------------------------------------------
if command -v ruff >/dev/null 2>&1; then
    RUFF_OUT="$(ruff check "${TARGETS[@]}" 2>&1)"
    RUFF_RC=$?
    if [ $RUFF_RC -eq 0 ]; then
        echo "[session-start] ruff: clean (no lint issues)."
    else
        echo "[session-start] ruff: LINT ISSUES DETECTED"
        # Indent ruff output for readability in the agent context.
        echo "$RUFF_OUT" | sed 's/^/  /'
        echo "  Action: run \`ruff check --fix src tests\` (or fix manually) before committing."
    fi
else
    echo "[session-start] ruff: not installed (run \`pip install -e .[dev]\` to enable this check)."
fi

# Always succeed: this hook is informational, not a hard gate.
exit 0
