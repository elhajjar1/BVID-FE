#!/usr/bin/env bash
# BVID-FE macOS local build script.
#
# Usage (from the repo root):
#   bash scripts/build_macos.sh
#
# Produces: dist/BVID-FE.app
# Requires: Python 3.9+, pip. Will create .venv if not present.

set -euo pipefail

# Move to repo root (parent of scripts/)
cd "$(dirname "$0")/.."

echo "=== BVID-FE macOS build ==="
echo

# 1. venv
if [ ! -x ".venv/bin/python" ]; then
    echo "[1/4] Creating .venv/"
    python3 -m venv .venv
else
    echo "[1/4] Reusing existing .venv/"
fi

# 2. deps
echo
echo "[2/4] Installing bvidfe[all] + pyinstaller..."
./.venv/bin/pip install --upgrade pip >/dev/null
./.venv/bin/pip install -e ".[all]"
./.venv/bin/pip install pyinstaller

# 3. smoke test
echo
echo "[3/4] Running test suite headlessly..."
QT_QPA_PLATFORM=offscreen BVIDFE_LOG_LEVEL=WARNING \
    ./.venv/bin/python -m pytest tests/ -q --tb=short

# 4. build
echo
echo "[4/4] Running PyInstaller..."
./.venv/bin/python -m PyInstaller BvidFE.spec --noconfirm --clean

echo
echo "=== BUILD COMPLETE ==="
echo
echo "Output: dist/BVID-FE.app"
echo
echo "To test locally: open dist/BVID-FE.app"
echo "To distribute:   zip -r BVID-FE-macOS.zip dist/BVID-FE.app"
echo "                 Recipients may need: xattr -rd com.apple.quarantine BVID-FE.app"
