@echo off
REM BVID-FE Windows local build script.
REM
REM Usage (from Windows cmd.exe or PowerShell, in the repo root):
REM   scripts\build_windows.bat
REM
REM Produces: dist\BVID-FE\BVID-FE.exe (a folder + exe — zip the folder to distribute).
REM Requires: Python 3.9+, pip, git. Will create a .venv if one is not present.

setlocal

REM Move to repo root (parent of scripts\)
cd /d "%~dp0\.."

echo === BVID-FE Windows build ===
echo.

REM 1. Ensure virtualenv exists
if not exist ".venv\Scripts\python.exe" (
    echo [1/4] Creating .venv\
    python -m venv .venv
    if errorlevel 1 (
        echo ERROR: python -m venv failed. Is Python 3.9+ on PATH?
        exit /b 1
    )
) else (
    echo [1/4] Reusing existing .venv\
)

REM 2. Install project + build deps
echo.
echo [2/4] Installing bvidfe[all] + pyinstaller (this may take a few minutes)...
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -e ".[all]"
.venv\Scripts\python.exe -m pip install pyinstaller
if errorlevel 1 (
    echo ERROR: pip install failed. See output above.
    exit /b 1
)

REM 3. Headless test-suite smoke check
echo.
echo [3/4] Running the test suite headlessly (QT_QPA_PLATFORM=offscreen)...
set QT_QPA_PLATFORM=offscreen
set BVIDFE_LOG_LEVEL=WARNING
.venv\Scripts\python.exe -m pytest tests\ -q --tb=short
if errorlevel 1 (
    echo ERROR: test suite failed. Aborting build.
    exit /b 1
)

REM 4. PyInstaller build
echo.
echo [4/4] Running PyInstaller (this will take 2-5 minutes)...
.venv\Scripts\python.exe -m PyInstaller BvidFE.spec --noconfirm --clean
if errorlevel 1 (
    echo ERROR: PyInstaller failed. See output above.
    exit /b 1
)

echo.
echo === BUILD COMPLETE ===
echo.
echo Output: dist\BVID-FE\BVID-FE.exe
echo.
echo To test locally: double-click dist\BVID-FE\BVID-FE.exe
echo To distribute:   zip the entire dist\BVID-FE\ folder and ship it.
echo.

endlocal
