@echo off
setlocal
cd /d "%~dp0"
title PlasmidGuard

rem --- find Python (the "py" launcher first; plain "python" can be a Microsoft Store stub) ---
set "PY="
py -3 --version >nul 2>nul && set "PY=py -3"
if not defined PY (
  python --version >nul 2>nul && set "PY=python"
)
if not defined PY (
  echo.
  echo Python was not found on this computer.
  echo Install Python 3.10 or newer from https://www.python.org/downloads/
  echo During install, tick "Add python.exe to PATH". Then double-click this file again.
  echo.
  pause
  exit /b 1
)

rem --- first-time setup: private environment inside this folder ---
if not exist ".venv\Scripts\python.exe" (
  echo First-time setup. This needs internet and takes 1-3 minutes...
  %PY% -m venv .venv
  if errorlevel 1 goto :fail
  ".venv\Scripts\python.exe" -m pip install --upgrade pip
  ".venv\Scripts\python.exe" -m pip install -r requirements.txt
  if errorlevel 1 goto :fail
)

rem --- skip Streamlit's one-time email prompt ---
if not exist "%USERPROFILE%\.streamlit\credentials.toml" (
  if not exist "%USERPROFILE%\.streamlit" mkdir "%USERPROFILE%\.streamlit"
  > "%USERPROFILE%\.streamlit\credentials.toml" echo [general]
  >> "%USERPROFILE%\.streamlit\credentials.toml" echo email = ""
)

echo.
echo Starting PlasmidGuard. A browser tab will open.
echo Keep this window open while you use it. Close it to stop PlasmidGuard.
echo.
".venv\Scripts\python.exe" -m streamlit run app.py --server.address 127.0.0.1
pause
exit /b 0

:fail
echo.
echo Setup failed. Check your internet connection and try again.
echo If it keeps failing, delete the ".venv" folder here and run this file again.
echo.
pause
exit /b 1
