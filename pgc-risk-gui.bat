@echo off
rem Double-click to open pgc-risk in your browser. It runs on this computer only.
rem Close this window to stop it.
setlocal
title pgc-risk
"%~dp0.venv\Scripts\python.exe" -m pgc_risk gui %*
if errorlevel 1 pause
endlocal
