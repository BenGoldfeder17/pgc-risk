@echo off
rem Double-click to open pgc-risk for this computer AND other devices on the same network.
rem The window shows the address to open on the other device and the access code to enter.
rem Other devices connect over HTTPS and must enter the code. Close this window to stop it.
setlocal
title pgc-risk (network)
"%~dp0.venv\Scripts\python.exe" -m pgc_risk gui --lan %*
if errorlevel 1 pause
endlocal
