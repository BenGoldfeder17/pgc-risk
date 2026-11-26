@echo off
rem Drag a genome file onto this file, or run:  pgc-risk.bat <genome file>
setlocal
if "%~1"=="" (
    echo Drag a genome file ^(23andMe / AncestryDNA / MyHeritage / FTDNA export or a VCF^) onto this file,
    echo or run:  pgc-risk.bat ^<genome file^>
    echo.
    "%~dp0.venv\Scripts\python.exe" -m pgc_risk --help
) else (
    "%~dp0.venv\Scripts\python.exe" -m pgc_risk %*
)
rem keep the window open when started by double-click or drag-and-drop
echo %cmdcmdline% | findstr /i /c:"%~nx0" >nul && pause
endlocal
