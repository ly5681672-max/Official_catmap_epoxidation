@echo off
setlocal

cd /d "%~dp0"

set "PYTHONPATH=H:\catmap\Official_catmap_epoxidation\catmap_official_source"

"G:\anaconda\python.exe" -X utf8 -u mkm_job.py

set "RUN_EXIT_CODE=%ERRORLEVEL%"

echo.
echo ====================================
echo CatMAP execution exit code: %RUN_EXIT_CODE%
echo ====================================

pause

exit /b %RUN_EXIT_CODE%