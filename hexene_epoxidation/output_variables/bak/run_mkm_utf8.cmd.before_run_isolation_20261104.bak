@echo off
setlocal

set "FULL_LOG=%~dp0alkene_epoxidation.log"

echo RUNNING - CatMAP calculation is in progress. Details: alkene_epoxidation.log
call :run_catmap > "%FULL_LOG%" 2>&1
set "RUN_EXIT_CODE=%ERRORLEVEL%"

>> "%FULL_LOG%" echo.
>> "%FULL_LOG%" echo ====================================
>> "%FULL_LOG%" echo CatMAP execution exit code: %RUN_EXIT_CODE%
>> "%FULL_LOG%" echo ====================================

echo FINISHED - CatMAP exit code: %RUN_EXIT_CODE%
pause

exit /b %RUN_EXIT_CODE%

:run_catmap
cd /d "%~dp0"

set "PYTHONPATH=H:\catmap\Official_catmap_epoxidation\catmap_official_source"

"G:\anaconda\python.exe" -X utf8 -u mkm_job.py
exit /b %ERRORLEVEL%
