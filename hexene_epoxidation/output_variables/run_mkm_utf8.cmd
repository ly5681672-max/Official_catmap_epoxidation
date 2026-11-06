@echo off
setlocal

set "SOURCE_DIR=%~dp0"
set "RUNS_ROOT=%SOURCE_DIR%runs"

if not exist "%RUNS_ROOT%" mkdir "%RUNS_ROOT%"
if errorlevel 1 (
    echo ERROR - cannot create runs folder: %RUNS_ROOT%
    pause
    exit /b 1
)

for /f %%I in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss_fff"') do set "RUN_ID=%%I"
if not defined RUN_ID (
    echo ERROR - cannot generate run ID.
    pause
    exit /b 1
)

set "RUN_DIR=%RUNS_ROOT%\run_%RUN_ID%"
if exist "%RUN_DIR%" (
    echo ERROR - run folder already exists: %RUN_DIR%
    pause
    exit /b 1
)
mkdir "%RUN_DIR%"
if errorlevel 1 (
    echo ERROR - cannot create run folder: %RUN_DIR%
    pause
    exit /b 1
)

for %%F in (alkene_epoxidation.mkm energies_dft_neb.txt alkene_epoxidation_dft.pkl mkm_job.py descriptor_output.py tof_export.py) do (
    if not exist "%SOURCE_DIR%%%F" (
        echo ERROR - required input is missing: %SOURCE_DIR%%%F
        pause
        exit /b 1
    )
    copy /Y "%SOURCE_DIR%%%F" "%RUN_DIR%\%%F" >nul
    if errorlevel 1 (
        echo ERROR - failed to snapshot input: %%F
        pause
        exit /b 1
    )
)

set "FULL_LOG=%RUN_DIR%\alkene_epoxidation.log"
echo RUNNING - CatMAP calculation is in progress.
echo RUN_ID=%RUN_ID%
echo RUN_DIR=%RUN_DIR%
echo FULL_LOG=%FULL_LOG%

pushd "%RUN_DIR%"
if errorlevel 1 (
    echo ERROR - cannot enter run folder: %RUN_DIR%
    pause
    exit /b 1
)

set "PYTHONPATH=H:\catmap\Official_catmap_epoxidation\catmap_official_source"
"G:\anaconda\python.exe" -X utf8 -u mkm_job.py > "%FULL_LOG%" 2>&1
set "RUN_EXIT_CODE=%ERRORLEVEL%"
popd

>> "%FULL_LOG%" echo.
>> "%FULL_LOG%" echo ====================================
>> "%FULL_LOG%" echo RUN_ID: %RUN_ID%
>> "%FULL_LOG%" echo CatMAP execution exit code: %RUN_EXIT_CODE%
>> "%FULL_LOG%" echo ====================================

echo FINISHED - CatMAP exit code: %RUN_EXIT_CODE%
echo Results and logs: %RUN_DIR%
pause
exit /b %RUN_EXIT_CODE%