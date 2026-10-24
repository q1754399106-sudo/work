@echo off
cd /d "%~dp0"
set PYTHONUTF8=1
set DEEPTUTOR_BACKEND_READY_TIMEOUT=300
set DEEPTUTOR_FRONTEND_READY_TIMEOUT=180
set PYTHONIOENCODING=utf-8
powershell.exe -NoProfile -Command "try { $r = Invoke-WebRequest -UseBasicParsing -Uri 'http://127.0.0.1:3782/' -TimeoutSec 3; if ($r.StatusCode -eq 200) { exit 0 } } catch { }; exit 1"
if not errorlevel 1 (
    start "" "http://127.0.0.1:3782/"
    exit /b 0
)
".venv\Scripts\python.exe" -m deeptutor start
if errorlevel 1 pause
