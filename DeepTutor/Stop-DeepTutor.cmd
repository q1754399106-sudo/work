@echo off
echo For console launches, press Ctrl+C in the DeepTutor console to stop.
echo This command checks instances started with the official --detach option.
cd /d "%~dp0"
set PYTHONUTF8=1
set PYTHONIOENCODING=utf-8
".venv\Scripts\python.exe" -m deeptutor stop
pause
