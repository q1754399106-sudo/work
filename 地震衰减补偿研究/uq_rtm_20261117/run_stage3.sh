#!/bin/bash
# Stage 3 (PLAN_STAGE3_ZH.md): C1-C4 shots in two processes, then noise estimate + migration + scores in one, then verdict.
cd "$(dirname "$0")"
M="/d/Program Files (x86)/MATLAB2024B/bin/matlab"
( while true; do powershell -NoProfile -Command "Get-Process MATLAB -ErrorAction SilentlyContinue | ForEach-Object { \$_.PriorityClass = 'High' }" >/dev/null 2>&1; sleep 300; done ) &
PRIO=$!
echo "=== stage 3 start $(date '+%F %T')"
"$M" -batch "cd scripts; fd_shots('C1',1:51,8); fd_shots('C3',1:51,8)" > logs/fd_C1_C3.log 2>&1 &
FD1=$!
"$M" -batch "cd scripts; fd_shots('C2',1:51,8); fd_shots('C4',1:51,8)" > logs/fd_C2_C4.log 2>&1 &
FD2=$!
wait $FD1 $FD2
echo "=== shots done $(date '+%F %T')"
"$M" -batch "cd scripts; run_stage2({'C1','C2','C3','C4'})"
echo "=== migration done $(date '+%F %T') exit $?"
"/c/Users/qq175/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe" scripts/judge_stage3.py
echo "=== stage 3 done $(date '+%F %T')"
kill $PRIO
