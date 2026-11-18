#!/bin/bash
# Stage 2 remainder after the restart (PLAN_STAGE2_ADD1_ZH.md): D4, velocity error (D1, D3), gamma-1 re-run of D1-D3.
cd "$(dirname "$0")"
M="/d/Program Files (x86)/MATLAB2024B/bin/matlab"
( while true; do powershell -NoProfile -Command "Get-Process MATLAB -ErrorAction SilentlyContinue | ForEach-Object { \$_.PriorityClass = 'High' }" >/dev/null 2>&1; sleep 300; done ) &
PRIO=$!
echo "=== restart $(date '+%F %T')"
"$M" -batch "cd scripts; run_stage2({'D4'}); run_stage2({'D1','D3'},[0.96 1.04]); rerun_g1({'D1','D2','D3'})"
echo "=== stage 2 done $(date '+%F %T') exit $?"
kill $PRIO
