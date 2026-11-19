#!/bin/bash
# Stage 4 field run (PLAN_STAGE4_ZH.md): fld_run in one MATLAB process at high priority.
cd "$(dirname "$0")"
( while true; do powershell -NoProfile -Command "Get-Process MATLAB -ErrorAction SilentlyContinue | ForEach-Object { \$_.PriorityClass = 'High' }" >/dev/null 2>&1; sleep 300; done ) &
PRIO=$!
echo "=== stage 4 start $(date '+%F %T')"
"/d/Program Files (x86)/MATLAB2024B/bin/matlab" -batch "cd scripts; fld_run"
echo "=== stage 4 done $(date '+%F %T') exit $?"
kill $PRIO
