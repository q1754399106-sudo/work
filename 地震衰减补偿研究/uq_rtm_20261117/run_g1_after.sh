#!/bin/bash
# Addendum 1: after the main stage-2 chain, add W1/UQW1 for D1 and D2 and re-score.
cd "$(dirname "$0")"
M="/d/Program Files (x86)/MATLAB2024B/bin/matlab"
until grep -q "=== stage 2 done" logs/stage2_all.log; do sleep 120; done
( while true; do powershell -NoProfile -Command "Get-Process MATLAB -ErrorAction SilentlyContinue | ForEach-Object { \$_.PriorityClass = 'High' }" >/dev/null 2>&1; sleep 300; done ) &
PRIO=$!
"$M" -batch "cd scripts; rerun_g1({'D1','D2'})"
echo "=== g1 done $(date '+%F %T') exit $?"
kill $PRIO
