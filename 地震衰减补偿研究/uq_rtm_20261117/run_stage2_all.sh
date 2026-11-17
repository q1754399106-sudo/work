#!/bin/bash
# Wait for the D1-D4 shots, smoke-test stage 2 on D1 (3 frequencies), then run stage 2 in full.
cd "$(dirname "$0")"
M="/d/Program Files (x86)/MATLAB2024B/bin/matlab"
until [ "$(ls cache/shots_D1 cache/shots_D2 cache/shots_D3 cache/shots_D4 2>/dev/null | grep -c '^s0')" -ge 204 ] && ! powershell -NoProfile -Command "if (Get-Process MATLAB -ErrorAction SilentlyContinue) {exit 0} else {exit 1}"; do sleep 60; done
( while true; do powershell -NoProfile -Command "Get-Process MATLAB -ErrorAction SilentlyContinue | ForEach-Object { \$_.PriorityClass = 'High' }" >/dev/null 2>&1; sleep 300; done ) &
PRIO=$!
echo "=== smoke $(date '+%F %T')"
UQRTM_SMOKE=1 "$M" -batch "cd scripts; run_stage2({'D1'})" || { echo "SMOKE FAILED"; kill $PRIO; exit 1; }
rm -f cache/data_D1.mat cache/tdata_D1.mat      # smoke datasets hold only 3 frequencies
echo "=== smoke ok $(date '+%F %T')"
"$M" -batch "cd scripts; run_stage2({'D1','D2','D3','D4'}); run_stage2({'D1','D3'},[0.96 1.04])"
echo "=== stage 2 done $(date '+%F %T') exit $?"
kill $PRIO
