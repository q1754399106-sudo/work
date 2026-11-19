#!/bin/bash
# Stage 4 post hoc check (PLAN_STAGE4_POSTHOC_ZH.md).
cd "$(dirname "$0")"
( while true; do powershell -NoProfile -Command "Get-Process MATLAB -ErrorAction SilentlyContinue | ForEach-Object { \$_.PriorityClass = 'High' }" >/dev/null 2>&1; sleep 300; done ) &
PRIO=$!
echo "=== post start $(date '+%F %T')"
"/d/Program Files (x86)/MATLAB2024B/bin/matlab" -batch "cd scripts; fld_post; fld_eval('images_fref25.mat','_fref25'); fld_eval('images_amp.mat','_amp')"
echo "=== post done $(date '+%F %T') exit $?"
kill $PRIO
