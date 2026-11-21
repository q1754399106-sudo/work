#!/bin/bash
# 2026-11-21 14:20: resume the four jobs killed by out-of-memory (other processes filled RAM around 14:11); per-(c, scene) resume.
cd "$(dirname "$0")/.."
PY='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
export OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONIOENCODING=utf-8 ANGLE35_OWNER=1
for j in a2:M25:60 a2:M27:60 a3:M28:60 a3:M29:60; do IFS=: read g m fp <<< "$j"
  "$PY" scripts/angle35.py $g $m $fp > logs/angle35r_${g}_${m}_f${fp}.log 2> logs/angle35r_${g}_${m}_f${fp}.err &
done; wait; echo ALL DONE
