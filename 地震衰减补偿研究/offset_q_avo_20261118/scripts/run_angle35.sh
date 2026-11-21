#!/bin/bash
# PLAN_ANGLE35_ZH.md: 16 jobs (gate A2 first, longest), 4 workers, 2 BLAS threads each
cd "$(dirname "$0")/.."
PY='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
export OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONIOENCODING=utf-8
jobs=""; for g in a2 a3; do for fp in 60 40; do for m in $( [ $g = a2 ] && echo M24 M25 M26 M27 || echo M28 M29 M30 M31 ); do jobs="$jobs $g:$m:$fp"; done; done; done
echo $jobs | tr ' ' '\n' | xargs -P 4 -I{} bash -c 'IFS=: read g m fp <<< "{}"; "$0" scripts/angle35.py $g $m $fp > logs/angle35_${g}_${m}_f${fp}.log 2> logs/angle35_${g}_${m}_f${fp}.err' "$PY"
echo ALL DONE
