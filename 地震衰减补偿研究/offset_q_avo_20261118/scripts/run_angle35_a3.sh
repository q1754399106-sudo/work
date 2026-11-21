#!/bin/bash
# Second runner (started 2026-11-21 ~10:00 for throughput): gate A3 jobs on 2 workers, logs angle35b_*; each job drops a .lock
# so that run_angle35.sh skips it when it reaches the a3 part of its queue.
cd "$(dirname "$0")/.."
PY='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
export OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONIOENCODING=utf-8
jobs=""; for fp in 60 40; do for m in M28 M29 M30 M31; do jobs="$jobs a3:$m:$fp"; done; done
echo $jobs | tr ' ' '\n' | xargs -P 2 -I{} bash -c 'IFS=: read g m fp <<< "{}"; L=logs/angle35_${g}_${m}_f${fp}.lock; [ -e $L ] && exit 0; touch $L; ANGLE35_OWNER=1 "$0" scripts/angle35.py $g $m $fp > logs/angle35b_${g}_${m}_f${fp}.log 2> logs/angle35b_${g}_${m}_f${fp}.err' "$PY"
echo ALL DONE
