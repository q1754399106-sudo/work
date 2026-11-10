#!/usr/bin/env bash
# Select (6 workers), pick the best (a, mu), evaluate it on the other 16 datasets (6 workers). Resumable.
cd "$(dirname "$0")/.." || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="C:/Users/qq175/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"; LP=../qunknown_comp_20261109/scripts/lowprio.py
for k in 0 1 2 3 4 5; do "$PY" $LP scripts/pilot_sc_tune.py select $k 6 > logs/pilot_sc_tune_sel6_$k.log 2> logs/pilot_sc_tune_sel6_$k.err & done; wait
read A MU < <("$PY" scripts/pilot_sc_tune.py best 2> logs/pilot_sc_tune_best.txt); echo "best a=$A mu=$MU" | tee -a logs/pilot_sc_tune_best.txt
for k in 0 1 2 3 4 5; do "$PY" $LP scripts/pilot_sc_tune.py eval $k 6 $A $MU > logs/pilot_sc_tune_eval_$k.log 2> logs/pilot_sc_tune_eval_$k.err & done; wait
echo "eval files: $(ls results/pilot_sc_tune/eval_* | wc -l)"
