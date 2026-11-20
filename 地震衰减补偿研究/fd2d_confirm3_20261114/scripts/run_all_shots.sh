#!/bin/bash
# Worker K of 4: shot chunks (model m = 1..6, 25 shots each, shots 0-99) taken round-robin; make_shots skips existing shots.
K=$1; cd "$(dirname "$0")"; i=0
for m in 1 2 3 4 5 6; do for a in 0 25 50 75; do
  if [ $((i % 4)) -eq $K ]; then python3 run_octave.py make_shots $m $a $((a+24)) || echo "FAILED make_shots $m $a"; fi; i=$((i+1))
done; done
echo "WORKER $K DONE"
