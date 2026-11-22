#!/bin/bash
# PLAN_BLOCKFLIP_CONFIRM_ZH.md: serial chain on C5-C8, one MATLAB process; each step is resumable and retried
# (10 min wait, up to 12 tries) after a failure such as out-of-memory. Log: logs/blockflip_confirm.log
cd "$(dirname "$0")"; L=../logs/blockflip_confirm.log
step() {  # step <done-file> <matlab command>
  for try in $(seq 1 12); do
    [ -n "$1" ] && [ -f "$1" ] && return 0
    echo "=== $2 try $try $(date '+%F %T') ===" >> $L
    matlab -batch "$2" >> $L 2>&1 && [ -z "$1" ] && return 0
    [ -n "$1" ] && [ -f "$1" ] && return 0
    sleep 600
  done
  echo "GAVE UP: $2" >> $L; return 1
}
C=../cache; A=../analysis
for m in C5 C6 C7 C8; do
  step "$C/shots_$m/s051.mat" "P=rtm_setup();fd_shots('$m',1:numel(P.ixs))" || exit 1
  step "$C/data_$m.mat" "prep_data('$m')" || exit 1
  step "$C/invblk_$m.mat" "est_block('$m',[],'',[1 300])" || exit 1
  step "$A/stage2_${m}_bf1.json" "run_block2('$m',1)" || exit 1
  step "$A/stage2_${m}_bf300.json" "run_block2('$m',300)" || exit 1
done
step "$A/blockflip_C_accuracy.json" "score_block_c" || exit 1
"/c/Users/qq175/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe" judge_block2.py C5 C6 C7 C8 > ../logs/blockflip_C_verdict.log 2>&1
echo "=== chain done $(date '+%F %T') ===" >> $L
