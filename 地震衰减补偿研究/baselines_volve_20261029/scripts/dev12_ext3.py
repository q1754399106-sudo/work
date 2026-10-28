"""Third grid extension of the experiment-12 development (PLAN_ZH.md section 5): after the second extension SC-a was
still best at the smallest weight (a = 1/256), so a goes down to 1/8192 (practically no sparsity term), and the
decoupled variant gets k = 1.15 next to k = 1.3 (k = 1.3 beat 1.6 at every a of the second extension).
Same job and file format as dev12_ext2.py; writes analysis/dev_ext3/<id>__<name>.json.

usage: dev12_ext3.py K N    worker K of N
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dev12_ext2 as E  # noqa: E402

E.A_ONE = [2.0 ** -9, 2.0 ** -11, 2.0 ** -13]
E.A_K, E.KS = [2.0 ** -9, 2.0 ** -11, 2.0 ** -13], [1.15, 1.3]
E.OUT = E.L.ROOT / "analysis" / "dev_ext3"
_jobs = E.jobs
E.jobs = lambda: _jobs() + [(d, 2.0 ** -8, 1.15) for d in sorted({j[0] for j in _jobs()})]  # a = 1/256 had k = 1.3 only

if __name__ == "__main__":
    E.main(int(sys.argv[1]), int(sys.argv[2]))
