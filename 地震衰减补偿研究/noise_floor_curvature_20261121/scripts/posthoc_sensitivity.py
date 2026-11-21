"""Post-hoc, descriptive (not in PLAN): on the line-mean spectra, multiply the M1 noise estimate of the deep window only
(or of both windows) by k and report the c2 drop of the line-mean ratio. Which k would reach R = 0.5?"""
import json, numpy as np
from pathlib import Path
R = json.loads((Path(__file__).resolve().parents[1] / "results" / "run.json").read_text(encoding="utf-8"))
F = np.array(R["F"]); FIT = (F >= 10) & (F <= 55); M = R["mean_spectra"]
P = {w: np.array(M[w]["P"]) for w in ("sh", "dp")}; N = {w: np.array(M[w]["N1"]) for w in ("sh", "dp")}
c2 = lambda v: float(np.polyfit(F[FIT], v[FIT], 2)[0])
c0 = c2(0.5 * np.log(P["dp"] / P["sh"])); print("line-mean ratio c2 %+.3e" % c0)
for mode in ("deep only", "both"):
    for k in (1, 2, 3, 5, 10, 20, 50):
        Ps = {w: np.maximum(P[w] - k * N[w] * (1 if (mode == "both" or w == "dp") else 1 / k), 0.1 * k * N[w]) for w in P}
        c = c2(0.5 * np.log(Ps["dp"] / Ps["sh"])); print("%-9s k=%3d  c2 %+.3e  drop %+.2f" % (mode, k, c, (c0 - c) / c0))
