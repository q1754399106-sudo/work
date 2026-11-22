"""Verdict of PLAN_BLOCKFLIP_ZH.md section 4 (as fixed in PLAN_BLOCKFLIP_STEP2_ZH.md) from
analysis/stage2_<D>_bf1.json and _bf300.json (deep window, E1). d = UQW1 - B+ with the same 1/snr;
M1 = mean over m < 1, M2 = m = 1, M3 = worst over m > 1 (fields M1_u1, M2_u1, M3_u1 of score_stage2).
-> analysis/blockflip_B.json"""
import json
from pathlib import Path
import numpy as np
import sys
ROOT = Path(__file__).resolve().parents[1]; MS = sys.argv[1:] or ["D1", "D2", "D3", "D4"]; BS = 300; SETS = ("m30", "m40", "l30")
OUT = "blockflip_B.json" if MS[0].startswith("D") else "blockflip_C.json"   # confirmation on C5-C8: same rules
METH = ("R0", "PH", "IQ30", "IQ70", "LP", "W0", "W1", "UQW1")


def load(m, b):
    return json.loads((ROOT / "analysis" / ("stage2_%s_bf%d.json" % (m, b))).read_text(encoding="utf-8"))


J = {(m, b): load(m, b) for m in MS for b in (1, BS)}
T = {(m, b, s): J[(m, b)][s]["deep"]["E1"] for m in MS for b in (1, BS) for s in SETS}
val = lambda m, b, s, q: float(np.atleast_1d(T[(m, b, s)][q])[0])
out = {"cells": {}}
for m in MS:
    for s in SETS:
        c = {}
        for b in (1, BS):
            t = T[(m, b, s)]
            c["B%d" % b] = dict(M1=val(m, b, s, "M1_u1"), M2=val(m, b, s, "M2_u1"), M3=val(m, b, s, "M3_u1"),
                                fc=J[(m, b)]["fc"][s], **{q: list(np.atleast_1d(t[q]).astype(float)) for q in METH + ("Bplus",)})
        out["cells"]["%s_%s" % (m, s)] = c
dl = lambda s, q: [val(m, BS, s, q) - val(m, 1, s, q) for m in MS]
v = {}
v["B1_l30_dM2_mean"] = float(np.mean(dl("l30", "M2_u1"))); v["B1_l30_dM1_mean"] = float(np.mean(dl("l30", "M1_u1")))
v["B1_pass"] = bool(v["B1_l30_dM2_mean"] >= 0.03 and v["B1_l30_dM1_mean"] >= 0.03)
ok = True
for s in ("m30", "m40"):
    v["B2_%s_dM1_mean" % s] = float(np.mean(dl(s, "M1_u1"))); ok &= v["B2_%s_dM1_mean" % s] >= -0.02
    v["B2_%s_M2_min" % s] = float(min(val(m, BS, s, "M2_u1") for m in MS)); ok &= v["B2_%s_M2_min" % s] >= -0.04
    v["B2_%s_M3_min" % s] = float(min(val(m, BS, s, "M3_u1") for m in MS)); ok &= v["B2_%s_M3_min" % s] >= -0.06
v["B2_pass"] = bool(ok); v["pass"] = bool(v["B1_pass"] and v["B2_pass"]); out["verdict"] = v
(ROOT / "analysis" / OUT).write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps(v, indent=1))
for k, c in out["cells"].items():
    print(k, "  ".join("B%d M1 %+.3f M2 %+.3f M3 %+.3f fc %.1f | UQW1@m1 %.3f B+@m1 %.3f" %
                       (b, c["B%d" % b]["M1"], c["B%d" % b]["M2"], c["B%d" % b]["M3"], c["B%d" % b]["fc"],
                        c["B%d" % b]["UQW1"][2], c["B%d" % b]["Bplus"][2]) for b in (1, BS)))
