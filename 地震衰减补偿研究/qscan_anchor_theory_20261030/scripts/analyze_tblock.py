"""Analysis of PLAN_TBLOCK_ZH.md: two-anchor kappa per contiguous block of live traces (3 blocks) on U1320A.

Selection per block: grid minimum of the summed J, vertex of the parabola through it and its neighbours; minimum at a
grid end -> unreliable. Descriptive: beta at v_a = 1 (v_sel - 1), median sigma, RMS ratio 3.2-3.7 s / 2.04-2.8 s,
15-trace sliding windows (step 5). -> analysis/tblock.json
"""
import json
from pathlib import Path
import numpy as np
import u1320a_kappa as uk

ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results20"


def select(v, J):
    jm = int(np.argmin(J))
    if jm in (0, len(v) - 1):
        return None
    a = np.polyfit(v[jm - 1:jm + 2], J[jm - 1:jm + 2], 2)
    return float(np.clip(-a[1] / (2 * a[0]), v[jm - 1], v[jm + 1])) if a[0] > 0 else float(v[jm])


def kap(rec, idx):
    s = {va: select(np.asarray(rec[va]["v"]), np.asarray(rec[va]["J_trace"])[:, idx].sum(axis=1)) for va in rec}
    k = None if s[2.0] is None or s[0.5] is None else (s[2.0] - s[0.5]) / 1.5
    return k, s


out = {}; u18 = json.loads((ROOT / "analysis" / "u1320a.json").read_text())
for line in uk.LINES:
    D, cm = uk.load(line); t = np.asarray(cm["t"])
    for w in uk.WINS:
        rec = {va: json.loads((R / f"{line}__{w}_va{va}.json").read_text()) for va in uk.VAS}
        live = np.flatnonzero(rec[1.0]["live"]); sig = np.asarray(rec[1.0]["sigma"])
        kall, sall = kap(rec, live)
        old = {va: json.loads((ROOT / "results18" / f"{line}__{w}_va{va}.json").read_text())["v_sel"] for va in uk.VAS}
        blocks = []
        for idx in np.array_split(live, 3):
            k, s = kap(rec, idx)
            deep = np.sqrt(np.mean(D[(t >= 3.2), :][:, idx] ** 2)); shal = np.sqrt(np.mean(D[(t <= 2.8), :][:, idx] ** 2))
            blocks.append(dict(traces=[int(idx[0]), int(idx[-1])], n=len(idx), kappa=k, v_sel=s,
                               beta_va1=(None if s[1.0] is None else s[1.0] - 1), sigma_median=float(np.median(sig[idx])),
                               rms_deep_over_shallow=float(deep / shal)))
        slide = []
        for i0 in range(0, len(live) - 14, 5):
            idx = live[i0:i0 + 15]; k, s = kap(rec, idx); slide.append(dict(first=int(idx[0]), kappa=k, v_va1=s[1.0]))
        out[f"{line}|{w}"] = dict(kappa_line=kall, kappa_results18=u18[f"{line}|{w}"]["kappa"],
                                  repro_max=max(abs(sall[va] - old[va]) for va in uk.VAS), n_live=len(live),
                                  all_converged=all(all(rec[va]["converged"]) for va in uk.VAS), blocks=blocks, sliding=slide)

ctrl = [b["kappa"] for b in out["3677-32|w3.00"]["blocks"]]; tgt = [b["kappa"] for b in out["4508-32|w3.00"]["blocks"]]
C0 = None if None in ctrl else all(abs(k - 0.816) <= 0.15 for k in ctrl)
if C0 is not True or None in tgt:
    B1 = B2 = None
else:
    B1 = (max(tgt) - min(tgt) > 0.20) and any(abs(k - 0.884) <= 0.15 for k in tgt); B2 = all(k <= 0.75 for k in tgt)
out["verdicts"] = dict(C0=C0, B1=B1, B2=B2)
(ROOT / "analysis" / "tblock.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for key, r in out.items():
    if key == "verdicts":
        continue
    print(key, "line %.3f (r18 %.3f) repro %.4f conv %s" % (r["kappa_line"], r["kappa_results18"], r["repro_max"], r["all_converged"]))
    for b in r["blocks"]:
        print("   traces %s kappa %s beta %s sigma %.3g rms %.3f" % (b["traces"], None if b["kappa"] is None else round(b["kappa"], 3),
              None if b["beta_va1"] is None else round(b["beta_va1"], 2), b["sigma_median"], b["rms_deep_over_shallow"]))
    print("   sliding", [None if s["kappa"] is None else round(s["kappa"], 2) for s in r["sliding"]])
print(out["verdicts"])
