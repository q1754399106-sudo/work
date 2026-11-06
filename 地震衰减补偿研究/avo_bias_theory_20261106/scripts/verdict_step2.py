"""Verdicts of STEP2_ZH.md -> results/step2_verdict.json"""
import gzip, json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"; R = []
for m in ("M1", "M10", "M11"):
    R += json.loads(gzip.open(OUT / f"step2_{m}.json.gz", "rt", encoding="utf-8").read())
K = ("eg", "cg", "wave", "mg", "ma", "mb"); M = ("method", "split", "direct"); idx = {}
for r in R:
    idx.setdefault((r["scene"], r["M"], r["snr"], r["qf"], r["m"]), []).append(r)
cells = sorted({k[:4] for k in idx}); T = {c: {m: {k: float(np.mean([x[k] for x in idx[c + (m,)]])) for k in K} for m in M} for c in cells}
for c in cells:
    for m in M:
        T[c][m]["lam"] = float(np.median([x["lam"] for x in idx[c + (m,)]])); T[c][m]["fb"] = idx[c + (m,)][0]["fallback"]
nmo = [c for c in cells if c[0].endswith("nmo")]; ali = [c for c in cells if c[0].endswith("aligned")]
med = lambda sel, m, k: float(np.median([T[c][m][k] for c in sel]))
U1 = med(nmo, "method", "eg") <= 1.1 * med(nmo, "split", "eg"); U2 = med(nmo, "method", "wave") <= med(nmo, "split", "wave"); U3 = med(nmo, "method", "cg") >= med(nmo, "split", "cg") - 0.03
fr = float(np.mean([T[c]["method"]["eg"] <= 1.2 * T[c]["split"]["eg"] for c in nmo])); U4 = fr >= 0.9; U5 = med(nmo, "method", "eg") <= med(nmo, "direct", "eg")
print(f"NMO cells {len(nmo)}")
print(f"U1 E_G median method {med(nmo, 'method', 'eg'):.3f}, split {med(nmo, 'split', 'eg'):.3f}, ratio {med(nmo, 'method', 'eg') / med(nmo, 'split', 'eg'):.3f} (<= 1.1) -> {'pass' if U1 else 'FAIL'}")
print(f"U2 wave median method {med(nmo, 'method', 'wave'):.3f}, split {med(nmo, 'split', 'wave'):.3f} -> {'pass' if U2 else 'FAIL'}")
print(f"U3 corr G method {med(nmo, 'method', 'cg'):.3f}, split {med(nmo, 'split', 'cg'):.3f} (>= split - 0.03) -> {'pass' if U3 else 'FAIL'}")
print(f"U4 E_G within 20% of split in {fr:.1%} (>= 90%) -> {'pass' if U4 else 'FAIL'}")
print(f"U5 E_G direct {med(nmo, 'direct', 'eg'):.3f} (method <= direct) -> {'pass' if U5 else 'FAIL'}")
print("OVERALL:", "PASS" if U1 and U2 and U3 and U4 and U5 else "FAIL")


def show(sel, title):
    print(f"\n{title} (median): method / split / direct")
    for k in K:
        print(f"  {k:5s} {med(sel, 'method', k):7.3f} {med(sel, 'split', k):7.3f} {med(sel, 'direct', k):7.3f}")


show(nmo, "NMO all"); show(ali, "aligned all")
for s in ("events", "thin", "dense"):
    show([c for c in nmo if c[0].startswith(s)], f"NMO {s}"); show([c for c in ali if c[0].startswith(s)], f"aligned {s}")
show([c for c in nmo if c[3] == 1.0], "NMO Q correct"); show([c for c in ali if c[3] == 1.0], "aligned Q correct"); show([c for c in nmo if c[3] != 1.0], "NMO Q wrong"); show([c for c in ali if c[3] != 1.0], "aligned Q wrong")
for mo in ("M1", "M10", "M11"):
    show([c for c in nmo if c[1] == mo], f"NMO {mo}")
print("\nFallback to smallest lambda (seed count over all cells):", {m: int(sum(T[c][m]["fb"] for c in cells)) for m in M}, "; median chosen lambda NMO:", {m: med(nmo, m, "lam") for m in M}, "aligned:", {m: med(ali, m, "lam") for m in M})
p = OUT / "step2_stress.json"
if p.exists():
    S = json.loads(p.read_text(encoding="utf-8")); print("Velocity stress, M1, Q correct, 18 dB: E_G / corr G / wave")
    for sc in ("events_nmo", "dense_nmo"):
        for vs in (0.98, 1.02):
            print(f"  {sc} v x{vs}: " + " | ".join(f"{m} {np.mean([r['eg'] for r in S if r['scene'] == sc and r['vs'] == vs and r['m'] == m]):.3f} {np.mean([r['cg'] for r in S if r['scene'] == sc and r['vs'] == vs and r['m'] == m]):.3f} {np.mean([r['wave'] for r in S if r['scene'] == sc and r['vs'] == vs and r['m'] == m]):.3f}" for m in M))
(OUT / "step2_verdict.json").write_text(json.dumps(dict(U1=bool(U1), U2=bool(U2), U3=bool(U3), U4=bool(U4), U5=bool(U5), within20=fr, cells={"|".join(map(str, c)): T[c] for c in cells}), indent=1), encoding="utf-8")
