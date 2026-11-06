"""Verdicts of STEP2B_ZH.md -> results/step2b_verdict.json"""
import gzip, json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"; R = []
for m in ("M12", "M13", "M14", "M15"):
    R += json.loads(gzip.open(OUT / f"step2b_{m}.json.gz", "rt", encoding="utf-8").read())
K = ("gac", "ga", "eg", "egl", "mg", "mglc", "cg", "wave", "lam", "ba", "bg"); M = ("sat", "split", "method", "direct"); idx = {}
for r in R:
    idx.setdefault((r["M"], r["fp"], r["scene"], r["snr"], r["qf"], r["m"]), []).append(r)
cells = sorted({k[:5] for k in idx}); T = {c: {m: {k: float(np.mean([x[k] for x in idx[c + (m,)]])) for k in K} for m in M} for c in cells}
for c in cells:
    for m in M:
        T[c][m]["fb"] = idx[c + (m,)][0]["fallback"]
ali = [c for c in cells if c[2].endswith("aligned")]; nmo = [c for c in cells if c[2].endswith("nmo")]; ev = lambda c: c[2].split("_")[0] in ("events", "thin")
AQ = [c for c in ali if c[4] == 1.0]; AQE = [c for c in AQ if ev(c)]; AW = [c for c in ali if c[4] != 1.0]; AWE = [c for c in AW if ev(c)]
med = lambda sel, m, k: float(np.median([T[c][m][k] for c in sel]))
V = {}
V["R1"] = med(AQE, "sat", "gac") <= 0.5 * med(AQE, "split", "gac"); fr2 = float(np.mean([T[c]["sat"]["gac"] < T[c]["split"]["gac"] for c in AQE])); V["R2"] = fr2 >= 0.8; V["R3"] = med(AQE, "sat", "ga") <= med(AQE, "split", "ga"); V["R4"] = med(AQ, "sat", "wave") < med(AQ, "split", "wave"); V["R5"] = med(AQ, "sat", "eg") <= 1.1 * med(AQ, "split", "eg")
V["R6a"] = med(AW, "sat", "wave") <= med(AW, "split", "wave"); V["R6b"] = med(AWE, "sat", "ga") <= 1.1 * med(AWE, "split", "ga"); V["R7a"] = med(nmo, "sat", "eg") <= 1.1 * med(nmo, "split", "eg"); V["R7b"] = med(nmo, "sat", "wave") <= med(nmo, "split", "wave") + 0.01
rat = [T[c]["sat"]["wave"] / T[c]["split"]["wave"] for c in cells]; V["R8"] = max(rat) <= 1.5; V["R9"] = med(AQE, "sat", "gac") <= med(AQE, "method", "gac") and med(AQE, "sat", "gac") <= med(AQE, "direct", "gac")
print(f"cells: aligned Q correct events+thin {len(AQE)}, aligned Q correct {len(AQ)}, aligned Q wrong {len(AW)}, events+thin {len(AWE)}, NMO {len(nmo)}, all {len(cells)}")
print(f"R1 noise-free per-event G/A median: sat {med(AQE, 'sat', 'gac'):.4f}, split {med(AQE, 'split', 'gac'):.4f} (<= 0.5 x split) -> {'pass' if V['R1'] else 'FAIL'}")
print(f"R2 sat better than split in {fr2:.1%} of cells (>= 80%) -> {'pass' if V['R2'] else 'FAIL'}")
print(f"R3 noisy per-event G/A: sat {med(AQE, 'sat', 'ga'):.4f}, split {med(AQE, 'split', 'ga'):.4f} (<=) -> {'pass' if V['R3'] else 'FAIL'}")
print(f"R4 wave aligned Q correct: sat {med(AQ, 'sat', 'wave'):.3f}, split {med(AQ, 'split', 'wave'):.3f} (<) -> {'pass' if V['R4'] else 'FAIL'}")
print(f"R5 E_G aligned Q correct: sat {med(AQ, 'sat', 'eg'):.3f}, split {med(AQ, 'split', 'eg'):.3f} (<= 1.1 x) -> {'pass' if V['R5'] else 'FAIL'}")
print(f"R6a wave aligned Q wrong: sat {med(AW, 'sat', 'wave'):.3f}, split {med(AW, 'split', 'wave'):.3f} (<=) -> {'pass' if V['R6a'] else 'FAIL'}")
print(f"R6b noisy G/A Q wrong: sat {med(AWE, 'sat', 'ga'):.3f}, split {med(AWE, 'split', 'ga'):.3f} (<= 1.1 x) -> {'pass' if V['R6b'] else 'FAIL'}")
print(f"R7a E_G NMO: sat {med(nmo, 'sat', 'eg'):.3f}, split {med(nmo, 'split', 'eg'):.3f} (<= 1.1 x) -> {'pass' if V['R7a'] else 'FAIL'}")
print(f"R7b wave NMO: sat {med(nmo, 'sat', 'wave'):.3f}, split {med(nmo, 'split', 'wave'):.3f} (<= split + 0.01) -> {'pass' if V['R7b'] else 'FAIL'}")
print(f"R8 worst wave ratio sat/split over all cells {max(rat):.3f} (<= 1.5) -> {'pass' if V['R8'] else 'FAIL'}")
print(f"R9 noise-free per-event G/A: method {med(AQE, 'method', 'gac'):.4f}, direct {med(AQE, 'direct', 'gac'):.4f}, sat not above either -> {'pass' if V['R9'] else 'FAIL'}")
print("OVERALL:", "PASS" if all(V.values()) else "FAIL")


def show(sel, title):
    print(f"\n{title} n={len(sel)}\n  {'':7s}" + "".join(f"{k:>8s}" for k in K))
    for m in M:
        print(f"  {m:7s}" + "".join((f"{med(sel, m, k):8.3f}" if k != "lam" else f"{med(sel, m, k):8.1e}") for k in K))


show(AQE, "aligned Q correct events+thin"); show(AQ, "aligned Q correct all"); show(AW, "aligned Q wrong"); show([c for c in ali if c[4] == 0.8], "aligned Q x0.8"); show([c for c in ali if c[4] == 1.2], "aligned Q x1.2"); show(nmo, "NMO all"); show([c for c in nmo if c[4] == 1.0], "NMO Q correct")
for s in ("events", "thin", "dense"):
    show([c for c in AQ if c[2].startswith(s)], f"aligned Q correct {s}")
for fp in (60.0, 40.0):
    show([c for c in AQ if c[1] == fp], f"aligned Q correct {fp} Hz"); show([c for c in AQE if c[1] == fp], f"aligned Q correct events+thin {fp} Hz")
for mo in ("M12", "M13", "M14", "M15"):
    show([c for c in AQ if c[0] == mo], f"aligned Q correct {mo}"); show([c for c in AQE if c[0] == mo], f"aligned Q correct events+thin {mo}")
for snr in (18, 12):
    show([c for c in AQE if c[3] == snr], f"aligned Q correct events+thin {snr} dB")
print("\nFallback seeds:", {m: int(sum(T[c][m]["fb"] for c in cells)) for m in M}); w = np.argsort(rat)[-3:]; print("worst wave ratio cells:", [(cells[i], round(rat[i], 3)) for i in w])
print("signed noise-free per-event G/A, aligned Q correct events+thin, mean over rows:")
for m in M:
    x = [r for r in R if r["m"] == m and r["qf"] == 1.0 and r["scene"] in ("events_aligned", "thin_aligned")]
    print(f"  {m:7s}" + " ".join(f"{np.mean([r[f'gac{j}'] for r in x]):+.4f}" for j in range(3)))
(OUT / "step2b_verdict.json").write_text(json.dumps(dict(verdict={k: bool(v) for k, v in V.items()}, fr2=fr2, worst_wave_ratio=max(rat), cells={"|".join(map(str, c)): T[c] for c in cells}), indent=1), encoding="utf-8")
