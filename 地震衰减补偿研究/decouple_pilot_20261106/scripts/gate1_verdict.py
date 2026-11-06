"""Applies the rules of GATE1_ZH.md to results/gate1_M*.json and prints the descriptive tables. -> results/gate1_summary.json"""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"
rows = []
for m in ("M1", "M2", "M3"):
    rows += json.loads((OUT / f"gate1_{m}.json").read_text(encoding="utf-8"))
KEYS = ("wave", "mg", "ma", "mb", "cg"); METHODS = ["split", "direct", "consistent40", "consistent60", "consistent80", "direct_disc", "covsplit_disc", "observed"]
idx = {}
for r in rows:
    idx.setdefault((r["scene"], r["M"], r["snr"], r["qf"], r["m"], r["lam"]), []).append(r)
cells = sorted({k[:4] for k in idx})


def best(cell, m):
    """The method's lambda with the smallest mean waveform error in the cell; means of all metrics there."""
    cand = [(k[5], v) for k, v in idx.items() if k[:4] == cell and k[4] == m]
    lam, v = min(cand, key=lambda z: np.mean([x["wave"] for x in z[1]]))
    return dict(lam=lam, **{k: float(np.mean([x[k] for x in v])) for k in KEYS})


T = {c: {m: best(c, m) for m in METHODS} for c in cells}
nmo = [c for c in cells if c[0].endswith("nmo")]; ali = [c for c in cells if c[0].endswith("aligned")]
ms = np.array([T[c]["split"]["mg"] for c in nmo]); md = np.array([T[c]["direct"]["mg"] for c in nmo])
frac = float(np.mean(ms < md)); ratio = float(np.median(ms / md)); diff = float(np.median(md - ms))
P1 = frac >= 0.8 and ratio <= 0.67; P2 = diff >= 0.05
print(f"NMO cells {len(nmo)}: split better in {frac:.1%}; median M_G ratio split/direct {ratio:.3f}; median difference direct - split {diff:+.3f}")
print(f"P1 {'pass' if P1 else 'fail'}; P2 {'pass' if P2 else 'fail'}; GATE 1 {'PASS' if P1 and P2 else 'FAIL'}\n")


def group(sel, label):
    print(label); print(f"{'method':15s} {'M_G':>6} {'|bA-1|':>7} {'|bG-1|':>7} {'corr G':>7} {'wave':>6}")
    for m in METHODS:
        v = {k: float(np.median([T[c][m][k] for c in sel])) for k in KEYS}; print(f"{m:15s} {v['mg']:6.3f} {v['ma']:7.3f} {v['mb']:7.3f} {v['cg']:7.3f} {v['wave']:6.3f}")
    print()


group(nmo, "Median over all NMO cells"); group(ali, "Median over all aligned cells")
for sc in ("events", "thin", "dense"):
    group([c for c in nmo if c[0].startswith(sc)], f"NMO, {sc}"); group([c for c in ali if c[0].startswith(sc)], f"aligned, {sc}")
group([c for c in nmo if c[3] == 1.0], "NMO, Q correct"); group([c for c in nmo if c[3] != 1.0], "NMO, Q wrong (x0.8, x1.2)")
print("Split better than direct by scene (M_G), fraction of cells:")
for sc in sorted({c[0] for c in cells}):
    s = [c for c in cells if c[0] == sc]; print(f"  {sc:16s} {np.mean([T[c]['split']['mg'] < T[c]['direct']['mg'] for c in s]):.2f}; median split {np.median([T[c]['split']['mg'] for c in s]):.3f}, direct {np.median([T[c]['direct']['mg'] for c in s]):.3f}; waveform split {np.median([T[c]['split']['wave'] for c in s]):.3f}, direct {np.median([T[c]['direct']['wave'] for c in s]):.3f}")
print("\nBy model, NMO cells (median M_G split / direct):")
for mo in ("M1", "M2", "M3"):
    s = [c for c in nmo if c[1] == mo]; print(f"  {mo}: {np.median([T[c]['split']['mg'] for c in s]):.3f} / {np.median([T[c]['direct']['mg'] for c in s]):.3f}; split better {np.mean([T[c]['split']['mg'] < T[c]['direct']['mg'] for c in s]):.2f}")
print("\nDense NMO, signed gains, median b_A, b_G per method at best waveform:")
sd = [c for c in nmo if c[0].startswith("dense")]
for m in ("split", "direct", "consistent60", "observed"):
    B = []
    for c in sd:
        lam = T[c][m]["lam"]; v = idx[c + (m, lam)]; B.append((np.mean([x["ba"] for x in v]), np.mean([x["bg"] for x in v])))
    print(f"  {m:13s} b_A {np.median([b[0] for b in B]):.3f}  b_G {np.median([b[1] for b in B]):.3f}")
(OUT / "gate1_summary.json").write_text(json.dumps(dict(frac=frac, ratio=ratio, diff=diff, P1=P1, P2=P2, cells={"|".join(map(str, c)): T[c] for c in cells}), indent=1), encoding="utf-8")
