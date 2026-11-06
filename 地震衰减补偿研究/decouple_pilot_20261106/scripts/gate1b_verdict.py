"""Applies GATE1B_ZH.md to results/gate1b_*.json.gz and prints the descriptive tables. -> results/gate1b_summary.json"""
import gzip, json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"; rows = []
for m in ("M1", "M2", "M3"):
    rows += json.loads(gzip.open(OUT / f"gate1b_{m}.json.gz", "rt", encoding="utf-8").read())
KEYS = ("wave", "eg", "mg", "ma", "mb", "cg", "tga"); METHODS = ["split_st", "direct_st", "dest_split", "split", "direct"]


def table(R, key_cols):
    idx = {}
    for r in R:
        idx.setdefault(tuple(r[k] for k in key_cols) + (r["m"], r["lam"]), []).append(r)
    cells = sorted({k[:len(key_cols)] for k in idx}); T = {}
    for c in cells:
        T[c] = {}
        for m in METHODS:
            cand = [(k[-1], v) for k, v in idx.items() if k[:len(key_cols)] == c and k[-2] == m]
            if not cand:
                continue
            lam, v = min(cand, key=lambda z: np.mean([x["wave"] for x in z[1]]))
            T[c][m] = dict(lam=lam, **{k: float(np.mean([x[k] for x in v])) for k in KEYS if k in v[0]})
    return T


T = table(rows, ("scene", "M", "snr", "qf")); nmo = [c for c in T if c[0].endswith("nmo")]; ali = [c for c in T if c[0].endswith("aligned")]
ev = [c for c in nmo if c[0] == "events_nmo"]
a, b = np.array([T[c]["split_st"]["tga"] for c in ev]), np.array([T[c]["split"]["tga"] for c in ev])
tr = float(np.median(a) / np.median(b)); dc = float(np.median([T[c]["split_st"]["cg"] for c in nmo]) - np.median([T[c]["split"]["cg"] for c in nmo])); R1 = tr <= 0.5 and dc >= 0.10
e1 = np.array([T[c]["split_st"]["eg"] for c in nmo]); e2 = np.array([T[c]["direct_st"]["eg"] for c in nmo])
fr = float(np.mean(e1 < e2)); ra = float(np.median(e1 / e2)); di = float(np.median(e2 - e1)); R2 = fr >= 0.8 and ra <= 0.8 and di >= 0.05
print(f"R1: template G/A median split_st/split {tr:.3f} (<= 0.5); corr G gain {dc:+.3f} (>= 0.10) -> {'pass' if R1 else 'fail'}")
print(f"R2: E_G split_st < direct_st in {fr:.1%} (>= 80%); median ratio {ra:.3f} (<= 0.8); median difference {di:+.3f} (>= 0.05) -> {'pass' if R2 else 'fail'}")
print("Outcome:", "PASS" if R1 and R2 else ("R1 only: stretch correction helps, no split advantage" if R1 else "R1 fails"), "\n")


def show(sel, title):
    print(title); print(f"{'method':11s} {'E_G':>6} {'tmpl G/A':>9} {'corr G':>7} {'M_G':>6} {'|bA-1|':>7} {'|bG-1|':>7} {'wave':>6}")
    for m in METHODS:
        c = [T[x][m] for x in sel if m in T[x]]
        if c:
            g = lambda k: float(np.median([y[k] for y in c])) if k in c[0] else float("nan")
            print(f"{m:11s} {g('eg'):6.3f} {g('tga'):9.3f} {g('cg'):7.3f} {g('mg'):6.3f} {g('ma'):7.3f} {g('mb'):7.3f} {g('wave'):6.3f}")
    print()


show(nmo, "All NMO cells (median)"); show([c for c in nmo if c[3] == 1.0], "NMO, Q correct"); show([c for c in nmo if c[3] != 1.0], "NMO, Q x0.8 and x1.2")
for s in ("events_nmo", "thin_nmo", "dense_nmo"):
    show([c for c in nmo if c[0] == s], s); show([c for c in nmo if c[0] == s and c[3] == 1.0], s + ", Q correct")
for s in ("events_aligned", "dense_aligned"):
    show([c for c in ali if c[0] == s], s + " (plain methods, Q correct): reference")
print("By model, NMO, E_G median split_st / direct_st, fraction better:")
for mo in ("M1", "M2", "M3"):
    s = [c for c in nmo if c[1] == mo]; print(f"  {mo}: {np.median([T[c]['split_st']['eg'] for c in s]):.3f} / {np.median([T[c]['direct_st']['eg'] for c in s]):.3f}, {np.mean([T[c]['split_st']['eg'] < T[c]['direct_st']['eg'] for c in s]):.2f}")
print("By scene, fraction split_st better than direct_st on E_G:", {s: round(float(np.mean([T[c]['split_st']['eg'] < T[c]['direct_st']['eg'] for c in nmo if c[0] == s])), 2) for s in ("events_nmo", "thin_nmo", "dense_nmo")})
print("Chosen lambda, M1, SNR 18, Q correct:", {m: T[("dense_nmo", "M1", 18, 1.0)][m]["lam"] for m in METHODS})
p = OUT / "gate1b_stress.json"
if p.exists():
    S = json.loads(p.read_text(encoding="utf-8")); TS = table(S, ("scene", "stress")); print("\nMismatch stress group, M1, Q correct, SNR 18 (E_G / template G/A / wave)")
    for c in sorted(TS):
        print(f"{c[0]:11s} {c[1]:6s} " + " | ".join(f"{m} {TS[c][m]['eg']:.3f}/{TS[c][m].get('tga', float('nan')):.3f}/{TS[c][m]['wave']:.3f}" for m in ("split_st", "direct_st", "split", "direct")))
(OUT / "gate1b_summary.json").write_text(json.dumps(dict(R1=bool(R1), R2=bool(R2), template_ratio=tr, corr_gain=dc, frac=fr, ratio=ra, diff=di, cells={"|".join(map(str, c)): T[c] for c in T}), indent=1), encoding="utf-8")
