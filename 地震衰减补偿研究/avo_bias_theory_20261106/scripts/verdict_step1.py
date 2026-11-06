"""Verdicts of STEP1_ZH.md on results/step1.json -> results/step1_verdict.json"""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"; R = json.loads((OUT / "step1.json").read_text(encoding="utf-8")); C = sorted({(r["M"], r["fp"], r["snr"]) for r in R})


def v(q, m, k, cells=C):
    return np.array([next(r[k] for r in R if (r["M"], r["fp"], r["snr"]) == c and r["qf"] == q and r["m"] == m) for c in cells])


bm, bs, bd = v(1.0, "method", "bias"), v(1.0, "split", "bias"), v(1.0, "direct", "bias"); wm, ws, wd = v(1.0, "method", "wave"), v(1.0, "split", "wave"), v(1.0, "direct", "wave")
S1 = bm.mean() <= 1.2 * bs.mean(); S1b = bm.mean() <= 0.25 * bd.mean(); S2 = (ws - wm).mean() >= 0.03 and np.mean(wm < ws) >= 0.75; S3 = wm.mean() <= wd.mean() + 0.02
Wm = np.r_[v(0.8, "method", "wave"), v(1.2, "method", "wave")]; Ws = np.r_[v(0.8, "split", "wave"), v(1.2, "split", "wave")]; Bm = np.r_[v(0.8, "method", "bias"), v(1.2, "method", "bias")]; Bs = np.r_[v(0.8, "split", "bias"), v(1.2, "split", "bias")]; S4 = Wm.mean() <= Ws.mean() and Bm.mean() <= 1.1 * Bs.mean() and not np.any(Wm > 1.5 * Ws)
print(f"S1 |bias| method {bm.mean():.3f}, split {bs.mean():.3f}, ratio {bm.mean() / bs.mean():.3f} (<= 1.2) -> {'pass' if S1 else 'FAIL'}")
print(f"S1b direct {bd.mean():.3f}, ratio {bm.mean() / bd.mean():.3f} (<= 0.25) -> {'pass' if S1b else 'FAIL'}")
print(f"S2 wave method {wm.mean():.3f}, split {ws.mean():.3f}, difference {(ws - wm).mean():+.3f} (>= 0.03), better in {np.mean(wm < ws):.0%} (>= 75%) -> {'pass' if S2 else 'FAIL'}")
print(f"S3 wave direct {wd.mean():.3f}, method - direct {wm.mean() - wd.mean():+.3f} (<= +0.02) -> {'pass' if S3 else 'FAIL'}")
print(f"S4 Q wrong: wave method {Wm.mean():.3f}, split {Ws.mean():.3f}; |bias| {Bm.mean():.3f}, split {Bs.mean():.3f} (ratio {Bm.mean() / Bs.mean():.3f} <= 1.1); worst wave ratio {np.max(Wm / Ws):.2f} (<= 1.5) -> {'pass' if S4 else 'FAIL'}")
print("OVERALL:", "PASS" if S1 and S1b and S2 and S3 and S4 else "FAIL")
print("\nCell details: M fp snr | Q correct method/split/direct |bias| ; wave")
for c, a, b, d, e, f, g in zip(C, bm, bs, bd, wm, ws, wd):
    print(f"{c[0]} {c[1]:.0f} {c[2]:2d} | {a:.3f} {b:.3f} {d:.3f} ; {e:.3f} {f:.3f} {g:.3f}")
for q in (0.8, 1.2):
    print(f"\nQ x{q} mean: |bias| method {v(q, 'method', 'bias').mean():.3f} split {v(q, 'split', 'bias').mean():.3f} direct {v(q, 'direct', 'bias').mean():.3f}; wave method {v(q, 'method', 'wave').mean():.3f} split {v(q, 'split', 'wave').mean():.3f} direct {v(q, 'direct', 'wave').mean():.3f}")
for lab, sel in (("fp 60", lambda c: c[1] == 60), ("fp 40", lambda c: c[1] == 40), ("SNR 18", lambda c: c[2] == 18), ("SNR 12", lambda c: c[2] == 12)):
    cs = [c for c in C if sel(c)]; print(f"{lab}, Q correct: |bias| {v(1.0, 'method', 'bias', cs).mean():.3f} {v(1.0, 'split', 'bias', cs).mean():.3f} {v(1.0, 'direct', 'bias', cs).mean():.3f}; wave {v(1.0, 'method', 'wave', cs).mean():.3f} {v(1.0, 'split', 'wave', cs).mean():.3f} {v(1.0, 'direct', 'wave', cs).mean():.3f}")
print("worst wave cells, Q wrong:", sorted(((round(float(r["wave"]), 2), r["M"], r["fp"], r["snr"], r["qf"]) for r in R if r["m"] == "method" and r["qf"] != 1.0), reverse=True)[:4], "split there:", [round(next(s["wave"] for s in R if s["m"] == "split" and (s["M"], s["fp"], s["snr"], s["qf"]) == c), 2) for c in [(r["M"], r["fp"], r["snr"], r["qf"]) for r in sorted((r for r in R if r["m"] == "method" and r["qf"] != 1.0), key=lambda r: -r["wave"])[:4]]])
(OUT / "step1_verdict.json").write_text(json.dumps(dict(S1=bool(S1), S1b=bool(S1b), S2=bool(S2), S3=bool(S3), S4=bool(S4), bias=[float(bm.mean()), float(bs.mean()), float(bd.mean())], wave=[float(wm.mean()), float(ws.mean()), float(wd.mean())]), indent=1), encoding="utf-8")
