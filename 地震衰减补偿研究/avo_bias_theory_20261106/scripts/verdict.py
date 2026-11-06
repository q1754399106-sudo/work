"""Verdicts of THEORY_ZH.md section 5 on results/test.json, test_disc.json, test_spec.json -> results/verdict.json"""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"
R = json.loads((OUT / "test.json").read_text(encoding="utf-8")); D = json.loads((OUT / "test_disc.json").read_text(encoding="utf-8")); P = json.loads((OUT / "test_spec.json").read_text(encoding="utf-8")); res = {}
bias = lambda r: r["meas"] - r["truth"]
t1 = []
for snr in (60, 18):
    for m in ("direct", "split"):
        for lab, sel in (("Q correct", lambda r: r["qf"] == 1.0), ("Q wrong", lambda r: r["qf"] != 1.0)):
            S = [r for r in R if r["snr"] == snr and r["m"] == m and sel(r)]; d = np.array([bias(r) for r in S]); p = np.array([r["pred"] - r["truth"] for r in S]); b = float(np.polyfit(p, d, 1)[0]); r2 = float(np.corrcoef(p, d)[0, 1] ** 2); ok = 0.9 <= b <= 1.1 and r2 >= 0.95; t1.append(ok)
            print(f"T1 SNR {snr} {m:6s} {lab:9s} n {len(S):3d}: slope {b:.3f} R2 {r2:.3f} rmse {np.sqrt(np.mean((d - p) ** 2)):.3f} {'ok' if ok else 'FAIL'}")
T1 = all(t1)
Q = [r for r in R if r["qf"] == 1.0]; key = lambda r: (r["M"], r["snr"], r["lam"], r["j"]); by = {}
for r in Q:
    by.setdefault(key(r), {})[r["m"]] = r
fs = np.array([abs(bias(v["fc_star"])) for v in by.values()]); sp = np.array([abs(bias(v["split"])) for v in by.values()]); di = np.array([abs(bias(v["direct"])) for v in by.values()])
frac = float(np.mean(fs < sp)); r_s = float(np.median(fs / sp)); r_d = float(np.median(fs / di)); T2 = frac >= 0.8 and r_s <= 0.5; T3 = r_d <= 0.25
W = {m: float(np.mean([r["wave"] for r in Q if r["m"] == m and r["snr"] == 18 and r["lam"] in (1e-3, 1e-2) and r["j"] == 0])) for m in ("fc_star", "split", "direct")}; T2b = W["fc_star"] - W["split"] <= 0.02
print(f"\nT2 |bias| fc* < split in {frac:.1%} of {len(fs)} cases (>= 80%), median ratio {r_s:.3f} (<= 0.5) -> {'pass' if T2 else 'FAIL'}"); print(f"T2b waveform at 18 dB, lambda 1e-3 and 1e-2: fc* {W['fc_star']:.3f}, split {W['split']:.3f}, direct {W['direct']:.3f}; fc* - split {W['fc_star'] - W['split']:+.3f} (<= +0.02) -> {'pass' if T2b else 'FAIL'}"); print(f"T3 median |bias| ratio fc*/direct {r_d:.3f} (<= 0.25) -> {'pass' if T3 else 'FAIL'}")
Dq = [r for r in D if r["qf"] == 1.0]; mb = {m: float(np.mean([abs(bias(r)) for r in Dq if r["m"] == m])) for m in ("fc_star", "split", "direct")}; mw = {m: float(np.mean([r["wave"] for r in Dq if r["m"] == m and r["j"] == 0])) for m in mb}
T4 = mb["fc_star"] <= 0.5 * mb["split"] and mb["fc_star"] <= 0.25 * mb["direct"] and mw["fc_star"] <= mw["split"] + 0.02
print(f"T4 discrepancy lambda, 18 dB, Q correct: mean |bias| fc* {mb['fc_star']:.3f}, split {mb['split']:.3f}, direct {mb['direct']:.3f}; waveform {mw['fc_star']:.3f}, {mw['split']:.3f}, {mw['direct']:.3f}; chosen lambda medians", {m: round(float(np.median([r['lam_med'] for r in Dq if r['m'] == m])), 7) for m in mb}, "->", 'pass' if T4 else 'FAIL')
print("\nOVERALL:", "PASS" if T1 and T2 and T2b and T3 and T4 else "FAIL", {"T1": T1, "T2": T2, "T2b": T2b, "T3": T3, "T4": T4})
print("\nDescriptive")
for snr in (60, 18):
    for qf in (0.8, 1.0, 1.2):
        print(f" SNR {snr} Q x{qf}: mean |bias| by method, lambda 1e-3/1e-2: " + "; ".join(f"{m} {np.mean([abs(bias(r)) for r in R if r['m'] == m and r['snr'] == snr and r['qf'] == qf and r['lam'] == 1e-3]):.3f}/{np.mean([abs(bias(r)) for r in R if r['m'] == m and r['snr'] == snr and r['qf'] == qf and r['lam'] == 1e-2]):.3f}" for m in ("fc_star", "split", "direct")))
print(" discrepancy lambda, Q wrong: " + "; ".join(f"Q x{qf}: " + ", ".join(f"{m} |bias| {np.mean([abs(bias(r)) for r in D if r['m'] == m and r['qf'] == qf]):.3f} wave {np.mean([r['wave'] for r in D if r['m'] == m and r['qf'] == qf and r['j'] == 0]):.3f}" for m in ("fc_star", "split", "direct")) for qf in (0.8, 1.2)))
print(" per model, 18 dB, lambda 1e-3, Q correct:", {M: {m: round(float(np.mean([abs(bias(r)) for r in Q if r['M'] == M and r['m'] == m and r['snr'] == 18 and r['lam'] == 1e-3])), 3) for m in ('fc_star', 'split', 'direct')} for M in ('M4', 'M5')})
print(" wrong wavelet in fc*: mean |bias| by assumed peak x0.8, x1, x1.2, lambda 1e-3/1e-2:", {fp: [round(float(np.mean([abs(r['bias']) for r in P if r['fp'] == fp and r['lam'] == l])), 3) for l in (1e-3, 1e-2)] for fp in (0.8, 1.0, 1.2)}, "mean fc", {fp: round(float(np.mean([r['fc'] for r in P if r['fp'] == fp])), 1) for fp in (0.8, 1.0, 1.2)})
(OUT / "verdict.json").write_text(json.dumps(dict(T1=bool(T1), T2=bool(T2), T2b=bool(T2b), T3=bool(T3), T4=bool(T4), frac=frac, ratio_split=r_s, ratio_direct=r_d, wave=W, disc_bias=mb, disc_wave=mw), indent=1), encoding="utf-8")
