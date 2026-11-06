"""Development of the stabilised method on used models M1-M5 (exploratory). SNR 18, Q x0.8, 1, 1.2, 20 seeds.
Variants at fixed lambda and at a no-truth lambda: gamma = 1, 1.25, 1.5, 2, none. Lambda rules for the no-truth case:
 L1 gather residual of the method itself; L2 gather residual of plain one-step ridge (the lambda then used by the method)."""
import json, sys, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, GEO, SIN2, NX, n, T, I, DT, a0_ax, atten, MODELS, synth
import theory as th
from measure import reflector_series, template, ga, ridge
import stable

OUT = Path(__file__).resolve().parents[1] / "results"; LX = 1e-3; GAMMAS = (1.0, 1.25, 1.5, 2.0, None); GRID = np.logspace(-6, 0, 13)[::-1]; QF = (0.8, 1.0, 1.2); SEEDS = 20; SNR = 18
rows = []; t0 = time.monotonic()
for Mname in sys.argv[1:] or ("M1", "M2", "M3", "M4", "M5"):
    Qt = MODELS[Mname]; idx, R = reflector_series(); dp, tru = synth(idx, R, Qt, False); clean = dp[:n]; rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); sigma = rms / 10 ** (SNR / 20); target2 = float(np.sum(sigma ** 2) * n)
    rng = np.random.default_rng(55000 + 1000 * list(MODELS).index(Mname)); Dat = clean[:, :, None] + sigma[None, :, None] * rng.standard_normal((n, NX, SEEDS)); aT = template(tru); gT = ga(aT)
    for qf in QF:
        Qm = Qt * qf; a0, ax = a0_ax(Qm); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; GtG = [G.T @ G for G in Gd]; GtD = [Gd[i].T @ Dat[:, i, :] for i in range(NX)]
        Y = np.stack([np.linalg.solve(G.T @ G + LX * I, G.T @ Dat[:, i, :]) for i, G in enumerate(Gx)], axis=1); res = {}; cols = {}

        def rec(m, lam, O):
            cols.setdefault(m, {})[float(lam)] = O
            r = sum(np.sum((Gd[i] @ O[:, i, :] - Dat[:, i, :]) ** 2, axis=0) for i in range(NX)); res.setdefault(m, {})[float(lam)] = r

        for lam in GRID:
            fs = np.array([th.fc_star(a0[k], lam) for k in range(n)]); F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); rec("split", lam, np.stack([F0 @ Y[:, i, :] for i in range(NX)], axis=1))
            rec("direct", lam, np.stack([np.linalg.solve(GtG[i] + lam * I, GtD[i]) for i in range(NX)], axis=1)); w0 = np.exp(-2 * np.pi * fs[:, None] * ax)
            for gm in GAMMAS:
                w = w0 if gm is None else np.maximum(w0, stable.bound(a0, ax, lam, gm)); rec(f"g{gm}", lam, np.stack([np.linalg.solve(GtG[i] + lam * np.diag(w[:, i]), GtD[i]) for i in range(NX)], axis=1))
        methods = list(cols)

        def measure(O):
            A = template(O); return ga(A.mean(axis=2)), float(np.mean(np.sqrt(np.sum((O - tru[:, :, None]) ** 2, axis=(0, 1)) / np.sum(tru ** 2))))

        for m in methods:
            for lam, O in cols[m].items():
                meas, wave = measure(O); rows.append(dict(M=Mname, qf=qf, m=m, rule="fixed", lam=lam, wave=wave, bias=[float(v) for v in meas - gT]))
            lams = sorted(cols[m], reverse=True)
            for rule, source in (("L1", m), ("L2", "direct")):
                O = np.zeros_like(Dat); ch = []
                for s in range(SEEDS):
                    pick = next((l for l in sorted(cols[source], reverse=True) if res[source][l][s] <= target2), min(lams)); O[:, :, s] = cols[m][pick][:, :, s]; ch.append(pick)
                meas, wave = measure(O); rows.append(dict(M=Mname, qf=qf, m=m, rule=rule, lam=float(np.median(ch)), wave=wave, bias=[float(v) for v in meas - gT]))
        print(Mname, qf, f"{time.monotonic() - t0:.0f}s", flush=True)
(OUT / "dev_stable.json").write_text(json.dumps(rows), encoding="utf-8")
ms = ["split", "direct"] + [f"g{g_}" for g_ in GAMMAS]
for rule in ("L1", "L2"):
    print(f"\nRule {rule} (no truth): mean |bias| / mean wave / max wave over models")
    for qf in QF:
        print(f" Q x{qf}: " + " | ".join(f"{m} {np.mean([np.abs(r['bias']).mean() for r in rows if r['rule'] == rule and r['qf'] == qf and r['m'] == m]):.3f}/{np.mean([r['wave'] for r in rows if r['rule'] == rule and r['qf'] == qf and r['m'] == m]):.3f}/{max(r['wave'] for r in rows if r['rule'] == rule and r['qf'] == qf and r['m'] == m):.2f}" for m in ms))
print("\nFixed lambda, Q correct: mean |bias| / wave")
for lam in (1e-4, 1e-3, 1e-2):
    print(f" {lam:.0e}: " + " | ".join(f"{m} {np.mean([np.abs(r['bias']).mean() for r in rows if r['rule'] == 'fixed' and r['qf'] == 1.0 and r['m'] == m and abs(r['lam'] - lam) / lam < 1e-6]):.3f}/{np.mean([r['wave'] for r in rows if r['rule'] == 'fixed' and r['qf'] == 1.0 and r['m'] == m and abs(r['lam'] - lam) / lam < 1e-6]):.3f}" for m in ms))
print("\nPer model, Q correct, rule L1 and L2 (|bias|/wave/lam):")
for M in sorted({r["M"] for r in rows}):
    for rule in ("L1", "L2"):
        print(f" {M} {rule}: " + " | ".join(f"{m} {np.abs(next(r for r in rows if r['M'] == M and r['rule'] == rule and r['qf'] == 1.0 and r['m'] == m)['bias']).mean():.3f}/{next(r for r in rows if r['M'] == M and r['rule'] == rule and r['qf'] == 1.0 and r['m'] == m)['wave']:.3f}" for m in ms))
