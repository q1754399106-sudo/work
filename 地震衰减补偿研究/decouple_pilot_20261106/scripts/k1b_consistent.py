"""Pilot K1b (exploratory): can one-step inversion get the split's AVO with its own waveform, by making the damping
consistent across offsets?  One-step ridge: (G'G + lam I)^-1 G'. A far trace is more attenuated, so at the same lam its
high frequencies are cut earlier. Offset-consistent damping: lam w(t, x), w = exp(-2 pi fc ax(t, x)), which moves each trace's
cut-off at frequency fc back to where the zero-offset trace has it. Also checks the repository's earlier variants:
one-step and covariance split with lam per trace from the discrepancy principle (known noise sigma), which is what the
earlier "split and one-step are nearly the same" comparison used."""
import json
from pathlib import Path
import numpy as np
from core import C, ray, atten, ricker, reference

OUT = Path(__file__).resolve().parents[1] / "results"; n = C.nt; OFF = C.offset; I = np.eye(n)
LAMS = np.logspace(-6, 0, 13); LX = 1e-3; SEEDS = 12; SNRS = (18, 12); QF = (0.8, 1.0, 1.2); FC = (40.0, 60.0, 80.0); GRID = np.logspace(-7, 0, 41)
truth = np.stack([reference(x)[0] for x in OFF], axis=1); SIN2 = np.stack([reference(x)[2] for x in OFF], axis=1)
Gtrue = [atten(ray(C.t, x)[0] + ray(C.t, x)[1]) for x in OFF]; clean = np.stack([Gtrue[i] @ truth[:, i] for i in range(len(OFF))], axis=1)


def amps(D):
    out = np.zeros((3, D.shape[1]))
    for j, te in enumerate(C.event_t):
        k = np.abs(C.t - te) <= 0.07; w = ricker(C.t[k] - te); out[j] = (w @ D[k]) / (w @ w)
    return out


def fit(a):
    return np.array([np.polyfit(SIN2[j], a[j], 1)[::-1] for j in range(3)])


aT = amps(truth); fT = fit(aT)


def score(D):
    a = amps(D); f = fit(a)
    return dict(wave=float(np.linalg.norm(D - truth) / np.linalg.norm(truth)), ga=float(np.mean(np.abs(f[:, 1] / f[:, 0] - fT[:, 1] / fT[:, 0]))),
                fn=float(np.mean(np.abs(a[:, -1] / a[:, 0] - aT[:, -1] / aT[:, 0]))), ia=float(np.mean(np.abs(f[:, 0] - fT[:, 0]) / np.abs(fT[:, 0]))))


def ridge(G, lam, w=None):
    d = lam * (np.ones(n) if w is None else w); return np.linalg.solve(G.T @ G + np.diag(d), G.T)


def discrepancy(G, y, target):
    """Largest lambda on GRID with residual <= target, as discrepancy_solve (per trace; SVD for speed)."""
    U, s, Vt = np.linalg.svd(G, full_matrices=False); b = U.T @ y; res0 = np.sum(y ** 2) - np.sum(b ** 2); best = None
    for lam in GRID[::-1]:
        r = np.sum((lam / (s ** 2 + lam) * b) ** 2) + res0
        if r <= target ** 2:
            best = lam; break
    best = GRID[0] if best is None else best
    return Vt.T @ (s / (s ** 2 + best) * b)


rows = []
for qf in QF:
    Q = C.Q * qf; a0 = ray(C.t, 0, Q)[0]; G0 = atten(a0); Ax = [ray(C.t, x, Q)[1] for x in OFF]; Gx = [atten(a) for a in Ax]; Gd = [atten(a0 + a) for a in Ax]
    Bx = [ridge(g, LX) for g in Gx]; ops = {}
    for lam in LAMS:
        ops[("split", lam)] = [ridge(G0, lam) @ b for b in Bx]
        ops[("direct", lam)] = [ridge(g, lam) for g in Gd]
        for fc in FC:
            ops[(f"consistent{fc:.0f}", lam)] = [ridge(g, lam, np.exp(-2 * np.pi * fc * a)) for g, a in zip(Gd, Ax)]
    print("operators", qf, flush=True)
    Ls = [np.linalg.cholesky(b @ b.T + 1e-6 * I) for b in Bx]
    for snr in SNRS:
        sig = np.linalg.norm(clean, axis=0) / np.sqrt(n) / 10 ** (snr / 20)
        for seed in range(SEEDS):
            d = clean + sig * np.random.default_rng(seed).standard_normal(clean.shape)
            for (m, lam), B in ops.items():
                D = np.stack([B[i] @ d[:, i] for i in range(len(OFF))], axis=1); rows.append(dict(qf=qf, snr=snr, seed=seed, m=m, lam=float(lam), **score(D)))
            D1 = np.stack([discrepancy(Gd[i], d[:, i], sig[i] * np.sqrt(n)) for i in range(len(OFF))], axis=1)
            D2 = np.stack([discrepancy(np.linalg.solve(Ls[i], G0), np.linalg.solve(Ls[i], Bx[i] @ d[:, i]), sig[i] * np.sqrt(n)) for i in range(len(OFF))], axis=1)
            rows.append(dict(qf=qf, snr=snr, seed=seed, m="direct_disc", lam=0.0, **score(D1))); rows.append(dict(qf=qf, snr=snr, seed=seed, m="covsplit_disc", lam=0.0, **score(D2)))
    print("done", qf, flush=True)
(OUT / "k1b_rows.json").write_text(json.dumps(rows), encoding="utf-8")
ms = ["split", "direct", "consistent40", "consistent60", "consistent80"]; K = ("wave", "ga", "fn", "ia"); summary = []
print("\nEach method at its own best mean waveform lambda (truth used). wave = waveform error, G/A = gradient/intercept error, F/N = far/near error, I = intercept rel. error")
for snr in SNRS:
    for qf in QF:
        print(f"--- SNR {snr}, Q x {qf}")
        for m in ms + ["direct_disc", "covsplit_disc"]:
            R = [r for r in rows if r["qf"] == qf and r["snr"] == snr and r["m"] == m]; lams = sorted({r["lam"] for r in R})
            per = {l: np.mean([r["wave"] for r in R if r["lam"] == l]) for l in lams}; l = min(per, key=per.get); S = [r for r in R if r["lam"] == l]; v = {k: float(np.mean([r[k] for r in S])) for k in K}
            print(f"{m:15s} lam {l:7.1e} wave {v['wave']:.3f} G/A {v['ga']:.3f} F/N {v['fn']:.3f} I {v['ia']:.3f}", flush=True); summary.append(dict(snr=snr, qf=qf, m=m, lam=l, **v))
print("\nTrade-off at SNR 18, Q correct: best G/A among lambdas with waveform within 0.02 of each method's best")
for m in ms:
    R = [r for r in rows if r["qf"] == 1.0 and r["snr"] == 18 and r["m"] == m]; per = {l: (np.mean([r["wave"] for r in R if r["lam"] == l]), np.mean([r["ga"] for r in R if r["lam"] == l])) for l in sorted({r["lam"] for r in R})}
    b = min(v[0] for v in per.values()); l, (w, g) = min(((l, v) for l, v in per.items() if v[0] <= b + 0.02), key=lambda z: z[1][1]); print(f"{m:15s} lam {l:7.1e} wave {w:.3f} G/A {g:.3f}")
print("\nFull lambda curves at SNR 18, Q correct (wave / G/A)")
for lam in LAMS:
    print(f"{lam:8.1e} " + " ".join(f"{m[:6]} {np.mean([r['wave'] for r in rows if r['qf'] == 1.0 and r['snr'] == 18 and r['m'] == m and abs(r['lam'] - lam) < 1e-14]):.3f}/{np.mean([r['ga'] for r in rows if r['qf'] == 1.0 and r['snr'] == 18 and r['m'] == m and abs(r['lam'] - lam) < 1e-14]):.3f}" for m in ms))
(OUT / "k1b_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
