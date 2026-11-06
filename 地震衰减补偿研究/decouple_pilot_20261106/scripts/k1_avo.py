"""Pilot K1 (exploratory, no pre-registered rule): does the split protect AVO better than one-step inversion?
Mechanism tested: one-step damping at one lambda acts harder on far traces, because they are more attenuated, which can tilt
the amplitude against offset; the split equalises the offsets first and then damps every trace with the same operator.
Same gather model as the original method. Output scored against truth per event: intercept and gradient of template amplitude
against sin^2 theta, gradient/intercept ratio, far/near ratio, and the whole-gather waveform error."""
import json, sys
from pathlib import Path
import numpy as np
from core import C, ray, atten, ricker, reference

OUT = Path(__file__).resolve().parents[1] / "results"; OUT.mkdir(exist_ok=True)
LAMS = np.logspace(-6, 0, 19); LX = 1e-3; SEEDS = range(20); SNRS = (18, 12); QF = (0.8, 1.0, 1.2)
OFF = C.offset; n = C.nt


def svd_op(G):
    U, s, Vt = np.linalg.svd(G, full_matrices=False); return U, s, Vt


def apply(op, y, lam):
    U, s, Vt = op; return Vt.T @ ((s / (s ** 2 + lam))[:, None] * (U.T @ y))


def amps(D):
    """Template amplitude per event and offset; for regularised outputs the template is the output's own zero-offset
    wavelet shape is not assumed: a least-squares scale on the true Ricker, as event_metrics.m does."""
    out = np.zeros((3, D.shape[1]))
    for j, te in enumerate(C.event_t):
        k = np.abs(C.t - te) <= 0.07; w = ricker(C.t[k] - te); out[j] = (w @ D[k]) / (w @ w)
    return out


def avo_fit(a, x2):
    return np.array([np.polyfit(x2[j], a[j], 1)[::-1] for j in range(3)])     # intercept, gradient per event


truth = np.stack([reference(x)[0] for x in OFF], axis=1); SIN2 = np.stack([reference(x)[2] for x in OFF], axis=1)
aT = amps(truth); fT = avo_fit(aT, SIN2)
Gtrue = [atten(np.sum(ray(C.t, x), axis=0)[:] if False else ray(C.t, x)[0] + ray(C.t, x)[1]) for x in OFF]
clean = np.stack([Gtrue[i] @ truth[:, i] for i in range(len(OFF))], axis=1)
ops = {}
for qf in QF:
    Q = C.Q * qf; a0 = ray(C.t, 0, Q)[0]; G0 = atten(a0); op0 = svd_op(G0); ox, od = [], []
    for x in OFF:
        ax = ray(C.t, x, Q)[1]; ox.append(svd_op(atten(ax))); od.append(svd_op(atten(a0 + ax)))
    ops[qf] = dict(op0=op0, ox=ox, od=od); print("operators ready", qf, flush=True)


def metrics(D):
    a = amps(D); f = avo_fit(a, SIN2)
    ratio = f[:, 1] / f[:, 0]; ratioT = fT[:, 1] / fT[:, 0]
    return dict(wave=float(np.linalg.norm(D - truth) / np.linalg.norm(truth)),
                grad_ratio_err=float(np.mean(np.abs(ratio - ratioT))),
                gradient_rel=float(np.mean(np.abs(f[:, 1] - fT[:, 1]) / np.abs(fT[:, 1]))),
                intercept_rel=float(np.mean(np.abs(f[:, 0] - fT[:, 0]) / np.abs(fT[:, 0]))),
                far_near_err=float(np.mean(np.abs(a[:, -1] / a[:, 0] - aT[:, -1] / aT[:, 0]))))


rows = []
for snr in SNRS:
    sig = np.linalg.norm(clean, axis=0) / np.sqrt(n) / 10 ** (snr / 20)
    for seed in SEEDS:
        d = clean + sig * np.random.default_rng(seed).standard_normal(clean.shape)
        for qf in QF:
            o = ops[qf]
            y = np.stack([apply(o["ox"][i], d[:, [i]], LX)[:, 0] for i in range(len(OFF))], axis=1)   # step 1
            G0y = y                                                                               # step 1 output
            for lam in LAMS:
                split = apply(o["op0"], y, lam)
                direct = np.stack([apply(o["od"][i], d[:, [i]], lam)[:, 0] for i in range(len(OFF))], axis=1)
                rows.append(dict(snr=snr, seed=seed, qf=qf, lam=float(lam), method="split", **metrics(split)))
                rows.append(dict(snr=snr, seed=seed, qf=qf, lam=float(lam), method="direct", **metrics(direct)))
    print("snr", snr, "done", flush=True)
(OUT / "k1_rows.json").write_text(json.dumps(rows), encoding="utf-8")


def sel(**kw):
    return [r for r in rows if all(r[k] == v for k, v in kw.items())]


print("\nMean over 20 seeds. Each method at its own lambda minimising the mean waveform error (truth used), so both are at their best waveform.")
print(f"{'SNR':>3} {'Q':>4} {'method':>7} {'lambda':>8} {'wave':>6} {'G/A err':>8} {'grad rel':>9} {'int rel':>8} {'far/near':>9}")
summary = []
for snr in SNRS:
    for qf in QF:
        best = {}
        for m in ("split", "direct"):
            per = {l: np.mean([r["wave"] for r in sel(snr=snr, qf=qf, method=m, lam=float(l))]) for l in LAMS}; l = min(per, key=per.get); best[m] = float(l)
            R = sel(snr=snr, qf=qf, method=m, lam=float(l)); v = {k: float(np.mean([r[k] for r in R])) for k in ("wave", "grad_ratio_err", "gradient_rel", "intercept_rel", "far_near_err")}
            print(f"{snr:3d} {qf:4.1f} {m:>7} {l:8.1e} {v['wave']:6.3f} {v['grad_ratio_err']:8.3f} {v['gradient_rel']:9.3f} {v['intercept_rel']:8.3f} {v['far_near_err']:9.3f}", flush=True)
            summary.append(dict(snr=snr, qf=qf, method=m, lam=float(l), **v))
        A = sel(snr=snr, qf=qf, method="split", lam=best["split"]); B = sel(snr=snr, qf=qf, method="direct", lam=best["direct"])
        d = [a["grad_ratio_err"] - b["grad_ratio_err"] for a, b in zip(A, B)]
        print(f"     paired G/A error split - direct: mean {np.mean(d):+.3f}, split better in {int(np.sum(np.array(d) < 0))}/{len(d)} seeds; waveform split - direct {np.mean([a['wave'] - b['wave'] for a, b in zip(A, B)]):+.3f}", flush=True)
print("\nSame lambda for both, SNR 18, Q correct: lambda, wave split/direct, G/A error split/direct, far/near split/direct")
for l in LAMS[::3]:
    a = sel(snr=18, qf=1.0, method="split", lam=float(l)); b = sel(snr=18, qf=1.0, method="direct", lam=float(l)); g = lambda R, k: np.mean([r[k] for r in R])
    print(f"{l:8.1e}  {g(a, 'wave'):.3f}/{g(b, 'wave'):.3f}  {g(a, 'grad_ratio_err'):.3f}/{g(b, 'grad_ratio_err'):.3f}  {g(a, 'far_near_err'):.3f}/{g(b, 'far_near_err'):.3f}", flush=True)
(OUT / "k1_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
