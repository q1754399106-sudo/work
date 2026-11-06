"""Post hoc diagnostic of gate 1 (does not change the verdict). K1 found split G/A error 0.04 against direct 0.38 with the
template metric; gate 1 finds M_G 0.62 against 0.59 on events_aligned. Is that the metric, or the forward model?
Same data and lambda choice as gate 1 (own best waveform), Q correct, SNR 18: template G/A per event, with sin^2 from the
event, next to M_G, corr G and the gains. Also scene orig_like: Shuey events through the gate-1 forward."""
import numpy as np
import gate1 as g
from gate1 import C, GEO, SIN2, NX, n, T, I, atten, ricker, avo, bandpass, score_batch, a0_ax, MODELS, LAMS, LX, solve_all, synth, to_t0, reflectors
from core import reference

EV = np.round(C.event_t / g.DT).astype(int)


def template(D, keep):
    amp = np.zeros((3, NX))
    for j, te in enumerate(C.event_t):
        k = np.abs(T - te) <= 0.07; w = ricker(T[k] - te); amp[j] = (w @ D[k]) / (w @ w)
    fit = []
    for j in range(3):
        ok = keep[EV[j]]; fit.append(np.polyfit(SIN2[EV[j]][ok], amp[j][ok], 1)[::-1])
    return np.array(fit), amp


def run(Mname, kind, nmo, shuey=False, lam_fixed=None):
    Q = MODELS[Mname]; a0, ax = a0_ax(Q); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]
    Bx = [np.linalg.solve(G.T @ G + LX * I, G.T) for G in Gx]; cols = []
    for seed in range(10):
        rng = np.random.default_rng(777 + seed); idx, R = reflectors(f"{kind}_x", rng)
        if shuey:
            idx = EV.copy(); R = np.array([[reference(x)[1][j] for x in C.offset] for j in range(3)])
        dp, tru = synth(idx, R, Q, nmo); rms = np.sqrt(np.sum(dp[:int(2.4 / g.DT)] ** 2, axis=0) / n)
        y, keep = to_t0(dp + rms / 10 ** (18 / 20) * rng.standard_normal(dp.shape), nmo); cols.append((y, tru))
    Dat = np.stack([c[0] for c in cols], axis=2); Tru = np.stack([c[1] for c in cols], axis=2); At, Bt, ok = avo(bandpass(Tru), keep); Y1 = np.stack([Bx[i] @ Dat[:, i, :] for i in range(NX)], axis=1); out = {}
    for m in ("split", "direct"):
        best = None
        for lam in LAMS:
            if m == "split":
                F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); O = np.stack([F0 @ Y1[:, i, :] for i in range(NX)], axis=1)
            else:
                O = np.stack([solve_all(Gd[i], np.ones(n), lam, Dat[:, i, :]) for i in range(NX)], axis=1)
            s = score_batch(O, Tru, keep, At, Bt, ok); w = s["wave"].mean()
            if best is None or w < best[0]:
                best = (w, lam, O, s)
        w, lam, O, s = best; ta = []
        for c in range(O.shape[2]):
            fT, _ = template(Tru[:, :, c], keep); fE, _ = template(O[:, :, c], keep)
            ta.append(np.mean(np.abs(fE[:, 1] / fE[:, 0] - fT[:, 1] / fT[:, 0])))
        out[m] = dict(lam=lam, wave=w, tGA=float(np.mean(ta)), mg=float(s["mg"].mean()), cg=float(s["cg"].mean()), ba=float(s["ba"].mean()), bg=float(s["bg"].mean()))
    return out


print("post hoc: Q correct, SNR 18, 10 seeds, each method at its own best waveform")
print(f"{'case':34s} {'method':7s} {'lam':>8} {'wave':>6} {'tmpl G/A':>9} {'M_G':>6} {'corr G':>7} {'b_A':>6} {'b_G':>6}")
for label, a in (("M1 Shuey events, gate-1 forward", ("M1", "events", False, True)), ("M1 events Zoeppritz aligned", ("M1", "events", False)), ("M2 events Zoeppritz aligned", ("M2", "events", False)), ("M3 events Zoeppritz aligned", ("M3", "events", False)),
                 ("M1 Shuey events, NMO", ("M1", "events", True, True)), ("M1 events Zoeppritz NMO", ("M1", "events", True)), ("M2 events Zoeppritz NMO", ("M2", "events", True)), ("M3 events Zoeppritz NMO", ("M3", "events", True))):
    r = run(*a)
    for m, v in r.items():
        print(f"{label:34s} {m:7s} {v['lam']:8.1e} {v['wave']:6.3f} {v['tGA']:9.3f} {v['mg']:6.3f} {v['cg']:7.3f} {v['ba']:6.3f} {v['bg']:6.3f}", flush=True)
print("\ntrue gradient/intercept per event (Zoeppritz fit over offsets; Shuey):")
for j in range(3):
    s = np.array([g.rpp(j, np.array([GEO[i]["p"][EV[j]] * C.vp[j]]))[0] for i in range(NX)]); x = SIN2[EV[j]]; p = np.polyfit(x, s, 1); sh = np.array([reference(xx)[1][j] for xx in C.offset]); q = np.polyfit(x, sh, 1)
    print(f"event {j}: Zoeppritz A {p[1]:+.3f} G {p[0]:+.3f}; Shuey A {q[1]:+.3f} G {q[0]:+.3f}; sin2 range {x.min():.3f}-{x.max():.3f}")
