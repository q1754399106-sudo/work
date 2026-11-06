"""Post hoc diagnostic (does not change any verdict). Was my gate-1 reading, "uncorrected stretch is the main cause", right, or is
it the stretch MUTE, which removes far offsets at shallow times and leaves an ill-conditioned gradient fit?
Same forward, M1-M3, Q correct, 10 seeds. Configurations:
  aligned_full   no stretch, all offsets
  aligned_mute   no stretch, the stretch mute's mask applied to data, truth and fit
  nmo_mute       real NMO stretch, mute at 1.3 (gate 1)
  nmo_nomute     real NMO stretch, mute at 2.5 (stretch kept on far offsets)
Each at noise-free (1e4 dB-like: SNR 60) and 18 dB. Reports E_G, correlation of gradient, template G/A per event (event 3 at 0.90 s
keeps nearly all offsets), kept offsets per event."""
import numpy as np
import gate1 as g
from gate1 import C, GEO, SIN2, NX, n, T, I, DT, NSYN, MODELS, LAMS, LX, a0_ax, atten, ricker, avo, bandpass, solve_all, synth, reflectors
from gate1b import score, EV, nmo_cols, txs


def keep_at(mute):
    return np.stack([GEO[i]["tx"] <= mute * np.maximum(T, DT) for i in range(NX)], axis=1)


K13 = keep_at(1.3); K25 = keep_at(2.5); KALL = np.ones((n, NX), bool)
print("kept offsets at the three events, mute 1.3:", K13[EV].sum(axis=1).tolist(), "of", NX, "; sin^2 max kept:", [round(float(SIN2[EV[j]][K13[EV[j]]].max()), 3) for j in range(3)], "; full:", [round(float(SIN2[EV[j]].max()), 3) for j in range(3)])


def template_ev(D, keep):
    out = []
    for j, te in enumerate(C.event_t):
        k = np.abs(T - te) <= 0.07; w = ricker(T[k] - te); a = (w @ D[k]) / (w @ w); ok = keep[EV[j]]
        out.append(np.polyfit(SIN2[EV[j]][ok], a[ok], 1)[::-1] if ok.sum() >= 3 else np.array([np.nan, np.nan]))
    return np.array(out)


def run(Mname, scene, cfg, snr):
    Q = MODELS[Mname]; a0, ax = a0_ax(Q); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]
    Bx = [np.linalg.solve(G.T @ G + LX * I, G.T) for G in Gx]; nmo = cfg.startswith("nmo"); keep = dict(aligned_full=KALL, aligned_mute=K13, nmo_mute=K13, nmo_nomute=K25)[cfg]; D, Tr = [], []
    for seed in range(10):
        rng = np.random.default_rng(555 + seed); idx, R = reflectors(scene, rng); dp, tru = synth(idx, R, Q, nmo); rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); d = dp + rms / 10 ** (snr / 20) * rng.standard_normal(dp.shape)
        y = np.stack([nmo_cols(d[:, [i]], GEO[i]["tx"], keep[:, i])[:, 0] for i in range(NX)], axis=1) if nmo else d[:n] * keep
        D.append(y); Tr.append(tru)
    Dat = np.stack(D, axis=2); Tru = np.stack(Tr, axis=2); At, Bt, ok = avo(bandpass(Tru), keep); Y1 = np.stack([Bx[i] @ Dat[:, i, :] for i in range(NX)], axis=1); res = {}
    for m in ("split", "direct"):
        best = None
        for lam in LAMS:
            if m == "split":
                F = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); O = np.stack([F @ Y1[:, i, :] for i in range(NX)], axis=1)
            else:
                O = np.stack([solve_all(Gd[i], np.ones(n), lam, Dat[:, i, :]) for i in range(NX)], axis=1)
            s = score(O, Tru, keep, At, Bt, ok, False)
            if best is None or s["wave"].mean() < best[0]:
                best = (s["wave"].mean(), lam, O, s)
        w, lam, O, s = best; te = []
        for c in range(O.shape[2]):
            fT = template_ev(Tru[:, :, c], keep); fE = template_ev(O[:, :, c], keep); te.append(np.abs(fE[:, 1] / fE[:, 0] - fT[:, 1] / fT[:, 0]))
        res[m] = dict(eg=float(s["eg"].mean()), cg=float(s["cg"].mean()), wave=float(w), te=np.nanmean(te, axis=0))
    return res


if __name__ == "__main__":
    for scene in ("events_x", "dense_x"):
        for Mname in ("M1", "M2", "M3"):
            print(f"\n{scene[:-2]} {Mname}: config, snr | split E_G corr wave [template G/A events 1,2,3] | direct same")
            for cfg in ("aligned_full", "aligned_mute", "nmo_mute", "nmo_nomute"):
                for snr in (60, 18):
                    r = run(Mname, scene, cfg, snr); a, b = r["split"], r["direct"]; ta = " ".join(f"{v:.2f}" for v in a["te"]); tb = " ".join(f"{v:.2f}" for v in b["te"])
                    print(f"{cfg:13s} {snr:2d} | {a['eg']:6.3f} {a['cg']:6.3f} {a['wave']:.3f} [{ta}] | {b['eg']:6.3f} {b['cg']:6.3f} {b['wave']:.3f} [{tb}]", flush=True)
