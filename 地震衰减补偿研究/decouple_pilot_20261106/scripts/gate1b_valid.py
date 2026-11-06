"""Validity check of gate 1b (post hoc): with the exact operator and almost no noise, split_st and direct_st should recover the gather.
If they do not, the 1b verdict may reflect a fault in the implementation rather than the method. M1, Q correct, 3 seeds.
Prints for each SNR and small lambda: waveform error, E_G, corr G, template G/A per event, including a much finer lambda range."""
import numpy as np
import gate1 as g
from gate1 import C, GEO, NX, n, T, I, DT, MODELS, LX, avo, bandpass, synth, reflectors, a0_ax, atten
from gate1b import operators, ridge_apply, score, nmo_cols, EV
from gate1b_diag import template_ev, K13

Q = MODELS["M1"]; op = operators(Q); G0 = atten(op["a0"]); sig = op["sig"]; keep = K13
print("largest/smallest singular value, F_s far offset:", np.round(np.linalg.svd(op["Fs"][-1], compute_uv=False)[[0, 200, 400, 600]], 8).tolist(), "; F0/sig:", np.round(np.linalg.svd(op["F0"] / sig, compute_uv=False)[[0, 200, 400, 600]], 8).tolist())
for scene in ("events_x", "dense_x"):
    for snr in (120, 60, 18):
        D, Tr = [], []
        for seed in range(3):
            rng = np.random.default_rng(555 + seed); idx, R = reflectors(scene, rng); dp, tru = synth(idx, R, Q, True); rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); d = dp + rms / 10 ** (snr / 20) * rng.standard_normal(dp.shape)
            D.append(np.stack([nmo_cols(d[:, [i]], GEO[i]["tx"], keep[:, i])[:, 0] for i in range(NX)], axis=1)); Tr.append(tru)
        Dat = np.stack(D, axis=2); Tru = np.stack(Tr, axis=2); At, Bt, ok = avo(bandpass(Tru), keep); print(f"\n{scene[:-2]} SNR {snr}")
        for lam in (1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3):
            F0inv = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); out = {}
            Rh = [ridge_apply(op["Fs"][i], lam, Dat[:, i, :] / sig) for i in range(NX)]
            out["direct_st"] = np.stack([op["W"] @ Rh[i] for i in range(NX)], axis=1)
            Rx = [ridge_apply(op["Fs"][i], LX, Dat[:, i, :] / sig) for i in range(NX)]
            out["split_st"] = np.stack([F0inv @ (op["F0"] @ Rx[i]) for i in range(NX)], axis=1)
            Rl = [ridge_apply(op["Fs"][i], lam, Dat[:, i, :] / sig) for i in range(NX)]
            out["split_st_lx=lam"] = np.stack([F0inv @ (op["F0"] @ Rl[i]) for i in range(NX)], axis=1)
            line = f"lam {lam:7.0e}"
            for m, O in out.items():
                s = score(O, Tru, keep, At, Bt, ok, False); te = np.nanmean([np.abs(template_ev(O[:, :, c], keep)[:, 1] / template_ev(O[:, :, c], keep)[:, 0] - template_ev(Tru[:, :, c], keep)[:, 1] / template_ev(Tru[:, :, c], keep)[:, 0]) for c in range(O.shape[2])], axis=0)
                line += f" | {m} w {s['wave'].mean():.3f} E {s['eg'].mean():.3f} c {s['cg'].mean():.3f} [{te[0]:.2f} {te[1]:.2f} {te[2]:.2f}]"
            print(line, flush=True)
