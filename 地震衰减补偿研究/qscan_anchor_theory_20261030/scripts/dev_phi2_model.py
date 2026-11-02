"""Development (before PLAN_BOUND_ZH.md is frozen): candidate pre-inversion weights phi2(t) on the profiles already
measured (base, S, D, D2, M; anchor-window and truncated-window geometries). Development seed 9940 only.

phi2 candidates are built from the local spectra of the clean signal s(f, t) and of the noise n(f, t) of one generated
dataset (c = 1), in 0.2 s windows sliding by 0.1 s, interpolated to the analysis samples:
  U   uniform (the published formula)
  P   sum_f f^2 s                    (f^2-weighted signal power; tried before in 5h)
  F   sum_f f^2 s/(s+n)              (Fisher information of a log-spectrum slope, saturating in SNR)
  F2  sum_f f^2 [s/(s+n)]^2          (Fisher information with noise on both sides)
  B   sum_f s/(s+n)                  (effective bandwidth)
kappa = a0 sum(phi a) / sum(phi a^2) over the analysis window.
-> analysis/dev_phi2_model.json, printed table.
"""
import json
from env15 import np, lib, semi, gen

STRUCT, NOISE, MULT, SEED = 10227, [10149, 10273, 10355], 0.5, 9940
PROF = {"base": None, "S": [30.0, 200.0, 200.0, 200.0], "D": [250.0, 250.0, 40.0, 40.0],
        "D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0]}
# measured K (two/three anchors, lambda 1/2): profile, analysis window, anchor centre -> K
MEAS = [("base", (1.2, 3.0), 1.4, 0.66), ("base", (1.2, 3.0), 2.0, 0.98), ("base", (1.2, 3.0), 2.6, 1.12),
        ("base", (1.2, 2.2), 1.4, 0.81), ("S", (1.2, 3.0), 1.4, 0.856), ("D", (1.2, 3.0), 1.4, 0.598),
        ("D", (1.2, 2.2), 1.4, 0.80), ("D", (1.2, 2.6), 1.4, 0.74), ("D2", (1.2, 3.0), 1.4, 0.624),
        ("M", (1.2, 3.0), 1.4, 0.614)]
TW = semi.TW; DT = float(TW[1] - TW[0])


def local_spec(X, t0, t1, nfft=256):
    i = (TW >= t0) & (TW < t1); w = np.hanning(int(i.sum()))[:, None]
    F = np.fft.rfft(X[i] * w, nfft, axis=0); return np.mean(np.abs(F) ** 2, axis=1)


def weights(x, nz):
    f = np.fft.rfftfreq(256, DT); band = (f >= 5) & (f <= 120); centres, rows = [], []
    for t0 in np.arange(1.2, 2.81, 0.1):
        s = local_spec(x, t0, t0 + 0.2)[band]; n = local_spec(nz, t0, t0 + 0.2)[band]; ff = f[band]
        r = s / (s + n); centres.append(t0 + 0.1)
        rows.append(dict(P=np.sum(ff ** 2 * s), F=np.sum(ff ** 2 * r), F2=np.sum(ff ** 2 * r ** 2), B=np.sum(r)))
    c = np.array(centres)
    return {k: np.interp(TW, c, np.array([r[k] for r in rows])) for k in rows[0]}


def kappa(a, phi, win, t_anchor, cm):
    m = (TW >= win[0] - 1e-9) & (TW <= win[1] + 1e-9); a0 = float(lib.absorption(cm, [t_anchor])[0])
    return a0 * np.sum(phi[m] * a[m]) / np.sum(phi[m] * a[m] ** 2)


if __name__ == "__main__":
    base = gen.prepare(STRUCT, NOISE); out = {}
    for p, prof in PROF.items():
        cm = dict(base["cm"]);
        if prof is not None:
            cm["Q"] = np.asarray(prof, float)
        ctx = dict(base, cm=cm, w_in=lib.wavelets(base["D"], TW, cm)["w_in"])
        g = gen.signal(ctx, 1.0, SEED); nz = MULT * g["level"] * ctx["noise"][0]
        W = weights(g["x"], nz); W["U"] = np.ones_like(TW); a = lib.absorption(cm, TW); out[p] = dict(cm=cm, a=a, W=W)
    rows = []
    print(f"{'profile':6s} {'win':9s} {'anc':4s} {'meas':>5s} " + " ".join(f"{k:>6s}" for k in ("U", "P", "F", "F2", "B")))
    for p, win, ta, K in MEAS:
        o = out[p]; ks = {k: float(kappa(o["a"], o["W"][k], win, ta, o["cm"])) for k in ("U", "P", "F", "F2", "B")}
        rows.append(dict(profile=p, window=win, anchor=ta, measured=K, **ks))
        print(f"{p:6s} {win[0]}-{win[1]:<5} {ta:<4} {K:5.3f} " + " ".join(f"{ks[k]:6.3f}" for k in ("U", "P", "F", "F2", "B")))
    for k in ("U", "P", "F", "F2", "B"):
        e = np.array([r[k] - r["measured"] for r in rows]); print(k, "max|err|", round(float(np.abs(e).max()), 3), "rms", round(float(np.sqrt((e ** 2).mean())), 3))
    prof_a = {p: dict(a_top=float(o["a"][0]), a_mid=float(np.interp(2.1, TW, o["a"])), a_end=float(o["a"][-1])) for p, o in out.items()}
    json.dump(dict(seed=SEED, rows=rows, absorption=prof_a), open("../analysis/dev_phi2_model.json", "w"), indent=1)
