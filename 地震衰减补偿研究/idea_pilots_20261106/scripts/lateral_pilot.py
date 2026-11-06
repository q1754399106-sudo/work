"""Exploratory kill test on used data (Z group, structures A, B, C, all used in experiment 15).
Can same-stratum spectra measure the LATERAL change of t* without knowing colour, source or wavelet?
For horizon windows, per trace position: dip-corrected cross-power of +-HW neighbouring pairs on the flattened section,
ln amplitude relative to the section mean over the same stratum, slope over 10-60 Hz; dt*_est = -slope / pi.
Truth: true t* at the same stratum minus its section median. Reported per dataset and window."""
import json, sys
sys.path.insert(0, r"D:\ChatGPT\研究\地震衰减补偿研究\mu_field_20261106\scripts")
from env18 import L, np, semi, lib, git_text, RES

sq = L.sq; TW = semi.TW; DT = semi.DT
META = json.loads(git_text("semisynth_fault_q_20261030/analysis/datasets.json"))
DATA = RES / "semisynth_fault_q_20261030" / "cache" / "data"
HW = 8; WINS = [(1.8, 2.4), (2.1, 2.7), (2.4, 3.0)]; F = (10.0, 60.0)


def flatten(D, tb):
    n = D.shape[0]; X = np.fft.rfft(D, 2 * n, axis=0); f = np.fft.rfftfreq(2 * n, DT)
    return np.fft.irfft(X * np.exp(2j * np.pi * f[:, None] * tb[None, :]), 2 * n, axis=0)[:n]


def lateral(D, tau, a, b):
    w = semi.win(TW, a, b); tb = np.median(tau[w], axis=0); ok = (a + tb >= TW[0]) & (b + tb <= TW[-1])
    W = flatten(D, tb)[w] * lib.taper(w.sum())[:, None]
    X = np.fft.rfft(W, sq.NFFT, axis=0); f = np.fft.rfftfreq(sq.NFFT, DT)
    C = np.real(X[:, :-1] * np.conj(X[:, 1:]) * np.exp(-2j * np.pi * f[:, None] * sq.pair_lags(W)[None, :]))
    nx = D.shape[1]; band = (f >= F[0]) & (f <= F[1]); mean = lib.smooth(C[:, ok[:-1] & ok[1:]].mean(axis=1), 5); est = np.full(nx, np.nan)
    for j in range(HW, nx - HW):
        p = np.arange(j - HW, j + HW)
        if not (ok[p].all() and ok[p + 1].all()):
            continue
        c = lib.smooth(C[:, p].mean(axis=1), 5); g = band & (c > 0) & (mean > 0)
        if g.sum() > 10:
            est[j] = -np.polyfit(f[g], 0.5 * np.log(c[g] / mean[g]), 1)[0] / np.pi
    return est, tb, w


rows = []
for key in "ABC":
    for c0 in (0.75, 1.05, 1.4):
        for r in (0, 1):
            d = f"Z{key}_c{c0}_r{r}"; a = {k: v.astype(float) for k, v in np.load(DATA / f"{d}.npz").items()}
            tau = semi.structure(semi.local_lags(a["D"]))
            for wn in WINS:
                est, tb, w = lateral(a["D"], tau, *wn)
                tm = 0.5 * sum(wn); truth = np.array([np.interp(tm + tb[j], TW, a["tstar"][:, j]) for j in range(a["D"].shape[1])])
                g = np.isfinite(est); e = est[g] - np.median(est[g]); t = truth[g] - np.median(truth[g])
                gas = a["gas"][w].max(axis=0)[g] > 0.1
                rows.append(dict(id=d, win=wn, corr=float(np.corrcoef(e, t)[0, 1]), sd_true_ms=float(t.std() * 1e3), rmse_ms=float(np.sqrt(np.mean((e - t) ** 2)) * 1e3),
                                 slope=float(np.polyfit(t, e, 1)[0]), n=int(g.sum()), gas_true_ms=float(t[gas].mean() * 1e3) if gas.any() else None, gas_est_ms=float(e[gas].mean() * 1e3) if gas.any() else None))
            print(d, flush=True)
R = rows
for wn in WINS:
    s = [x for x in R if x["win"] == wn]
    print(f"window {wn}: median corr {np.median([x['corr'] for x in s]):.2f} (range {min(x['corr'] for x in s):.2f} to {max(x['corr'] for x in s):.2f}), true lateral sd {np.median([x['sd_true_ms'] for x in s]):.1f} ms, rmse {np.median([x['rmse_ms'] for x in s]):.1f} ms, regression slope {np.median([x['slope'] for x in s]):.2f}")
for key in "ABC":
    s = [x for x in R if x["id"][1] == key]; g = [x for x in s if x["gas_true_ms"] is not None]
    print(key, f"median corr {np.median([x['corr'] for x in s]):.2f}; under/around the lens (window touches it): true {np.median([x['gas_true_ms'] for x in g]):+.1f} ms, estimated {np.median([x['gas_est_ms'] for x in g]):+.1f} ms, n {len(g)}")
print("noise dependence (window 2.4-3.0 s):", {c: round(float(np.median([x['corr'] for x in R if x['win'] == WINS[2] and f'_c{c}_' in x['id']])), 2) for c in (0.75, 1.05, 1.4)})
