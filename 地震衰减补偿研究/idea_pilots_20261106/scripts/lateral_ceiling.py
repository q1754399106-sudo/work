"""Exploratory payoff ceiling on used data (Z group, structures A, B, C; c0 0.75, 1.05; r 0, 1). 70 dB inverse Q with:
  CONST  scalar c_eff_cum x model t*  (what experiment 15's oracle did)
  LAT    CONST + true LATERAL t* change only (t*(t, j) minus its median over traces, truth)
  TRUE   the full true t*(t, j)
  EST    CONST + lateral change estimated from same-stratum spectra (no truth), piecewise linear in time
Scores: S against the ideal target; local trace-window correlation below the gas lens and beside it (as experiment 15's masks)."""
import json, sys
sys.path.insert(0, r"D:\ChatGPT\研究\地震衰减补偿研究\mu_field_20261106\scripts")
sys.path.insert(0, r"C:\Users\qq175\AppData\Local\Temp\claude\cb810a85168e9b1c\0777d816-dbf9-4ae1-b342-2f76f5565fae\scratchpad")
from env18 import L, np, semi, lib, git_text, RES

lib12 = L.lib12; sq = L.sq; TW = semi.TW; DT = semi.DT
META = json.loads(git_text("semisynth_fault_q_20261030/analysis/datasets.json"))
DATA = RES / "semisynth_fault_q_20261030" / "cache" / "data"
cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); G = 70; NFFT = 1024; HW = 8; WINS = [(1.5, 2.1), (1.8, 2.4), (2.1, 2.7), (2.4, 3.0)]
f = np.fft.rfftfreq(NFFT, DT); LF = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); S2 = semi.gain_sigma2(G)
LAG = np.arange(len(TW))


def iq_tstar(D, ts):
    """lib12.iq / semi.inverse_q, with a t* profile per trace (ts: (nt, nx))."""
    out = np.empty_like(D); n = len(TW)
    for j in range(D.shape[1]):
        beta = np.exp(-np.pi * np.outer(ts[:, j], f)); lam = (beta + S2) / (beta ** 2 + S2)
        H = lam * np.exp(-2j * np.outer(ts[:, j], f) * LF); H[:, 0] = H[:, 0].real; H[:, -1] = H[:, -1].real
        M = np.fft.irfft(H, NFFT, axis=1); out[:, j] = M[np.arange(n)[:, None], (np.arange(n)[:, None] - LAG[None, :]) % NFFT] @ D[:, j]
    return out


def flatten(D, tb):
    n = D.shape[0]; X = np.fft.rfft(D, 2 * n, axis=0); ff = np.fft.rfftfreq(2 * n, DT)
    return np.fft.irfft(X * np.exp(2j * np.pi * ff[:, None] * tb[None, :]), 2 * n, axis=0)[:n]


def lateral_est(D):
    """Lateral t* change per trace at window centres, section median removed; interpolated in time, smoothed along x."""
    tau = semi.structure(semi.local_lags(D)); nx = D.shape[1]; cen = []; vals = []
    for a, b in WINS:
        w = semi.win(TW, a, b); tb = np.median(tau[w], axis=0); ok = (a + tb >= TW[0]) & (b + tb <= TW[-1])
        W = flatten(D, tb)[w] * lib.taper(w.sum())[:, None]; X = np.fft.rfft(W, sq.NFFT, axis=0); ff = np.fft.rfftfreq(sq.NFFT, DT)
        C = np.real(X[:, :-1] * np.conj(X[:, 1:]) * np.exp(-2j * np.pi * ff[:, None] * sq.pair_lags(W)[None, :]))
        band = (ff >= 10) & (ff <= 60); mean = lib.smooth(C[:, ok[:-1] & ok[1:]].mean(axis=1), 5); e = np.full(nx, np.nan)
        for j in range(HW, nx - HW):
            p = np.arange(j - HW, j + HW)
            if ok[p].all() and ok[p + 1].all():
                c = lib.smooth(C[:, p].mean(axis=1), 5); g = band & (c > 0) & (mean > 0)
                if g.sum() > 10:
                    e[j] = -np.polyfit(ff[g], 0.5 * np.log(c[g] / mean[g]), 1)[0] / np.pi
        g = np.isfinite(e); e = np.interp(np.arange(nx), np.where(g)[0], e[g]) - np.median(e[g]); cen.append(0.5 * (a + b)); vals.append(e)
    V = np.array(vals); first = np.interp(1.2, [1.2, cen[0]], [0, 1])
    return np.stack([np.interp(TW, [1.2] + cen, np.r_[0.0, V[:, j]]) for j in range(nx)], axis=1)


def masks(a):
    gas = a["gas"]; cols = gas.max(axis=0) > 0.5; base = int(np.max(np.where(gas[:, cols].max(axis=1) > 0.5)[0]))
    rows = np.zeros(len(TW), bool); rows[base + 12:min(base + 101, len(TW))] = True
    return rows[:, None] & cols[None, :], rows[:, None] & (gas.max(axis=0) < 0.01)[None, :]


def loc(o, ideal, m):
    return float(np.corrcoef(o[m], ideal[m])[0, 1])


rows = []
for key in "ABC":
    for c0 in (0.75, 1.05):
        for r in (0, 1):
            d = f"Z{key}_c{c0}_r{r}"; a = {k: v.astype(float) for k, v in np.load(DATA / f"{d}.npz").items()}; c = META[d]["c_eff_cum"]
            base = lib.absorption(dict(cm, Q=Q * c), TW)[:, None] * np.ones((1, a["D"].shape[1])); T = a["tstar"]; lat = T - np.median(T, axis=1, keepdims=True)
            sets = dict(CONST=base, LAT=base + lat, TRUE=T, EST=base + lateral_est(a["D"]))
            below, beside = masks(a); rec = dict(id=d)
            for k, ts in sets.items():
                o = iq_tstar(a["D"], ts); rec[k] = dict(S=lib12.score(o, a["ideal"])["S"], below=loc(o, a["ideal"], below), beside=loc(o, a["ideal"], beside))
            if key == "A" and c0 == 0.75 and r == 0:
                chk = lib12.iq(a["D"], cm, Q, c, c, G); o = iq_tstar(a["D"], base); print("check against lib12.iq, relative difference", float(np.linalg.norm(o - chk) / np.linalg.norm(chk)), flush=True)
            rows.append(rec); print(d, {k: round(rec[k]["S"], 3) for k in sets}, "below", {k: round(rec[k]["below"], 3) for k in sets}, flush=True)
print()
for scope, sel in (("all", lambda x: True), *[(f"structure {k}", (lambda k: lambda x: x["id"][1] == k)(k)) for k in "ABC"]):
    s = [x for x in rows if sel(x)]
    print(scope, " ".join(f"{m}: S {np.mean([x[m]['S'] for x in s]):.3f} below {np.mean([x[m]['below'] for x in s]):.3f} beside {np.mean([x[m]['beside'] for x in s]):.3f} |" for m in ("CONST", "LAT", "TRUE", "EST")))
for m in ("LAT", "TRUE", "EST"):
    dS = np.array([x[m]["S"] - x["CONST"]["S"] for x in rows]); dB = np.array([x[m]["below"] - x["CONST"]["below"] for x in rows])
    print(m, "minus CONST: S", f"{dS.mean():+.3f} ({(dS > 0).sum()}/{len(dS)} positive), below gas {dB.mean():+.3f} ({(dB > 0).sum()}/{len(dB)} positive)")
json.dump(rows, open(r"C:\Users\qq175\AppData\Local\Temp\claude\cb810a85168e9b1c\0777d816-dbf9-4ae1-b342-2f76f5565fae\scratchpad\lateral_ceiling.json", "w"), indent=1)
