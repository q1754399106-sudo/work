"""Noise-floor explanation of the field spectral-ratio curvature (PLAN_ZH.md). Descriptive diagnostic, no compensation.

Per Volve inline (105 lines of field_baselines_volve_20261029/analysis/lines.json):
  original  c2 of L = 0.5 ln(P_deep/P_shallow), 10-55 Hz (identical to sid_uq_20261120 step0 Z2a)
  M1        dip-aligned adjacent-trace coherent part: S = |sum a conj(b_aligned)| / npair, N = Ppair - S
  M2        lateral second-difference residual: N = mean|x_i - (x~_{i-1}+x~_{i+1})/2|^2 / 1.5
  M3        tail window 4.0-4.4 s mean power (descriptive)
  debias    Ps = max(P - N, eps N), eps = 0.1; c2_deb; closed-form split of dc2 into deep and shallow terms.
The Volve cube is read (read-only) from the main checkout.
usage: run.py [nproc]
"""
import sys, json, subprocess, time
from pathlib import Path
from multiprocessing import Pool

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
MAIN = Path(r"D:/ChatGPT/研究/地震衰减补偿研究/field_baselines_volve_20261029")
sys.path.insert(0, str(MAIN / "scripts"))
import lib13 as L  # noqa: E402
from lib13 import np, lib, semi, DT  # noqa: E402

NF = 512; F = np.fft.rfftfreq(NF, DT); FIT = (F >= 10) & (F <= 55)
T_SH, T_DP, T_TL = (1.2, 1.6), (2.3, 2.8), (4.0, 4.4); EPS = 0.1; MAXLAG = 3; NS = 1126
TFULL = np.arange(NS) * DT


def lines():
    txt = subprocess.run(["git", "show", "HEAD:地震衰减补偿研究/field_baselines_volve_20261029/analysis/lines.json"],
                         cwd=REPO, capture_output=True, check=True).stdout.decode("utf-8")
    return json.loads(txt)["lines"]


def read_full(ln):
    """Inline ln, same traces as lib13.read_inline but the whole 0-4.5 s trace; dead traces are zero."""
    h, x, y, il, xl = L.cube(); D = np.zeros((NS, (L.XL1 - L.XL0) // L.XSTEP + 1))
    for j in np.where((il == ln) & (xl >= L.XL0) & (xl <= L.XL1) & ((xl - L.XL0) % L.XSTEP == 0))[0]:
        D[:, (xl[j] - L.XL0) // L.XSTEP] = semi.ibm(h[j, 240:].copy().view(">u4"))
    return D


def spec_all(D, a, b):
    """Complex spectra (NF/2+1, all traces, dead included) of the hann-tapered window a-b s; and the live mask."""
    k = (TFULL >= a - 1e-9) & (TFULL <= b + 1e-9)
    return np.fft.rfft(D[k] * np.hanning(k.sum())[:, None], NF, axis=0)


def lag(a, b):
    """Sub-sample lag tau (samples) maximising the circular cross-correlation of a and b within +-MAXLAG; b(t)~a(t-tau)."""
    c = np.fft.irfft(a * np.conj(b), NF)
    lags = np.arange(-MAXLAG, MAXLAG + 1); v = c[lags % NF]; i = int(np.argmax(v))
    if 0 < i < len(v) - 1:
        d = v[i - 1] - 2 * v[i] + v[i + 1]
        off = 0.5 * (v[i - 1] - v[i + 1]) / d if d < 0 else 0.0
    else:
        off = 0.0
    return -(lags[i] + off)


def align(b, tau):
    """Shift b by tau samples so that it lines up with the reference (b(t) ~ a(t - tau) -> b(t + tau))."""
    return b * np.exp(1j * 2 * np.pi * np.arange(len(b)) / NF * tau)


def sm(v):
    return lib.smooth(v, 5)


def noise_m1(X, live):
    idx = np.where(live)[0]; pairs = [(i, i + 1) for i in idx if i + 1 < X.shape[1] and live[i + 1]]
    cr = np.zeros(X.shape[0], complex); pp = np.zeros(X.shape[0]); taus = []
    for i, j in pairs:
        a, b = X[:, i], X[:, j]; t = lag(a, b); taus.append(t); bb = align(b, t)
        cr += a * np.conj(bb); pp += 0.5 * (np.abs(a) ** 2 + np.abs(b) ** 2)
    n = len(pairs); S = np.abs(sm(cr.real) + 1j * sm(cr.imag)) / n; Pp = sm(pp) / n
    return S, Pp - S, Pp, float(np.median(np.abs(taus)))


def noise_m2(X, live):
    acc = np.zeros(X.shape[0]); n = 0
    for i in range(1, X.shape[1] - 1):
        if not (live[i - 1] and live[i] and live[i + 1]):
            continue
        a = X[:, i]; l = align(X[:, i - 1], lag(a, X[:, i - 1])); r = align(X[:, i + 1], lag(a, X[:, i + 1]))
        acc += np.abs(a - 0.5 * (l + r)) ** 2; n += 1
    return sm(acc / n) / 1.5


def c2(v):
    return float(np.polyfit(F[FIT], v[FIT], 2)[0])


def one(ln):
    D = read_full(ln); live = np.abs(D[(TFULL >= 1.2) & (TFULL <= 3.0)]).sum(0) > 0
    out = dict(line=int(ln), n_live=int(live.sum())); W = {}
    for name, w in (("sh", T_SH), ("dp", T_DP), ("tl", T_TL)):
        X = spec_all(D, *w); P = sm(np.mean(np.abs(X[:, live]) ** 2, axis=1)); W[name] = dict(P=P)
        if name != "tl":
            S, N1, Pp, tau = noise_m1(X, live); N2 = noise_m2(X, live)
            W[name].update(S=S, N1=N1, N2=N2, Pp=Pp); out["tau_med_" + name] = tau
    r0 = 0.5 * np.log(W["dp"]["P"] / W["sh"]["P"]); out["c2_orig"] = c2(r0)
    for m in ("N1", "N2"):
        Ps = {}
        for w in ("sh", "dp"):
            N = np.maximum(W[w][m], 0.0); raw = W[w]["P"] - N; Ps[w] = np.maximum(raw, EPS * N)
            out["clip_%s_%s" % (m, w)] = float(np.mean((raw < EPS * N)[FIT]))
            q = N / Ps[w]; t = 0.5 * np.log1p(q)
            out["dc2_%s_%s" % (m, w)] = c2(t) if w == "dp" else -c2(t)
        rd = 0.5 * np.log(Ps["dp"] / Ps["sh"]); out["c2_deb_" + m] = c2(rd)
        out["dc2_pred_" + m] = out["dc2_%s_dp" % m] + out["dc2_%s_sh" % m]
    spectra = {w: {k: v.tolist() for k, v in W[w].items()} for w in W}
    return out, spectra


def summ(v):
    v = np.asarray(v); return dict(median=float(np.median(v)), iqr=[float(x) for x in np.percentile(v, [25, 75])], frac_pos=float(np.mean(v > 0)),
                                   n_pos=int(np.sum(v > 0)), n=int(len(v)))


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4; LN = lines(); t0 = time.time()
    with Pool(nproc) as pool:
        res = pool.map(one, LN, chunksize=4)
    rows = [r[0] for r in res]; specs = [r[1] for r in res]
    print("lines %d, %.0f s" % (len(rows), time.time() - t0))
    g = lambda k: np.array([r[k] for r in rows])  # noqa: E731
    S0 = summ(g("c2_orig")); print("original c2: median %+.3e IQR [%+.3e, %+.3e], c2>0 %d/%d" % (S0["median"], *S0["iqr"], S0["n_pos"], S0["n"]))
    out = dict(original=S0, lines=rows)
    for m in ("N1", "N2"):
        Sd = summ(g("c2_deb_" + m)); R = (S0["median"] - Sd["median"]) / S0["median"]
        dp, sh, pr = g("dc2_%s_dp" % m), g("dc2_%s_sh" % m), g("dc2_pred_" + m); obs = g("c2_orig")
        cc = np.corrcoef(pr, obs)[0, 1]; R2 = float(cc ** 2)
        print("%s: debiased c2 median %+.3e IQR [%+.3e, %+.3e], c2>0 %d/%d, drop R = %.3f" % (m, Sd["median"], *Sd["iqr"], Sd["n_pos"], Sd["n"], R))
        print("    dc2 deep term median %+.3e, shallow term median %+.3e, total %+.3e; median(pred)/m0 %.3f; across-line R2(obs~pred) %.3f (r %+.2f)" % (
            np.median(dp), np.median(sh), np.median(pr), np.median(pr) / S0["median"], R2, cc))
        print("    clipped fraction (10-55 Hz) median: shallow %.3f deep %.3f; lines with >20%% clipped: shallow %d deep %d" % (
            np.median(g("clip_%s_sh" % m)), np.median(g("clip_%s_dp" % m)), int(np.sum(g("clip_%s_sh" % m) > 0.2)), int(np.sum(g("clip_%s_dp" % m) > 0.2))))
        out[m] = dict(debiased=Sd, drop=float(R), dc2_deep_median=float(np.median(dp)), dc2_shallow_median=float(np.median(sh)),
                      dc2_pred_median=float(np.median(pr)), ratio_pred_obs=float(np.median(pr) / S0["median"]), R2=R2, r=float(cc))
    R1 = out["N1"]["drop"]; out["verdict"] = "holds" if R1 >= 0.5 else "fails"
    print("\nVERDICT (M1, R >= 0.50):", out["verdict"], "R = %.3f" % R1)
    # mean spectra over lines, descriptive table
    mean = {w: {k: np.mean([s[w][k] for s in specs], 0) for k in specs[0][w]} for w in specs[0]}
    fs = np.arange(8, 73, 4); at = lambda v, f: float(v[int(np.argmin(np.abs(F - f)))])  # noqa: E731
    tab = {}
    for w in ("sh", "dp"):
        tab[w] = dict(snr_M1={int(f): at(mean[w]["S"], f) / max(at(mean[w]["N1"], f), 1e-30) for f in fs},
                      snr_M2={int(f): at(mean[w]["P"] - mean[w]["N2"], f) / at(mean[w]["N2"], f) for f in fs},
                      tail_over_N1={int(f): at(mean["tl"]["P"], f) / max(at(mean[w]["N1"], f), 1e-30) for f in fs})
        print("%s  f: " % w + " ".join("%5d" % f for f in fs))
        for k in tab[w]:
            print("  %-13s" % k + " ".join("%5.2f" % tab[w][k][int(f)] for f in fs))
    out["mean_table"] = tab; out["tau_med"] = dict(sh=float(np.median(g("tau_med_sh"))), dp=float(np.median(g("tau_med_dp"))))
    print("median |tau| (samples): shallow %.2f deep %.2f" % (out["tau_med"]["sh"], out["tau_med"]["dp"]))
    out["mean_spectra"] = {w: {k: v.tolist() for k, v in mean[w].items()} for w in mean}; out["F"] = F.tolist()
    (ROOT / "results").mkdir(exist_ok=True); (ROOT / "results" / "run.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
