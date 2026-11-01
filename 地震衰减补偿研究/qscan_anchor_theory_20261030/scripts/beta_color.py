"""Reflectivity-colour test (BETA_ZH.md step 6): is beta the colour of the wavelet window relative to the deeper data?

usage: beta_color.py    worker; jobs claimed with results16/<job>.claim, resumable.
Theory: a reflectivity colour tilt exp(pi f tau(t)) enters the mismatch exactly like the wavelet anchor:
delta(t) = (u - v) a(t) - a0 (u - v_a) + (tau0 - tau(t)), so beta_colour = sum a (tau0 - tau) / sum a^2
(tau0: colour in the wavelet window 1.2-1.6 s; uniform weights). tau(t) is measured on the TRUE reflectivity:
lateral-mean amplitude spectrum of 0.4 s windows (centres 1.4-2.8 s every 0.1 s), ln A = const + pi f tau over 8-60 Hz.
Data: experiment-14 recipe (structure 10227, noise line 10149 x 1), c = 1.0.
  L<seed>  block-bootstrap F-11 T2 log reflectivity (gen.signal), new seeds 10001-10012;
  W<seed>  white control: the same generator with the log replaced by a white Bernoulli-Gaussian series
           (semi.log_reflectivity patched; same length, density 0.3), seeds 10101-10106.
Scans: fixed anchor c_a 0.7 / 1.0 / 1.4, default wavelet window, lambda 1/2 (qprofile.scan).
-> results16/<set><seed>.json (beta from kappa_hat of the anchor pair; tau profile; beta_pred)
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from qprofile import scan

OUT = ROOT / "results16"; bm.MAXITER = 20000
JOBS = [f"L{s}" for s in range(10001, 10013)] + [f"W{s}" for s in range(10101, 10107)]
CENT = np.round(np.arange(1.4, 2.8 + 1e-9, 0.1), 3); DV = 1 / 0.7 - 1 / 1.4
_orig = semi.log_reflectivity


def white_log(seed):
    t, rc = _orig(); rng = np.random.default_rng(seed + 777)
    w = rng.normal(size=len(rc)) * (rng.random(len(rc)) < 0.3); w *= np.std(rc) / np.std(w)
    return lambda: (t, w)


def tau_profile(R):
    N = 1024; f = np.fft.rfftfreq(N, semi.FINE); b = (f >= 8) & (f <= 60); out = []
    for tc in CENT:
        i = (semi.TF >= tc - 0.2) & (semi.TF < tc + 0.2); X = R[i] * np.hanning(i.sum())[:, None]
        A = np.sqrt(np.mean(np.abs(np.fft.rfft(X, N, axis=0)) ** 2, axis=1))
        out.append(float(np.polyfit(np.pi * f[b], np.log(A[b]), 1)[0]))
    return np.array(out)


def kurt_shallow(R):
    """Kurtosis of the 6-60 Hz band-limited true reflectivity in 1.2-1.8 s (development correlate of beta)."""
    i = (semi.TF >= 1.2) & (semi.TF < 1.8); F = np.fft.rfft(R[i], axis=0); f = np.fft.rfftfreq(i.sum(), semi.FINE)
    F[(f < 6) | (f > 60)] = 0; x = np.fft.irfft(F, i.sum(), axis=0).ravel(); return float(np.mean(x ** 4) / np.mean(x ** 2) ** 2)


def beta_pred(tau, cm):
    a = lib.absorption(cm, CENT); return float(np.sum(a * (tau[0] - tau)) / np.sum(a ** 2))


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    for job in JOBS:
        path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic(); seed = int(job[1:])
        semi.log_reflectivity = white_log(seed) if job[0] == "W" else _orig
        g = gen.signal(ctx, 1.0, seed); semi.log_reflectivity = _orig
        D = g["x"] + g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm); v = {}
        for a in (0.7, 1.0, 1.4):
            wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * a)); r = scan(D, cm, Q, pp, lambda cc: wv)
            v[str(a)] = dict(v_sel=r["v_sel"], edge=r["edge"], ok=r["neighbourhood_converged"])
        k = (v["0.7"]["v_sel"] - v["1.4"]["v_sel"]) / DV; beta = v["1.0"]["v_sel"] - (1 - k) * 1.0 - k
        tau = tau_profile(g["R"])
        rec = dict(job=job, set=job[0], seed=seed, kappa=k, beta=beta, beta_pred=beta_pred(tau, cm), kurt_sh=kurt_shallow(g["R"]), tau=tau.tolist(),
                   centres=CENT.tolist(), scans=v, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"beta {beta:.3f} pred {rec['beta_pred']:.3f} kappa {k:.3f}", f"{rec['seconds']:.0f}s", flush=True)
