"""Stage 2b (PLAN_STAGE2B_ZH.md): UQ-IQ on the independent 2-D viscoacoustic finite-difference data H1-H6
(fd2d_20261025: GSLS attenuation, fault diffraction, known Ricker source, time axis 0-1.8 s).

uqiq.py is imported unchanged (operator_rows). New code here only generalises the time axis and the anchor:
  iq_t   = lib12.iq line for line with TW replaced by a time axis argument (checked equal to lib12.iq on TW);
  uq_t   = uqiq.uqiq with the time axis and the anchor as arguments (sig = 0 checked equal to iq_t).
In the FD data the source wavelet is known (processing restores the Ricker spectrum at the source), so the anchor
is t = 0: the whole t* is uncertain (DESIGN_ZH.md 2.1: uncertainty starts where the wavelet is known).

usage: fdcheck.py smoke        equivalence checks only (no scores)
       fdcheck.py run K N      worker K of N over (model, level, realisation) -> results/stage2b/<H>_<lvl>_r<r>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqiq import ROOT, BASE, np, lib, semi, lib12, TW, G_IQ, operator_rows, load  # noqa: E402
from stage0 import M_ALL, SEL  # noqa: E402
import scipy.io as sio  # noqa: E402

FD = BASE / "fd2d_20261025"; MODELS = [f"H{i}" for i in range(1, 7)]
DT = float(lib12.DT); SIG = 0.447; SIG_DESC = (0.245, 0.894); WIN = (1.0, 1.75)
LEVELS = {"14dB": 14.0, "08dB": 8.0}; R = (1, 2, 3)
OUT = ROOT / "results" / "stage2b"


def load_fd(name):
    m = sio.loadmat(FD / "models" / f"{name}.mat", squeeze_me=True, struct_as_record=False)
    t = sio.loadmat(FD / "truth" / f"{name}_truth.mat", squeeze_me=True, struct_as_record=False)
    c = m["config"]; cm = dict(Q=np.asarray(c.Q, float), vp=np.asarray(c.vp, float), h=np.asarray(c.h, float), fref=float(c.fref))
    return np.asarray(m["clean"], float), np.asarray(t["target"], float), cm, np.asarray(c.t, float), float(c.fpeak)


def noise(D0, snr_db, seed):
    """noise_replica.m (AR(1) 0.45, 10% trace-common, Frobenius SNR) with numpy's generator instead of MATLAB's."""
    rng = np.random.default_rng(seed); u = rng.standard_normal(D0.shape)
    for i in range(1, len(u)):
        u[i] += 0.45 * u[i - 1]
    u = 0.9 * u + 0.1 * u.mean(axis=1, keepdims=True)
    return u * (np.linalg.norm(D0) / (10 ** (snr_db / 20) * np.linalg.norm(u)))


def ricker(fp, n=100):
    t = np.arange(-n, n + 1) * DT; a = (np.pi * fp * t) ** 2; return (1 - 2 * a) * np.exp(-a)


def iq_t(D, cm, Q, c_amp, c_ph, g_db, tw, nfft=1024):
    """lib12.iq with the time axis tw as an argument (identical statements)."""
    Q = np.asarray(Q, float); n = len(tw); f = np.fft.rfftfreq(nfft, DT); lag = np.arange(n)
    tp = lib.absorption(dict(cm, Q=Q * c_ph), tw)
    ta = None if c_amp is None else lib.absorption(dict(cm, Q=Q * c_amp), tw)
    s = None if c_amp is None else semi.gain_sigma2(g_db); M = np.empty((n, n))
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    for i in range(n):
        H = np.exp(-2j * f * tp[i] * lf)
        if ta is not None:
            beta = np.exp(-np.pi * f * ta[i]); H = H * (beta + s) / (beta ** 2 + s)
        H[0] = H[0].real; H[-1] = H[-1].real
        M[i] = np.fft.irfft(H, nfft)[(i - lag) % nfft]
    return M @ D


def uq_t(D, cm, Q, c, sig, tw, t_anchor, law="lognormal", g_db=G_IQ, nfft=1024):
    """uqiq.uqiq with the time axis and the anchor as arguments (operator_rows unchanged)."""
    Q = np.asarray(Q, float); n = len(tw); f = np.fft.rfftfreq(nfft, DT); lag = np.arange(n)
    x_hat = lib.absorption(dict(cm, Q=Q * c), tw); s = semi.gain_sigma2(g_db)
    dx_hat = x_hat - lib.absorption(dict(cm, Q=Q * c), [t_anchor])[0]
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    H = operator_rows(x_hat, dx_hat, sig, f, lf, s, law); M = np.empty((n, n))
    for i in range(n):
        h = H[i].copy(); h[0] = h[0].real; h[-1] = h[-1].real
        M[i] = np.fft.irfft(h, nfft)[(i - lag) % nfft]
    return M @ D


def score_fd(o, target, tw):
    """S = (M1 + M1_lag) / 2 as lib15.score, on the window WIN of the FD axis (shift +-20 ms, 1 ms step)."""
    w = (tw >= WIN[0] - 1e-9) & (tw <= WIN[1] + 1e-9)
    z = float(np.median(semi.colcorr(o[w], target[w])))
    lag = max(float(np.median(semi.colcorr(shift(o, s)[w], target[w])))
              for s in np.arange(-20, 21) / 1000)
    return 0.5 * (z + lag), z, lag


def shift(X, s):
    """run.shift on an axis starting at 0 (same statements; the axis origin cancels)."""
    t = np.arange(len(X)) * DT
    return np.stack([np.interp(t + s, t, a) for a in X.T], axis=1)


def records(name, level, r):
    D0, T, cm, tw, fp = load_fd(name); Q = cm["Q"]; i = MODELS.index(name) + 1
    seed = 20261113000 + 1000 * i + 10 * list(LEVELS).index(level) + r
    D = D0 + noise(D0, LEVELS[level], seed); w_out = ricker(fp); out = []
    add = lambda method, m, sig, X, **kw: out.append(dict(method=method, m=m, sig=sig, S=score_fd(X, T, tw)[0], **kw))
    for m in M_ALL:
        ref = iq_t(D, cm, Q, m, m, G_IQ, tw); add("IQ70", m, None, ref)
        x0 = uq_t(D, cm, Q, m, 0.0, tw, 0.0)
        add("UQ0", m, 0.0, x0, rel_diff=float(np.abs(x0 - ref).max() / np.abs(ref).max()))
        add("UQ", m, SIG, uq_t(D, cm, Q, m, SIG, tw, 0.0))
        for q in SIG_DESC:
            add("UQ", m, q, uq_t(D, cm, Q, m, q, tw, 0.0))
        add("UQ_anchor09", m, SIG, uq_t(D, cm, Q, m, SIG, tw, 0.9))        # wrong anchor (descriptive)
        add("VDB", m, None, lib12.tvw(iq_t(D, cm, Q, None, m, 0, tw), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"]))
        add("C1", m, None, iq_t(D, cm, Q, m, m, 30, tw))
        c2 = m * np.exp(SIG); add("C2", m, SIG, iq_t(D, cm, Q, c2, c2, G_IQ, tw))
    return out, seed


if __name__ == "__main__" and sys.argv[1] == "smoke":
    t0 = time.monotonic(); arr, meta = load("ZA_c0.75_r0"); cmz = semi.curtain_cm(); Qz = np.asarray(cmz["Q"], float)
    a = lib12.iq(arr["D"], cmz, Qz, 0.8, 0.8, G_IQ); b = iq_t(arr["D"], cmz, Qz, 0.8, 0.8, G_IQ, TW)
    a2 = lib12.iq(arr["D"], cmz, Qz, None, 0.8, 0); b2 = iq_t(arr["D"], cmz, Qz, None, 0.8, 0, TW)
    print("iq_t vs lib12.iq on TW: max|d| %.1e (amp) %.1e (phase only)" % (np.abs(a - b).max(), np.abs(a2 - b2).max()))
    import uqiq
    u1 = uqiq.uqiq(arr["D"], cmz, Qz, 0.8, 0.447, law="lognormal"); u2 = uq_t(arr["D"], cmz, Qz, 0.8, 0.447, TW, 1.4)
    print("uq_t vs uqiq on TW, anchor 1.4, sig 0.447: max|d| %.1e" % np.abs(u1 - u2).max())
    sh = shift(arr["D"], 0.012); sh0 = lib12.run.shift(arr["D"], 0.012)
    print("shift vs run.shift: max|d| %.1e" % np.abs(sh - sh0).max())
    for name in MODELS:
        D0, T, cm, tw, fp = load_fd(name)
        x = uq_t(D0, cm, cm["Q"], 1.0, 0.0, tw, 0.0); y = iq_t(D0, cm, cm["Q"], 1.0, 1.0, G_IQ, tw)
        print(name, "Q", np.round(cm["Q"], 1), "t* at 1.75 s %.4f" % lib.absorption(cm, [1.75])[0],
              "UQ0 vs IQ70 rel %.1e" % (np.abs(x - y).max() / np.abs(y).max()))
    print("%.0fs" % (time.monotonic() - t0))

if __name__ == "__main__" and sys.argv[1] == "run":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(h, lv, r) for h in MODELS for lv in LEVELS for r in R]
    for j, (h, lv, r) in enumerate(jobs):
        path = OUT / f"{h}_{lv}_r{r}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); rec, seed = records(h, lv, r)
        path.write_text(json.dumps(dict(model=h, level=lv, r=r, seed=seed, records=rec, seconds=time.monotonic() - t0), indent=1),
                        encoding="utf-8")
        print(h, lv, r, len(rec), f"{time.monotonic() - t0:.0f}s", flush=True)
