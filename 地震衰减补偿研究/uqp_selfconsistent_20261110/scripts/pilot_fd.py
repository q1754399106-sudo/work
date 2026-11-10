"""Exploratory pilot on the used FD data H1-H6 (no pre-registered rules, NOTE_ZH.md): how far are W0 / UQ-W from a
cell-wise Wiener ceiling, and does the self-consistent operator (uqp.py) close the gap?

2 x 2 design at every Q multiplier m (assumed Q = m x true Q; true attenuation multiplier u = m):
                      measured local SNR (uqw.py)     reference-level signal power (uqp.py)
  prior (sig 0.447)   UQW                             UQA
  posterior           PUQW                            PUQ      <- UQ-P
plug-in (sig 0) versions W0, A0 (assumed Q) and PW0, PA0 (fitted u). Ceiling ORA: true Q, noise power of the
actual noise, signal power of the elastic target, both per 0.4 s window.

usage: pilot_fd.py diag           noise-free attenuation diagnostics (no scores)
       pilot_fd.py run K N LEVEL  worker K of N; LEVEL 14dB (seeds of stage 2b / W) or 30dB (new seeds)
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402
from uqp import np, lib, semi, ROOT, DT, SIG0, UGRID  # noqa: E402
import fdcheck as fd  # noqa: E402
from uqw import uqw, snr_rows, G_IQ  # noqa: E402
from stage0 import M_ALL  # noqa: E402

OUT = ROOT / "results" / "pilot_fd"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ)
LEVELS = {"14dB": (14.0, 20261113000), "30dB": (30.0, 20261121000)}
KEYS = ("W0", "UQW", "A0", "UQA", "PW0", "PA0", "PUQW", "PUQ")


def source_log_power(fp, f):
    return np.log(np.abs(np.fft.rfft(fd.ricker(fp), 2 * (len(f) - 1))) ** 2 + 1e-300)


def diag():
    """Noise-free: attenuation multiplier of the FD data against the layered constant-Q model, with the Ricker
    spectrum as source (white reflectivity), with a free common spectrum, and with the elastic target's window
    spectrum as source (colour-free)."""
    for name in fd.MODELS:
        D0, T, cm, tw, fp = fd.load_fd(name); x = lib.absorption(cm, tw)
        f, cen, C, N, P = uqp.win_spectra(D0, tw); PT = uqp.win_spectra(T, tw)[4]
        wt = uqp.weights(f, P); tab = uqp.atten_table(tw, x, x, f, cen)
        a = uqp.surface_fit(P, wt, tab, g=source_log_power(fp, f)); b = uqp.surface_fit(P, wt, tab, g=np.log(PT))
        c = uqp.surface_fit(P, wt, tab)
        print(name, "Q", np.round(cm["Q"], 1), "t*(1.75) %.4f" % lib.absorption(cm, [1.75])[0],
              " u: white %.3f  free spectrum %.3f  colour-free %.3f" % (a["u"], c["u"], b["u"]), flush=True)


def run(name, level, r):
    D0, T, cm, tw, fp = fd.load_fd(name); Q = cm["Q"]; db, base = LEVELS[level]
    seed = base + 1000 * (fd.MODELS.index(name) + 1) + r
    E = fd.noise(D0, db, seed); D = D0 + E; x1 = lib.absorption(cm, tw)
    f, cen, C, N, P = uqp.win_spectra(D, tw); wt = uqp.weights(f, P)      # free common spectrum (diag: blue reflectivity)
    inv = snr_rows(D, tw, F_OP); sc = lambda X: float(fd.score_fd(X, T, tw)[0]); rec = []
    ap = lambda x, sig, rho: sc(uqp.apply_rows(uqp.operator(x, x, sig, rho, F_OP, cm["fref"], S), D))
    one = lambda tab, u: (lambda i: uqp.surface_fit(C, wt, tab[i:i + 1], ugrid=UGRID[i:i + 1]))(int(np.argmin(np.abs(UGRID - u))))
    # ceiling: true Q, power of the actual noise, signal power of the elastic target (window level)
    PT = uqp.win_spectra(T, tw)[4]; PE = uqp.win_spectra(E, tw)[4]
    rho_or = uqp.rho_rows(f, cen, PE, np.log(PT), tw, F_OP); ora = ap(x1, 0.0, rho_or)
    ora_u = {str(u): ap(x1 * u, 0.0, rho_or) for u in (0.7, 0.85, 1.2, 1.4)}
    for m in M_ALL:
        x = x1 / m; tab = uqp.atten_table(tw, x, x, f, cen)
        fit = uqp.surface_fit(C, wt, tab); up, sp, se = uqp.posterior(fit, C, wt, tab); uf = fit["u"]
        d = dict(m=m, u_fit=uf, se=se, u_post=up, sig_post=sp, n_rel=int((fit["rel"] >= 3).sum()))
        d["UQW"] = sc(uqw(D, cm, Q, m, SIG0, tw, 0.0, 0.0, inv_snr=inv)); d["W0"] = sc(uqw(D, cm, Q, m, 0.0, tw, 0.0, 0.0, inv_snr=inv))
        rho1 = uqp.rho_rows(f, cen, N, uqp.log_s0(one(tab, 1.0)), tw, F_OP)
        d["A0"] = ap(x, 0.0, rho1); d["UQA"] = ap(x, SIG0, rho1)
        d["PW0"] = sc(uqw(D, cm, Q, m / uf, 0.0, tw, 0.0, 0.0, inv_snr=inv))
        d["PA0"] = ap(x * uf, 0.0, uqp.rho_rows(f, cen, N, uqp.log_s0(fit), tw, F_OP))
        d["PUQW"] = sc(uqw(D, cm, Q, m / up, sp, tw, 0.0, 0.0, inv_snr=inv))
        d["PUQ"] = ap(x * up, sp, uqp.rho_rows(f, cen, N, uqp.log_s0(one(tab, up)), tw, F_OP))
        rec.append(d)
    return dict(model=name, level=level, r=r, seed=seed, ORA=ora, ORA_u=ora_u, records=rec)


if __name__ == "__main__":
    if sys.argv[1] == "diag":
        diag(); sys.exit(0)
    k, n, level = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(h, r) for h in fd.MODELS for r in fd.R]
    for j, (h, r) in enumerate(jobs):
        path = OUT / f"{h}_{level}_r{r}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); out = run(h, level, r); out["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(out, indent=1), encoding="utf-8"); q = out["records"][4]
        print(h, level, r, "ORA %.3f" % out["ORA"], " ".join("%s %.3f" % (key, q[key]) for key in KEYS),
              "u_fit(m=1) %.2f se %.2f" % (q["u_fit"], q["se"]), f"{out['seconds']:.0f}s", flush=True)
