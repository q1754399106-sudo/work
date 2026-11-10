"""Pre-freeze checks of PLAN_CONFIRM_MW_ZH.md on development data only (exploratory, NOTE_ZH.md section 13).

  P1  fd   H1-H6, 14 and 30 dB, 3 realizations, 9 Q multipliers: dip-steered mwdip.Gather (K 9) against mw.py without
           steering (K 9) and the single-trace UQ-W on the same noise; 30 dB was never run in development.
  P2  semi 24 development datasets, cal and 4 x noise, 9 Q multipliers: STK = UQ-W applied to the dip-aligned 9-trace
           mean (SNR re-estimated on the stacked section) - "is UQ-MW just stacking?". UQW / UQMW_K9 of pilot_semi.py
           are the references. P1 also records STK.

usage: pilot_prefreeze.py fd K N | semi K N
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mw, mwdip  # noqa: E402
from uqp import np, lib, semi, ROOT, DT, SIG0  # noqa: E402
from uqw import uqw, snr_rows, G_IQ, TW, load  # noqa: E402
from lib17 import DEV, score  # noqa: E402
from stage0 import M_ALL  # noqa: E402
import fdcheck as fd  # noqa: E402

OUT = ROOT / "results" / "pilot_prefreeze"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); M = 4
LEVELS = {"14dB": (14.0, 20261113000), "30dB": (30.0, 20261121000)}


def stack(g):
    return np.mean([g.nb[j] for j in range(-M, M + 1)], axis=0)


def run_fd(name, level, r):
    D0, T, cm, tw, fp = fd.load_fd(name); db, base = LEVELS[level]; seed = base + 1000 * (fd.MODELS.index(name) + 1) + r
    D = D0 + fd.noise(D0, db, seed); x1 = lib.absorption(cm, tw); sc = lambda X: float(fd.score_fd(X, T, tw)[0])
    g = mwdip.Gather(D, tw, M); f, cen, P, C = mw.lag_spectra(D, tw); inv = snr_rows(D, tw, F_OP); Ds = stack(g); invs = snr_rows(Ds, tw, F_OP); rec = []
    for m in M_ALL:
        x = x1 / m; mod = g.model(x, x); ls0, q, N = mw.signal_model(f, cen, P, C, tw, x, x)
        rec.append(dict(m=m, UQW=sc(uqw(D, cm, cm["Q"], m, SIG0, tw, 0.0, 0.0, inv_snr=inv)),
                        STK=sc(uqw(Ds, cm, cm["Q"], m, SIG0, tw, 0.0, 0.0, inv_snr=invs)),
                        UQMW_K9_steer=sc(g.apply(x, x, SIG0, mod, F_OP, cm["fref"], S, M)), MW0_K9_steer=sc(g.apply(x, x, 0.0, mod, F_OP, cm["fref"], S, M)),
                        UQMW_K9_flat=sc(mw.uqmw(D, tw, x, x, SIG0, ls0, q, N, f, cen, F_OP, cm["fref"], S, M)),
                        MW0_K9_flat=sc(mw.uqmw(D, tw, x, x, 0.0, ls0, q, N, f, cen, F_OP, cm["fref"], S, M))))
    return dict(model=name, level=level, r=r, seed=seed, dip_abs_p95_ms=float(1e3 * np.quantile(np.abs(g.p), 0.95)), records=rec)


def run_semi(did, version, cm, Q):
    arr, meta = load(did); D = arr["signal"] + (1.0 if version == "cal" else 4.0) * arr["noise"]
    g = mwdip.Gather(D, TW, M); Ds = stack(g); invs = snr_rows(Ds, TW, F_OP)
    rec = [dict(m=m, STK=float(score(uqw(Ds, cm, Q, m * meta["c_eff_cum"], SIG0, TW, 1.4, 0.0, inv_snr=invs), arr["ideal"])["S"])) for m in M_ALL]
    return dict(id=did, version=version, records=rec)


if __name__ == "__main__":
    mode, k, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    if mode == "fd":
        jobs = [(h, lv, r) for lv in LEVELS for h in fd.MODELS for r in fd.R]
    else:
        cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); jobs = [(d, v) for v in ("cal", "nd") for d in DEV]
    for j, job in enumerate(jobs):
        path = OUT / (mode + "_" + "_".join(str(a) for a in job) + ".json")
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); out = run_fd(*job) if mode == "fd" else run_semi(*job, cm, Q); out["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(out, indent=1), encoding="utf-8"); print(job, {a: round(b, 3) for a, b in out["records"][4].items()}, f"{out['seconds']:.0f}s", flush=True)
