"""Exploratory development of dip-steered UQ-MW (mwdip.py) on the used semi-synthetic development set (24 datasets,
batches 1-2: dipping Volve structure, faults, complex Q). No pre-registered rules (NOTE_ZH.md section 9).

Versions: cal (calibrated noise) and nd (4 x noise, as the UQ-W confirmation). At every Q multiplier m
(c = m x c_eff_cum): single-trace references UQW, W0, IQ70 (recomputed here on the same data) and
MW0 / UQMW with K = 1, 5, 9. SR route: IQ70, W0, MW0 at the SR estimate of c (sq.estimate on the same data).
Descriptive: K 9 without dip steering at m = 1, median coherence q with and without steering.

usage: pilot_semi.py run K N
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mwdip  # noqa: E402
from uqp import np, lib, semi, ROOT, DT, SIG0  # noqa: E402
from uqw import uqw, snr_rows, G_IQ, TW, lib12, load  # noqa: E402
from lib17 import DEV, score, sq  # noqa: E402
from stage0 import M_ALL  # noqa: E402

OUT = ROOT / "results" / "pilot_semi"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); ANCHOR = 1.4; MS = (0, 2, 4)


def run(did, version, cm, Q):
    arr, meta = load(did); D = arr["signal"] + (1.0 if version == "cal" else 4.0) * arr["noise"]; ideal = arr["ideal"]
    sc = lambda X: float(score(X, ideal)["S"]); inv = snr_rows(D, TW, F_OP); g = mwdip.Gather(D, TW, max(MS)); fb = (g.f >= 10) & (g.f <= 60)

    def xs(c):
        x = lib.absorption(dict(cm, Q=Q * c), TW); return x, x - lib.absorption(dict(cm, Q=Q * c), [ANCHOR])[0]

    rec = []
    for m in M_ALL:
        c = m * meta["c_eff_cum"]; x, dx = xs(c); mod = g.model(x, dx); d = dict(m=m)
        d["UQW"] = sc(uqw(D, cm, Q, c, SIG0, TW, ANCHOR, 0.0, inv_snr=inv)); d["W0"] = sc(uqw(D, cm, Q, c, 0.0, TW, ANCHOR, 0.0, inv_snr=inv))
        d["IQ70"] = sc(lib12.iq(D, cm, Q, c, c, G_IQ))
        for k in MS:
            d[f"MW0_K{2 * k + 1}"] = sc(g.apply(x, dx, 0.0, mod, F_OP, cm["fref"], S, k))
            d[f"UQMW_K{2 * k + 1}"] = sc(g.apply(x, dx, SIG0, mod, F_OP, cm["fref"], S, k))
        rec.append(d)
    c_sr = float(sq.estimate(D, cm, Q, nboot=0)["c"]); x, dx = xs(c_sr); mod = g.model(x, dx)
    sr = dict(c_sr=c_sr, ratio=c_sr / meta["c_eff_cum"], IQ70=sc(lib12.iq(D, cm, Q, c_sr, c_sr, G_IQ)),
              W0=sc(uqw(D, cm, Q, c_sr, 0.0, TW, ANCHOR, 0.0, inv_snr=inv)))
    for k in MS:
        sr[f"MW0_K{2 * k + 1}"] = sc(g.apply(x, dx, 0.0, mod, F_OP, cm["fref"], S, k))
    g0 = mwdip.Gather(D, TW, max(MS), steer=False); x, dx = xs(meta["c_eff_cum"]); m0 = g0.model(x, dx); m1 = g.model(x, dx)
    diag = dict(q_steer=float(np.median(m1[1][:, fb])), q_nosteer=float(np.median(m0[1][:, fb])),
                MW0_K9_nosteer=sc(g0.apply(x, dx, 0.0, m0, F_OP, cm["fref"], S, 4)), dip_abs_med_ms=float(1e3 * np.median(np.abs(g.p))),
                dip_abs_p95_ms=float(1e3 * np.quantile(np.abs(g.p), 0.95)))
    return dict(id=did, version=version, c_eff_cum=meta["c_eff_cum"], records=rec, sr=sr, diag=diag)


if __name__ == "__main__":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [(d, v) for v in ("cal", "nd") for d in DEV]
    for j, (did, v) in enumerate(jobs):
        path = OUT / f"{did}_{v}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); out = run(did, v, cm, Q); out["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(out, indent=1), encoding="utf-8"); q = out["records"][4]
        print(did, v, " ".join("%s %.3f" % (a, b) for a, b in q.items() if a != "m"), "| SR", {a: round(b, 3) for a, b in out["sr"].items()},
              "|", {a: round(b, 3) for a, b in out["diag"].items()}, f"{out['seconds']:.0f}s", flush=True)
