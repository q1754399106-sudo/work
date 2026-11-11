"""Semi-synthetic development set with spatially white noise (exploratory, NOTE_ZH.md section 14).

The "real" noise of the semi-synthetic sets is semi.residual_noise: d_j - (d_{j-1} + d_{j+1}) / 2, a lateral second
difference. It has no laterally flat component (lag correlation about +0.12, -0.37, -0.22), so a 9-trace mean removes
it about 3 times better (in power) than it removes white noise, and it breaks the "noise uncorrelated between traces"
assumption of mw.py. Here the noise is semi.gaussian_noise: spatially white, same time-variant spectrum and rms as
4 x the residual noise. 24 development datasets, 9 Q multipliers: UQW, STK (aligned 9-trace mean + UQ-W), UQMW / MW0 K 9.

usage: pilot_white.py run K N
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mwdip  # noqa: E402
from uqp import np, lib, semi, ROOT, DT, SIG0  # noqa: E402
from uqw import uqw, snr_rows, G_IQ, TW, load  # noqa: E402
from lib17 import DEV, score  # noqa: E402
from stage0 import M_ALL  # noqa: E402

OUT = ROOT / "results" / "pilot_white"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); M = 4; ANCHOR = 1.4; SEED = 20261122000


def run(did, j, cm, Q):
    arr, meta = load(did); noi = semi.gaussian_noise(4.0 * arr["noise"], np.random.default_rng(SEED + j)); D = arr["signal"] + noi
    sc = lambda X: float(score(X, arr["ideal"])["S"]); g = mwdip.Gather(D, TW, M); inv = snr_rows(D, TW, F_OP)
    Ds = np.mean([g.nb[i] for i in range(-M, M + 1)], axis=0); invs = snr_rows(Ds, TW, F_OP); rec = []
    for m in M_ALL:
        c = m * meta["c_eff_cum"]; x = lib.absorption(dict(cm, Q=Q * c), TW); dx = x - lib.absorption(dict(cm, Q=Q * c), [ANCHOR])[0]; mod = g.model(x, dx)
        rec.append(dict(m=m, UQW=sc(uqw(D, cm, Q, c, SIG0, TW, ANCHOR, 0.0, inv_snr=inv)), STK=sc(uqw(Ds, cm, Q, c, SIG0, TW, ANCHOR, 0.0, inv_snr=invs)),
                        UQMW_K9=sc(g.apply(x, dx, SIG0, mod, F_OP, cm["fref"], S, M)), MW0_K9=sc(g.apply(x, dx, 0.0, mod, F_OP, cm["fref"], S, M))))
    r = (TW >= 2.3) & (TW <= 2.8)
    return dict(id=did, seed=SEED + j, snr_db=float(10 * np.log10(np.sum(arr["signal"][r] ** 2) / np.sum(noi[r] ** 2))), records=rec)


if __name__ == "__main__":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for j, did in enumerate(DEV):
        path = OUT / f"{did}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); out = run(did, j, cm, Q); out["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(out, indent=1), encoding="utf-8")
        mean = {a: round(float(np.mean([q[a] for q in out["records"]])), 3) for a in ("UQW", "STK", "UQMW_K9", "MW0_K9")}
        print(did, "mean over m", mean, "| m=1", {a: round(b, 3) for a, b in out["records"][4].items()}, f"{out['seconds']:.0f}s", flush=True)
