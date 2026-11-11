"""Why is the UQ-MW aperture too narrow at 4 x field noise? (exploratory, development data, NOTE_ZH.md section 14).

For each dataset and noise version, in the scoring window 2.3-2.8 s (dip-aligned neighbours, dip from the noisy data):
  rs[d], rn[d]   lateral correlation of the attenuated signal and of the noise at trace lag d = 0..8 (time domain)
  rd[d]          the same for the data
  q_fit          fitted geometric coherence gamma(2) / gamma(1) of mw.signal_model, median over 10-50 Hz
  n_ratio        estimated noise power / true noise power, 10-50 Hz (spectra of the window)
  w[j]           mean weight profile of MW0_K9 (m = 1) over 10-50 Hz in the window, normalised to sum 1

usage: diag_aperture.py ID [ID ...]
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mw, mwdip  # noqa: E402
from uqp import np, lib, semi, ROOT, DT  # noqa: E402
from uqw import G_IQ, TW, load  # noqa: E402

OUT = ROOT / "results" / "diag_aperture"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); M = 4; T0, T1, B0, B1 = 2.3, 2.8, 10.0, 50.0


def lat_corr(Z, p, r, lags=range(0, 2 * M + 1)):
    out = []
    for d in lags:
        A = mwdip.aligned(Z, TW, p, d)[r]; B = Z[r]; nx = Z.shape[1] - d
        out.append(float(np.sum(A[:, :nx] * B[:, :nx]) / np.sqrt(np.sum(A[:, :nx] ** 2) * np.sum(B[:, :nx] ** 2))))
    return out


def run(did, version, cm, Q):
    arr, meta = load(did); k = 1.0 if version == "cal" else 4.0; sig, noi = arr["signal"], k * arr["noise"]; D = sig + noi
    g = mwdip.Gather(D, TW, M); r = (TW >= T0) & (TW <= T1); c = meta["c_eff_cum"]
    x = lib.absorption(dict(cm, Q=Q * c), TW); dx = x - lib.absorption(dict(cm, Q=Q * c), [1.4])[0]; ls0, q, N = g.model(x, dx)
    wc = (g.cen >= T0) & (g.cen <= T1); fb = (g.f >= B0) & (g.f <= B1)
    _, _, Pn, _ = mw.lag_spectra(noi, TW); fo = (F_OP >= B0) & (F_OP <= B1)
    lf = np.log(np.maximum(F_OP, np.finfo(float).tiny) / float(cm["fref"])); A, Cinv, _ = uqp.expectations(x, dx, 0.0, F_OP, lf, "lognormal")
    rho = uqp.rho_rows(g.f, g.cen, N, ls0, TW, F_OP); qr = np.clip(mw.rows(g.f, g.cen, q, TW, F_OP), 0.05, 1.0); W = mw.weights_rows(qr, (S + rho) * Cinv, M)
    w = W[r][:, fo].mean(axis=(0, 1)); snr = 10 * np.log10(np.sum(sig[r] ** 2) / np.sum(noi[r] ** 2))
    return dict(id=did, version=version, snr_db=float(snr), rs=lat_corr(sig, g.p, r), rn=lat_corr(noi, g.p, r), rd=lat_corr(D, g.p, r),
                q_fit=float(np.median(q[wc][:, fb])), n_ratio=float(np.median(N[wc][:, fb] / Pn[wc][:, fb])), w=(w / w.sum()).tolist())


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); f3 = lambda v: " ".join("%.2f" % a for a in v)
    for did in sys.argv[1:]:
        for version in ("cal", "nd"):
            o = run(did, version, cm, Q); (OUT / f"{did}_{version}.json").write_text(json.dumps(o, indent=1), encoding="utf-8")
            print(did, version, "SNR %.1f dB  q_fit %.2f  N_est/N_true %.2f" % (o["snr_db"], o["q_fit"], o["n_ratio"]))
            print("  signal corr", f3(o["rs"])); print("  noise  corr", f3(o["rn"])); print("  data   corr", f3(o["rd"])); print("  weights    ", f3(o["w"]), flush=True)
