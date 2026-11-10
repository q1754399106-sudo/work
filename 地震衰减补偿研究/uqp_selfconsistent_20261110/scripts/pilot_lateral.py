"""Exploratory feasibility check (no pre-registered rules): how much does lateral averaging add on the used FD data at
14 dB, compared with the single-trace ceiling ORA? A boxcar mean over K adjacent traces before the single-trace
operators (no dip steering; the FD layering is nearly flat, one fault). usage: pilot_lateral.py"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402
from uqp import np, lib, semi, ROOT, DT, SIG0, UGRID  # noqa: E402
import fdcheck as fd  # noqa: E402
from uqw import uqw, snr_rows, G_IQ  # noqa: E402

F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); MS = (0.6, 1.0, 1.7); KS = (1, 3, 5, 9, 15); out = {}


def boxcar(D, k):
    if k == 1:
        return D
    pad = np.pad(D, ((0, 0), (k // 2, k // 2)), mode="edge"); c = np.cumsum(np.pad(pad, ((0, 0), (1, 0))), axis=1)
    return (c[:, k:] - c[:, :-k]) / k


for name in fd.MODELS:
    D0, T, cm, tw, fp = fd.load_fd(name); Q = cm["Q"]; x1 = lib.absorption(cm, tw); rec = {}
    D = D0 + fd.noise(D0, 14.0, 20261113000 + 1000 * (fd.MODELS.index(name) + 1) + 1)
    for k in KS:
        Dk = boxcar(D, k); inv = snr_rows(Dk, tw, F_OP); f, cen, C, N, P = uqp.win_spectra(Dk, tw); wt = uqp.weights(f, P)
        for m in MS:
            x = x1 / m; i1 = int(np.argmin(np.abs(UGRID - 1))); tab = uqp.atten_table(tw, x, x, f, cen, UGRID[i1:i1 + 1])
            rho = uqp.rho_rows(f, cen, N, uqp.log_s0(uqp.surface_fit(C, wt, tab, ugrid=UGRID[i1:i1 + 1])), tw, F_OP)
            rec[f"K{k}_m{m}"] = dict(UQW=float(fd.score_fd(uqw(Dk, cm, Q, m, SIG0, tw, 0.0, 0.0, inv_snr=inv), T, tw)[0]),
                                     A0=float(fd.score_fd(uqp.apply_rows(uqp.operator(x, x, 0.0, rho, F_OP, cm["fref"], S), Dk), T, tw)[0]))
        print(name, "K", k, " ".join("m%.1f UQW %.3f A0 %.3f |" % (m, rec[f"K{k}_m{m}"]["UQW"], rec[f"K{k}_m{m}"]["A0"]) for m in MS), flush=True)
    out[name] = rec
(ROOT / "results" / "pilot_lateral.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
