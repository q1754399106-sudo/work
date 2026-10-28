"""Part C of experiment 9 (PLAN_ZH.md): how much would the logged reflectivity colour bias SR at F-12 and F-14?

Logs only, no seismic. For every SR window (sq.windows(), 0.6 s every 0.3 s over 1.2-3.0 s) that lies inside the
well's Gardner-RC span (wq.well), the noise-free window spectrum is modelled as RC power x attenuation of a true
factor c: L_k(f) = 0.5 ln smooth(R_k, SM_C) + table_k(c). SR's fit (sq.mw_fit, double centring) over 10 Hz to
F_HI gives c_hat; the colour bias is c_hat(log RC) / c_hat(white RC). F_HI is the SR band limit of the inline
nearest to the well (experiment 10, lines.json, interpolated in inline number as tie.sr_factor).

usage: colour.py check   white-RC sanity only (exact recovery)
       colour.py run     writes results/colour.json
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq, tie  # noqa: E402
from wq import np, sq  # noqa: E402

WELLS = ["F-12", "F-14"]
C_TRUE = [0.77, 1.0, 1.3]
SM_C = 21


def sr_fhi(inline):
    L = json.loads((wq.IN / "specslope_q_volve_20261028" / "results" / "real" / "lines.json").read_text(encoding="utf-8"))
    k = sorted((int(a), v["f_hi"]) for a, v in L.items() if a != "10141")
    return float(np.interp(inline, [a for a, _ in k], [v for _, v in k]))


def bias(W, wins, f_hi, cm, Q, white=False):
    f = np.fft.rfftfreq(sq.NFFT, wq.DT); T = sq.mw_table(cm, Q, f, wins); sel = (f >= sq.F_LO) & (f <= f_hi)
    R = np.array([np.ones(len(f)) if white else wq.lib.smooth(wq.rc_power(W, a, b), SM_C) for a, b in wins])
    out = {}
    for c in C_TRUE:
        i = int(np.argmin(np.abs(sq.CGRID - c))); L = 0.5 * np.log(R) + T[i]
        out[str(c)] = sq.mw_fit(L, T, sel)[0] / float(sq.CGRID[i])
    return out


def main(mode):
    cm = wq.semi.curtain_cm(); Q = np.asarray(cm["Q"], float); P = wq.picks(); h, x, y, il, xl = tie.cube(); res = {}
    for n in WELLS:
        W = wq.well(n, P); wins = [(a, b) for a, b in sq.windows() if a >= W["t0"] and b <= W["t1"]]
        e, nn = np.mean(W["E"]), np.mean(W["N"]); line = int(il[int(np.argmin((x - e) ** 2 + (y - nn) ** 2))])
        fh = sr_fhi(line); r = dict(windows=wins, inline=line, f_hi=fh, white=bias(W, wins, fh, cm, Q, True))
        if mode == "run":
            r["log"] = bias(W, wins, fh, cm, Q)
        res[n] = r; print(n, {k: v for k, v in r.items() if k != "log"}, flush=True)
    if mode == "run":
        wq.ROOT.joinpath("results").mkdir(exist_ok=True)
        (wq.ROOT / "results" / "colour.json").write_text(json.dumps(res, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
