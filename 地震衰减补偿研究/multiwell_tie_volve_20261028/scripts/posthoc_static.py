"""Post-hoc control of experiment 9 (not pre-registered, descriptive): is the tie gain specific to the time-variant
Q compensation, or would a time-invariant spectral boost do the same?

STATIC_SR: the 30 dB inverse-Q operator of IQ30_SR frozen at the t* of 2.5 s (amplitude and dispersion phase),
applied to the whole trace as one filter. Same wells, traces, reference and tie as tie.py -> results/posthoc_static.json.
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq, tie  # noqa: E402
from wq import np, semi, lib, TW  # noqa: E402

T_FIX = 2.5


def static_iq(D, cm, Q, g_db, nfft=1024):
    c = dict(cm); c["Q"] = np.asarray(Q, float); ts = float(lib.absorption(c, [T_FIX])[0])
    f = np.fft.rfftfreq(nfft, wq.DT); s = semi.gain_sigma2(g_db); beta = np.exp(-np.pi * f * ts)
    H = (beta + s) / (beta ** 2 + s) * np.exp(-2j * f * ts * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])))
    H[0] = H[0].real; H[-1] = H[-1].real
    return np.fft.irfft(np.fft.rfft(D, nfft, axis=0) * H[:, None], nfft, axis=0)[:D.shape[0]]


def main():
    R = json.loads((wq.ROOT / "results" / "tie.json").read_text(encoding="utf-8"))
    P = wq.picks(); h, x, y, il, xl = tie.cube(); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    wt = semi.win(np.arange(1126) * wq.DT, semi.T0, semi.T1); rows = {}
    for n, r in R["wells"].items():
        g = tie.tie_log(n, P); w, idx = tie.path(g, x, y, il, xl); uniq = np.unique(idx); pos = {j: i for i, j in enumerate(uniq)}
        D = np.stack([semi.ibm(h[j, 240:].copy().view(">u4")) for j in uniq], axis=1)[wt]
        A = static_iq(D, cm, Q * r["c_sr"], tie.G_DB)
        tr = np.zeros(len(TW)); tr[w] = [A[ii, [pos[j] for j in q]].mean() for ii, q in zip(np.where(w)[0], idx)]
        ref = np.load(wq.ROOT / "cache" / f"tie_{n}.npz")["ref"]
        t2, sh = lib.tie(tr, ref, TW, semi.win(TW, g["t0"], g["t1"]))
        rows[n] = dict(STATIC_SR=t2, shift_ms=sh, IQ30_SR=r["scores"]["IQ30_SR"]["T2"], Input=r["scores"]["Input"]["T2"])
    G = R["groups"]
    d = [float(np.mean([rows[n]["IQ30_SR"] - rows[n]["STATIC_SR"] for n in m])) for m in G.values()]
    s = [float(np.mean([rows[n]["STATIC_SR"] - rows[n]["Input"] for n in m])) for m in G.values()]
    out = dict(t_fix=T_FIX, wells=rows, IQ30_SR_minus_STATIC_by_group=dict(zip(G, d)),
               STATIC_minus_Input_by_group=dict(zip(G, s)),
               mean_IQ30_SR_minus_STATIC=float(np.mean(d)), mean_STATIC_minus_Input=float(np.mean(s)))
    (wq.ROOT / "results" / "posthoc_static.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "wells"}, indent=1))


if __name__ == "__main__":
    main()
