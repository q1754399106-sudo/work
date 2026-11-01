"""Second survey (PLAN_U1320A_ZH.md): two-anchor kappa on the U1320A lines (Gulf of Mexico, NAMSS 1988 stacks).

usage: u1320a_kappa.py [--smoke]    worker; jobs claimed with results18/<job>.claim, resumable.
Inputs: field_vsp_u1320a_20261024/inputs/<line>.mat (experiment 4: section 1.9-3.7 s at 4 ms, 99 traces, model cm with
water Q 1e6 above the seafloor (1.992 s) and the VSP Q_SR = 266 below). Analysis window 2.04-3.70 s.
Wavelet: the lib.wavelets recipe with the shallow window [t1, t2] and absorption at its centre (2.04-2.44 s / 2.24 s,
or 2.80-3.20 s / 3.00 s), the output band from the coherent SNR of 3.0-3.7 s, at the fixed anchor c_a.
Lasso objective of lib.qscan (lambda 1/2, sigma: haar per trace, K: bm.predictor) on a grid of v = 1/c from -1.0 to
5.0 in steps of 0.25 (v <= 0: no or negative attenuation; the attenuation is weak, t*(3.7 s) = 6.4 ms at v = 1);
scan from small v upwards, stop 3 points past the running minimum; vertex of the parabola through the minimum.
Anchors c_a = 2.0, 1.0, 0.5 (v_a 0.5, 1, 2).
-> results18/<line>__<window>_va<v_a>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, ROOT
from scipy.io import loadmat

OUT = ROOT / "results18"; bm.MAXITER = 20000
SRC = ROOT.parent / "field_vsp_u1320a_20261024" / "inputs"
LINES = ("3677-32", "4508-32"); WINS = {"w2.24": (2.04, 2.44, 2.24), "w3.00": (2.80, 3.20, 3.00)}; VAS = (0.5, 1.0, 2.0)
VG = np.round(np.arange(-1.0, 5.0 + 1e-9, 0.25), 4); T0, T1 = 2.04, 3.70


def load(line):
    d = loadmat(SRC / f"{line}.mat", simplify_cells=True); cm = dict(d["cm"]); t = np.asarray(cm["t"], float)
    w = (t >= T0 - 1e-9) & (t <= T1 + 1e-9); D = np.asarray(d["D"], float)[w]
    cmw = dict(cm, t=t[w], nt=int(w.sum()), Q=np.asarray(cm["Q"], float), vp=np.asarray(cm["vp"], float),
               h=np.asarray(cm["h"], float))
    return D, cmw


def wavelets(D, tw, cm, t1, t2, tc, nfft=512):
    f, snr_s, P_s = lib.coherence_snr(D, tw, t1, t2, nfft); _, snr_d, _ = lib.coherence_snr(D, tw, 3.00, 3.70, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * lib.absorption(cm, [tc])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))]); fs = lib.snr_band(f, snr_s, fpk); fd = lib.snr_band(f, snr_d, fpk)

    def cut(fc):
        g = np.ones_like(f); r = (f > fc - 5) & (f < fc)
        g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fc - 5)) / 5)); g[f >= fc] = 0; g[f < 3] = 0
        return g

    a_in = amp * cut(fs); a_out = a_in * cut(min(fd, fs))

    def zp(a):
        w = np.fft.irfft(a, nfft); w = np.roll(w, 100)[:201]; return w / np.max(np.abs(w))

    w_in, w_out = zp(a_in), zp(a_out)
    w_out *= np.max(np.abs(np.fft.rfft(w_in, nfft))) / np.max(np.abs(np.fft.rfft(w_out, nfft)))
    return dict(w_in=w_in, w_out=w_out, fmax_shallow=float(fs), fmax_target=float(fd))


def model_q(Q, v):
    with np.errstate(divide="ignore"):
        q = np.where(Q >= 1e5, Q, Q / v) if v != 0 else np.where(Q >= 1e5, Q, np.inf)
    return q


def scan(D, cm, Q, wv, pp):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    J = {}; conv = {}; prev = None
    for j in range(len(VG)):
        if J and j >= min(J, key=J.get) + 4:
            break
        cc = dict(cm); cc["Q"] = model_q(Q, VG[j]); B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, 0.5, 0.0, live, start=prev); prev = X
        res = D - B @ X
        J[j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + 0.5 * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
        conv[j] = bool(info["converged"])
    jm = min(J, key=J.get); edge = jm in (0, len(VG) - 1) or (jm + 1 not in J); v = float(VG[jm])
    if not edge and jm - 1 in J:
        ii = [jm - 1, jm, jm + 1]; a = np.polyfit(VG[ii], [J[x] for x in ii], 2)
        if a[0] > 0:
            v = float(np.clip(-a[1] / (2 * a[0]), VG[jm - 1], VG[jm + 1]))
    else:
        edge = True
    return dict(v_sel=v, v_grid=float(VG[jm]), edge=bool(edge), ok=bool(all(conv.get(x, True) for x in (jm - 1, jm, jm + 1))),
                J={str(VG[k]): J[k] for k in sorted(J)})


def prepare(D):
    live = np.abs(D).sum(axis=0) > 0; K, _, _ = bm.predictor(D, live)
    sigma = np.array([sc11.haar_sigma(D[:, j]) for j in range(D.shape[1])]); return dict(live=live, K=K, sigma=sigma)


def run(line, wname, va):
    D, cm = load(line); Q = np.asarray(cm["Q"], float); pp = prepare(D); t1, t2, tc = WINS[wname]
    wv = wavelets(D, cm["t"], dict(cm, Q=model_q(Q, va)), t1, t2, tc); t0 = time.monotonic(); r = scan(D, cm, Q, wv, pp)
    return dict(line=line, window=wname, v_anchor=va, band=[wv["fmax_shallow"], wv["fmax_target"]], seconds=time.monotonic() - t0, **r)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    if "--smoke" in sys.argv:
        r = run("3677-32", "w2.24", 1.0); r.pop("J"); print(json.dumps(r, indent=1)); sys.exit()
    for line in LINES:
        for wname in WINS:
            for va in VAS:
                job = f"{line}__{wname}_va{va}"; path = OUT / f"{job}.json"
                if path.exists():
                    continue
                try:
                    os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
                except FileExistsError:
                    continue
                rec = run(line, wname, va)
                path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
                print(job, f"v {rec['v_sel']:.3f}", "edge" if rec["edge"] else "", f"{rec['seconds']:.0f}s", flush=True)
