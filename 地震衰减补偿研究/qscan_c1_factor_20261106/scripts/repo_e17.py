"""Repo pipeline (crit_kappa.py of qscan_anchor_theory_20261030, unchanged lasso solutions), adding eq. 17 of Cheng et al. 2025 (CJG).

usage: repo_e17.py dev | test     test: seeds 10201-10202 (same data as RESULTS_CRIT_ZH.md), anchors 0.7/1.4.
On every candidate: chengD (as crit_kappa.py, regression against results17), and eq. 17 of the output O = P X with the same rule as
scan_q.m of qscan_orig_solver_20261106: live-trace mean amplitude spectrum of the first and last T seconds of the analysis window,
band 5-90 Hz, (sum|D_s| - sum|D_d|)/sum|D_s|; from weak compensation the + to - crossing or the first minimum of |E1|.
"""
import sys, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; RR = HERE.parents[1]
sys.path.insert(0, str(RR / "qscan_anchor_theory_20261030" / "scripts"))
from env15 import np, lib, semi, sc11, bm, gen
OUT = HERE.parent / "results"; bm.MAXITER = 20000
QG = np.asarray(lib.QGRID, float); IDX = [j for j in range(len(QG)) if QG[j] >= 0.519]; VG = 1.0 / QG; T17 = 0.25; BAND = (5.0, 90.0)


def spec(O, live, t0, t1, nfft=256):
    i = semi.win(semi.TW, t0, t1); X = O[i][:, live] * np.hanning(i.sum())[:, None]
    return np.fft.rfftfreq(nfft, semi.DT), np.sqrt(np.mean(np.abs(np.fft.rfft(X, nfft, axis=0)) ** 2, axis=1))


def e17(O, live):
    n = int(round(T17 / semi.DT)); nf = int(round(1 / semi.DT)); X = O[:, live]
    f = np.fft.rfftfreq(nf, semi.DT); b = (f >= BAND[0]) & (f <= BAND[1])
    Ds = np.abs(np.fft.rfft(X[:n], nf, axis=0)).mean(axis=1); Dd = np.abs(np.fft.rfft(X[-n:], nf, axis=0)).mean(axis=1)
    return float((Ds[b].sum() - Dd[b].sum()) / Ds[b].sum())


def first_min(g, e):
    for k in range(len(e) - 1):
        if e[k] > 0 and e[k + 1] <= 0:
            return float(g[k] + (g[k + 1] - g[k]) * e[k] / (e[k] - e[k + 1])), True
        if k > 0 and abs(e[k]) < abs(e[k - 1]) and abs(e[k + 1]) >= abs(e[k]):
            a = np.polyfit(g[k - 1:k + 2], np.abs(e[k - 1:k + 2]), 2); v = float(g[k])
            if a[0] > 0:
                v = float(np.clip(-a[1] / (2 * a[0]), g[k - 1], g[k + 1]))
            return v, False
    return float("nan"), False


def centroid(w, dt):
    W = np.abs(np.fft.rfft(np.asarray(w, float).ravel(), 1024)); f = np.fft.rfftfreq(1024, dt); return float((f * W).sum() / W.sum())


def scan(D, cm, Q, pp, wv):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; rows = []; prev = None
    for j in sorted(IDX, reverse=True):
        cc = dict(cm); cc["Q"] = Q * QG[j]; B, P, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        X, info = bm.solve(B, np.eye(D.shape[0]), D, sigma, scale, K, 0.5, 0.0, live, start=prev); prev = X
        O = P @ X; f, As = spec(O, live, 1.2, 1.8); _, Ad = spec(O, live, 2.4, 3.0); b = (f >= 10) & (f <= 60)
        ch = float(abs(np.polyfit(f[b], np.log(np.maximum(Ad[b], 1e-30) / np.maximum(As[b], 1e-30)), 1)[0]))
        rows.append((float(VG[j]), ch, e17(O, live), bool(info["converged"])))
    g, ch, e, cv = map(np.array, zip(*rows))           # ascending v
    m = int(np.argmin(ch)); vc = float(g[m]); edge = m in (0, len(g) - 1)
    if not edge:
        a = np.polyfit(g[m - 1:m + 2], ch[m - 1:m + 2], 2)
        if a[0] > 0:
            vc = float(np.clip(-a[1] / (2 * a[0]), g[m - 1], g[m + 1]))
    v1, cross = first_min(g, e)
    return dict(v_chengD=vc, edge_chengD=bool(edge), v_e17=v1, cross_e17=cross, vgrid=g.tolist(), chengD=ch.tolist(), E1=e.tolist(), converged=cv.tolist())


if __name__ == "__main__":
    mode = sys.argv[1]; OUT.mkdir(exist_ok=True)
    ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    seeds = (9940,) if mode == "dev" else (10201, 10202); cs = (1.0,) if mode == "dev" else (0.7, 1.0, 1.4)
    for s in seeds:
        for c in cs:
            g = gen.signal(ctx, c, s); D = g["x"] + g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
            for a in (0.7, 1.4):
                t0 = time.monotonic(); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * a)); r = scan(D, cm, Q, pp, wv)
                r.update(dataset=f"S_c{c}_s{s}", u=1 / c, c_anchor=a, v_anchor=1 / a, f_in=centroid(wv["w_in"], semi.DT), seconds=time.monotonic() - t0)
                if mode == "test":
                    (OUT / f"repo_c{c}_s{s}__a{a}.json").write_text(json.dumps(r), encoding="utf-8")
                print(r["dataset"], "a", a, "chengD %.3f e17 %.3f cross %s f_in %.1f %.0fs" % (r["v_chengD"], r["v_e17"], r["cross_e17"], r["f_in"], r["seconds"]), "E1", np.round(r["E1"], 2).tolist(), flush=True)
