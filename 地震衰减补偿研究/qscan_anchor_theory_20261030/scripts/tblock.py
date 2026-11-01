"""Trace blocks on U1320A (PLAN_TBLOCK_ZH.md): full-grid scans that record the lasso objective per trace.

usage: tblock.py    worker; jobs claimed with results20/<job>.claim, resumable.
mu_rel = 0, so the traces are independent and J(v) = sum over traces; any block's selection follows from the
per-trace J. Data, wavelet (whole line), lambda 1/2, sigma as u1320a_kappa.py; grid v = -1..5 step 0.25, no early stop.
-> results20/<line>__<window>_va<v_a>.json  (v grid, per-trace J [grid x trace], live mask, converged)
"""
import os, json, time
from env15 import np, lib, bm, ROOT
import u1320a_kappa as uk

OUT = ROOT / "results20"; bm.MAXITER = 20000; LAM = 0.5


def run(line, wname, va):
    D, cm = uk.load(line); Q = np.asarray(cm["Q"], float); pp = uk.prepare(D); t1, t2, tc = uk.WINS[wname]
    wv = uk.wavelets(D, cm["t"], dict(cm, Q=uk.model_q(Q, va)), t1, t2, tc)
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    Jtr, conv = [], []; prev = None; t0 = time.monotonic()
    for v in uk.VG:
        B, _, scale = lib.dictionary(dict(cm, Q=uk.model_q(Q, v)), wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, LAM, 0.0, live, start=prev); prev = X
        res = D - B @ X
        jt = (0.5 * np.sum(res ** 2, axis=0) + LAM * np.sum(pen1 * np.abs(X), axis=0)) / s2
        Jtr.append(np.where(live, jt, 0.0).round(6).tolist()); conv.append(bool(info["converged"]))
    return dict(line=line, window=wname, v_anchor=va, v=uk.VG.tolist(), J_trace=Jtr, live=live.tolist(),
                sigma=sigma.tolist(), converged=conv, seconds=time.monotonic() - t0)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for line in uk.LINES:
        for wname in uk.WINS:
            for va in uk.VAS:
                job = f"{line}__{wname}_va{va}"; path = OUT / f"{job}.json"
                if path.exists():
                    continue
                try:
                    os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
                except FileExistsError:
                    continue
                rec = run(line, wname, va)
                path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
                print(job, f"{rec['seconds']:.0f}s", "all converged" if all(rec["converged"]) else "NOT all converged", flush=True)
