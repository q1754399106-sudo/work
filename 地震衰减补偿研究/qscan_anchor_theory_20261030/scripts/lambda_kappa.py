"""lambda manipulation on U1320A (PLAN_LAMBDA_ZH.md): full-grid scans at lambda 1/16, 1/4, 1, 2.

usage: lambda_kappa.py    worker; jobs claimed with results21/<job>.claim, resumable.
Data, wavelet (whole line), sigma as u1320a_kappa.py; lambda enters bm.solve and the L1 term of J; grid v = -1..7
step 0.25, no early stop; per-trace J recorded (as tblock.py). lambda 1/2 is results20/.
-> results21/<line>__<window>_lam<lam>_va<v_a>.json
"""
import os, json, time
from env15 import np, lib, bm, ROOT
import u1320a_kappa as uk

OUT = ROOT / "results21"; bm.MAXITER = 20000
LAMS = (0.0625, 0.25, 1.0, 2.0); VG = np.round(np.arange(-1.0, 7.0 + 1e-9, 0.25), 4)


def run(line, wname, lam, va):
    D, cm = uk.load(line); Q = np.asarray(cm["Q"], float); pp = uk.prepare(D); t1, t2, tc = uk.WINS[wname]
    wv = uk.wavelets(D, cm["t"], dict(cm, Q=uk.model_q(Q, va)), t1, t2, tc)
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    Jtr, conv = [], []; prev = None; t0 = time.monotonic()
    for v in VG:
        B, _, scale = lib.dictionary(dict(cm, Q=uk.model_q(Q, v)), wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live, start=prev); prev = X
        res = D - B @ X
        jt = (0.5 * np.sum(res ** 2, axis=0) + lam * np.sum(pen1 * np.abs(X), axis=0)) / s2
        Jtr.append(np.where(live, jt, 0.0).round(6).tolist()); conv.append(bool(info["converged"]))
    return dict(line=line, window=wname, lam=lam, v_anchor=va, v=VG.tolist(), J_trace=Jtr, live=live.tolist(),
                converged=conv, seconds=time.monotonic() - t0)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for lam in LAMS:
        for line in uk.LINES:
            for wname in ("w3.00", "w2.24"):
                for va in uk.VAS:
                    job = f"{line}__{wname}_lam{lam}_va{va}"; path = OUT / f"{job}.json"
                    if path.exists():
                        continue
                    try:
                        os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
                    except FileExistsError:
                        continue
                    rec = run(line, wname, lam, va)
                    path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
                    print(job, f"{rec['seconds']:.0f}s", "all converged" if all(rec["converged"]) else "NOT all converged", flush=True)
