"""Lasso worker: claims ../results/lasso_<job>.claim, resumable. usage: run_lasso.py"""
import os, json, time, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from lasso_scan import np, scan
from scipy.io import loadmat
OUT = HERE.parent / "results"; C = HERE.parent / "cache"; F0 = 30.0; TS = 1.0
FW = {"known": lambda fp: F0, "W3_0": lambda fp: fp, "W3_1": lambda fp: 1 / np.sqrt(1 / fp**2 - np.pi * TS * 0.01 / (2 * fp)),
      "F0.8": lambda fp: 0.8 * F0, "F1.25": lambda fp: 1.25 * F0}
for s in range(1201, 1207):
    for Qt in (50, 70, 100, 140):
        for cond, ffw in FW.items():
            name = f"lasso_Q{Qt}_s{s}__{cond}"; path = OUT / f"{name}.json"
            if path.exists():
                continue
            try:
                os.close(os.open(OUT / f"{name}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
            except FileExistsError:
                continue
            m = loadmat(C / f"d_Q{Qt}_s{s}.mat"); fp = float(m["fp"].ravel()[0]); fw = float(ffw(fp)); t0 = time.monotonic()
            S = scan(m["dd"], fw); S.update(solver="lasso", Qt=Qt, seed=s, cond=cond, fp=fp, fw=fw, seconds=time.monotonic() - t0)
            path.write_text(json.dumps(S), encoding="utf-8"); (OUT / f"{name}.claim").unlink()
            print(name, "fw %.1f vJ %.4f v1 %.4f (%s) vd %.4f v2 %.4f conv %.2f %.0fs" % (fw, S["vJ"], S["v1"], S["cross"], S["vd"], S["v2"], np.mean(S["converged"]), S["seconds"]), flush=True)
