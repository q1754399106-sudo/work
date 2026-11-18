"""Residual-principle lambda of one-step Tikhonov at 30 and 24 dB (PLAN_LAMHI_ZH.md), Q correct, M12-M15, step 2b data
generator and selection loop unchanged (legacy2b/ = avo_bias_theory_20261106/scripts at HEAD, copied). 18 and 12 dB are
re-run as a reproduction check against results/lam_scan.json. D at the chosen lambda: lookup in results/lam_scan.json
(closed form, same half-decade grid). usage: lam_hi.py run | verdict"""
import sys, json, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE / "legacy2b"))
import step2b as s  # noqa: E402

SNRS = (30, 24, 18, 12)


def run():
    s.SNRS = SNRS; rows = []; t0 = time.monotonic()
    for M in ["M12", "M13", "M14", "M15"]:
        for fp in s.FPS:
            s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2
            for scene in s.SCENES:
                D, Cl, Tr, sig, meta, keep = s.columns(M, scene, fp); nc = D.shape[2]; a0, ax = s.a0_ax(s.NEW[M])
                Gd = [s.atten(a0 + ax[:, i]) for i in range(s.NX)]; GtG = [G.T @ G for G in Gd]; GtR = [Gd[i].T @ D[:, i, :] for i in range(s.NX)]
                target = np.sum(sig ** 2, axis=0) * s.n; ch = [None] * nc
                for lam in s.GRID:
                    O = np.stack([np.linalg.solve(GtG[i] + lam * s.I, GtR[i]) for i in range(s.NX)], axis=1)
                    res = sum(np.sum(((Gd[i] @ O[:, i, :] - D[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(s.NX))
                    for c in range(nc):
                        if ch[c] is None and res[c] <= target[c]:
                            ch[c] = float(lam)
                    if all(v is not None for v in ch):
                        break
                for c in range(nc):
                    rows.append(dict(M=M, fp=fp, scene=scene, snr=meta[c][0], seed=meta[c][1], lam=ch[c] if ch[c] is not None else float(s.GRID[-1]), fallback=ch[c] is None))
                print(M, fp, scene, "%.0fs" % (time.monotonic() - t0), flush=True)
    (ROOT / "results" / "lam_hi.json").write_text(json.dumps(rows))


def verdict():
    R = json.loads((ROOT / "results" / "lam_hi.json").read_text()); L = json.loads((ROOT / "results" / "lam_scan.json").read_text())
    grid = np.log10(L["grid"])

    def Dof(fp, lam):
        k = int(np.argmin(np.abs(grid - np.log10(lam))))
        if abs(grid[k] - np.log10(lam)) > 0.01:
            return np.nan  # outside the scanned grid
        return L["D"]["sig0.447_fp%g" % fp][k]
    out = {"fallbacks": int(sum(r["fallback"] for r in R)), "cells": {}}
    for snr in SNRS:
        for geom in ("aln", "nmo"):
            for fp in (60.0, 40.0, None):
                rr = [r for r in R if r["snr"] == snr and r["scene"].endswith("nmo") == (geom == "nmo") and (fp is None or r["fp"] == fp)]
                ll = np.log10([r["lam"] for r in rr]); dd = np.array([Dof(r["fp"], r["lam"]) for r in rr])
                out["cells"]["%d_%s_%s" % (snr, geom, "all" if fp is None else "%g" % fp)] = dict(
                    n=len(rr), log10lam_median=float(np.median(ll)), q25=float(np.percentile(ll, 25)), q75=float(np.percentile(ll, 75)),
                    nan_D=int(np.isnan(dd).sum()), absD_median=float(np.nanmedian(np.abs(dd))))
    a24 = out["cells"]["24_aln_all"]["absD_median"]; a30 = out["cells"]["30_aln_all"]["absD_median"]
    out["rule"] = "proceed_gate2_ge24dB" if a24 >= 0.15 else ("stop" if a30 < 0.092 else "report_user_ge30dB")
    (ROOT / "results" / "lam_hi_verdict.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"run": run, "verdict": verdict}[sys.argv[1]]()
