"""Descriptive: closed-form D(lambda) = dGA(UQ sig) - dGA(Tikhonov), u = 1, M12-M15, on a half-decade lambda grid
(gate 1 used only 1e-4, 1e-3, 1e-2; closed form validated there, P1 slope 0.974). Also the step-2b residual-principle lambda."""
import sys, json, gzip, subprocess
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate1 as g

ROOT = Path(__file__).resolve().parents[1]
LAMS = 10.0 ** np.arange(-4.5, -0.99, 0.5)
out = {"grid": LAMS.tolist(), "D": {}, "tik": {}, "uq": {}}
for sig in (0.447, 0.25):
    for fp in g.FPS:
        S2 = g.spec(fp) ** 2; rows = []
        for M in g.GATE_MODELS:
            a_ev, s2, co = g.paths(g.MODELS[M])
            EG0, EG20 = g.moments(a_ev, 0.0); EG1, EG21 = g.moments(a_ev, sig)
            for lam in LAMS:
                H0 = np.conj(EG0) / (EG20 + lam); H1 = np.conj(EG1) / (EG21 + lam)
                for j in range(3):
                    d0 = g.dga(g.gamma_closed(a_ev[j], 1.0, H0[j], S2), s2[j], co[j]); d1 = g.dga(g.gamma_closed(a_ev[j], 1.0, H1[j], S2), s2[j], co[j])
                    rows.append((lam, d0, d1))
        r = np.array(rows); k = "sig%.3f_fp%g" % (sig, fp)
        out["D"][k] = [float(np.median(r[r[:, 0] == l, 2] - r[r[:, 0] == l, 1])) for l in LAMS]
        out["tik"][k] = [float(np.median(np.abs(r[r[:, 0] == l, 1]))) for l in LAMS]
        out["uq"][k] = [float(np.median(np.abs(r[r[:, 0] == l, 2]))) for l in LAMS]
# residual-principle lambda from step 2b (Q correct), read from git
lamsel = {}
for M in g.GATE_MODELS:
    b = subprocess.run(["git", "show", "HEAD:地震衰减补偿研究/avo_bias_theory_20261106/results/step2b_%s.json.gz" % M], capture_output=True, cwd=ROOT).stdout
    for r in json.loads(gzip.decompress(b)):
        if r["qf"] == 1.0:
            lamsel.setdefault("%s_fp%g_snr%d_%s" % (r["m"], r["fp"], r["snr"], "nmo" if r["scene"].endswith("nmo") else "aln"), []).append(np.log10(r["lam"]))
out["step2b_log10lam"] = {k: dict(median=float(np.median(v)), q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)), n=len(v)) for k, v in sorted(lamsel.items())}
(ROOT / "results" / "lam_scan.json").write_text(json.dumps(out, indent=1))
print("log10 lam:", " ".join("%6.1f" % x for x in np.log10(LAMS)))
for name in ("D", "tik", "uq"):
    for k, v in out[name].items():
        print("%-4s %-18s" % (name, k), " ".join("%6.3f" % x for x in v))
for k, v in out["step2b_log10lam"].items():
    if k.startswith("direct"):
        print(k, v)
