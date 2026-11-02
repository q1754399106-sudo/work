"""PLAN_BETASR_ZH.md: beta referenced to SR on the colour-test semi-synthetics (results16 L set), as on the field.
Data rebuilt as beta_color.py (gen.prepare(10227, [10149, 10273, 10355]), c = 1.0, D = x + level * noise[0]).
-> analysis/betasr.json"""
import json
from pathlib import Path
from env15 import np, gen, sq, semi, ROOT
import beta_color as bc

DV = 1 / 0.7 - 1 / 1.4
ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float); rows = {}
for job in bc.JOBS:
    r = json.loads((ROOT / "results16" / f"{job}.json").read_text()); seed = int(job[1:])
    semi.log_reflectivity = bc.white_log(seed) if job[0] == "W" else bc._orig
    g = gen.signal(ctx, 1.0, seed); semi.log_reflectivity = bc._orig
    D = g["x"] + g["level"] * ctx["noise"][0]; csr = float(sq.estimate(D, cm, Q)["c"])
    v = {k: x["v_sel"] for k, x in r["scans"].items()}; k = (v["0.7"] - v["1.4"]) / DV
    rows[job] = dict(set=job[0], kappa=k, c_sr=csr, beta_true=v["1.0"] - (1 - k) - k, beta_sr=v["1.0"] - (1 - k) / csr - k,
                     beta_true_stored=r["beta"])
    print(job, "c_sr %.3f beta_true %.3f beta_sr %.3f" % (csr, rows[job]["beta_true"], rows[job]["beta_sr"]), flush=True)
L = [x for x in rows.values() if x["set"] == "L"]; bs = np.array([x["beta_sr"] for x in L]); bt = np.array([x["beta_true"] for x in L])
sb = np.array([1 / x["c_sr"] - 1 for x in L]); FIELD = 0.278
q10, q90 = np.quantile(bs, [0.1, 0.9])
verdict = "consistent" if q10 <= FIELD <= q90 else ("field larger" if FIELD > bs.max() else "upper tail")
out = dict(rows=rows, L_beta_sr=dict(median=float(np.median(bs)), q10=float(q10), q90=float(q90), min=float(bs.min()), max=float(bs.max())),
           L_beta_true=dict(median=float(np.median(bt)), min=float(bt.min()), max=float(bt.max())),
           corr_sr_true=float(np.corrcoef(bs, bt)[0, 1]), corr_srbias_betatrue=float(np.corrcoef(sb, bt)[0, 1]),
           c_sr_median=float(np.median([x["c_sr"] for x in L])), field_beta_sr=FIELD, verdict=verdict)
(ROOT / "analysis" / "betasr.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print({k: v for k, v in out.items() if k != "rows"})
