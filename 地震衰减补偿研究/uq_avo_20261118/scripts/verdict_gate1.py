"""Verdict for gate 1 (PLAN_GATE1_ZH.md section 4) plus descriptive tables."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
R = json.loads((ROOT / "results" / "gate1.json").read_text())
GM = ["M12", "M13", "M14", "M15"]; SIG = 0.447


def pick(**k):
    return [r for r in R if all(r[a] == b for a, b in k.items())]


def key(r):
    return (r["model"], r["fp"], r["lam"], r["ev"])


def D_table(kind_hi="uq", sig_hi=SIG, models=GM, field="dga_num", u=1.0):
    lo = {key(r): r[field] for r in pick(kind="uq", sig=0.0, u=u) if r["model"] in models}
    hi = {key(r): r[field] for r in pick(kind=kind_hi, sig=sig_hi, u=u) if r["model"] in models}
    return {k: hi[k] - lo[k] for k in hi}


out = {}
D = D_table(); d = np.array(list(D.values()))
out["n_cells"] = len(d)
out["G1_median_absD"] = float(np.median(np.abs(d))); out["G1_pass"] = bool(out["G1_median_absD"] >= 0.092)
out["G2_frac_neg"] = float(np.mean(d < 0)); out["G2_pass"] = bool(out["G2_frac_neg"] >= 0.8)
pts = [r for r in pick(kind="uq", sig=SIG) if r["model"] in GM and "dga_num" in r]
x = np.array([r["dga_num"] for r in pts]); y = np.array([r["dga_closed"] for r in pts])
sl = float(np.polyfit(x, y, 1)[0]); r2 = float(np.corrcoef(x, y)[0, 1] ** 2)
out["P1_n"] = len(pts); out["P1_slope"] = sl; out["P1_R2"] = r2; out["P1_pass"] = bool(0.9 <= sl <= 1.1 and r2 >= 0.95)
pts0 = [r for r in pick(kind="uq", sig=0.0, u=1.0) if r["model"] in GM]
x0 = np.array([r["dga_num"] for r in pts0]); y0 = np.array([r["dga_closed"] for r in pts0])
out["P1_tik_slope"] = float(np.polyfit(x0, y0, 1)[0]); out["P1_tik_R2"] = float(np.corrcoef(x0, y0)[0, 1] ** 2)

# descriptive
desc = {}
def med(vals):
    return float(np.median(vals)) if len(vals) else None
for grp, f in (("fp", lambda k: k[1]), ("lam", lambda k: k[2]), ("ev", lambda k: k[3]), ("model", lambda k: k[0])):
    g = {}
    for k, v in D.items():
        g.setdefault(str(f(k)), []).append(v)
    desc["D_by_" + grp] = {a: dict(median=med(b), median_abs=med(np.abs(b))) for a, b in g.items()}
desc["absdGA_tik_u1"] = med([abs(r["dga_num"]) for r in pick(kind="uq", sig=0.0, u=1.0) if r["model"] in GM])
desc["absdGA_uq_u1"] = med([abs(r["dga_num"]) for r in pick(kind="uq", sig=SIG, u=1.0) if r["model"] in GM])
D25 = np.array(list(D_table(sig_hi=0.25, field="dga_num").values()))
desc["sig0.25_median_absD"] = med(np.abs(D25)); desc["sig0.25_frac_neg"] = float(np.mean(D25 < 0))
DM1 = np.array(list(D_table(models=["M1"]).values()))
desc["M1_median_absD"] = med(np.abs(DM1)); desc["M1_frac_neg"] = float(np.mean(DM1 < 0))
# Q error (closed form, all u)
ue = {}
for u in [0.5, 0.7, 0.8, 1.0, 1.2, 1.4, 2.0]:
    row = {}
    for name, kind, sig in (("tik", "uq", 0.0), ("uq", "uq", SIG), ("oc0", "oc", 0.0), ("oc", "oc", SIG)):
        v = [r["dga_closed"] for r in pick(kind=kind, sig=sig, u=u) if r["model"] in GM]
        row[name] = dict(median_abs=med(np.abs(v)), median=med(v))
    ue[str(u)] = row
desc["by_u_closed"] = ue
for name, kind, sig in (("tik", "uq", 0.0), ("uq", "uq", SIG), ("oc0", "oc", 0.0), ("oc", "oc", SIG)):
    v = [abs(r["dga_closed"]) for r in R if r["kind"] == kind and r["sig"] == sig and r["model"] in GM]
    desc["mean_abs_over_u_" + name] = float(np.mean(v))
# OC at u=1 numeric and near-offset gain
for name, kind, sig in (("oc0", "oc", 0.0), ("oc", "oc", SIG)):
    rows = [r for r in pick(kind=kind, sig=sig, u=1.0) if r["model"] in GM]
    desc[name + "_absdGA_u1_num"] = med([abs(r["dga_num"]) for r in rows])
    base = {key(r): r for r in pick(kind="uq", sig=sig, u=1.0) if r["model"] in GM}
    desc[name + "_near_gain_ratio"] = med([r["gam_num"][0] / base[key(r)]["gam_num"][0] for r in rows])
out["descriptive"] = desc
(ROOT / "results" / "gate1_verdict.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
