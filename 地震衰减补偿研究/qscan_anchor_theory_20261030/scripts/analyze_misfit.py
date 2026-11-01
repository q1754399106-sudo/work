"""Analysis of PLAN_MISFIT_ZH.md: paired change of the lambda -> 0 beta per variant. -> analysis/misfit.json"""
import json, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results24"
IDS = [f"c{c}_s{s}" for s in (9971, 9972) for c in (0.7, 1.0, 1.4)]


def b0(b16, b2):
    return b16 - (b2 - b16) / 7.0


def row(var, did):
    r = json.loads((R / f"{var}_{did}.json").read_text()); u = r["u"]
    if var == "R0":
        v = {lam: json.loads((ROOT / "results23" / f"S_{did}_lam{lam}.json").read_text()) for lam in (0.0625, 0.5)}
        v = {str(k): x for k, x in v.items()}
    else:
        v = r["scans"]
    edge = any(v[k]["edge"] for k in ("0.0625", "0.5"))
    bl, bh = v["0.0625"]["v_sel"] - u, v["0.5"]["v_sel"] - u
    return dict(edge=edge, beta16=bl, beta_half=bh, beta0=b0(bl, bh), sr_bias=1 / r["c_sr"] - u,
                beta_sr=v["0.5"]["v_sel"] - 1 / r["c_sr"], check=r.get("check_max_abs_diff_rel"))


out = {}
base = {d: row("R0", d) for d in IDS}; out["R0"] = dict(rows=base)
for var in ("P1", "P2", "P3", "M1"):
    rows = {d: row(var, d) for d in IDS}; ok = [d for d in IDS if not rows[d]["edge"] and not base[d]["edge"]]
    dif = lambda k: st.median(rows[d][k] - base[d][k] for d in ok) if ok else None
    db0 = dif("beta0"); cls = None if db0 is None else ("main" if db0 >= 0.15 else "partial" if db0 >= 0.05 else "not")
    out[var] = dict(rows=rows, n=len(ok), d_beta0=db0, d_beta_half=dif("beta_half"), d_sr_bias=dif("sr_bias"),
                    d_beta_sr=dif("beta_sr"), verdict=cls, sr_route=(dif("beta_sr") is not None and dif("beta_sr") >= 0.15 and cls != "main"))
(ROOT / "analysis" / "misfit.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("R0 check max", max(x["check"] for x in base.values()), " beta0 median %.3f  sr_bias median %.3f" %
      (st.median(x["beta0"] for x in base.values()), st.median(x["sr_bias"] for x in base.values())))
for var in ("P1", "P2", "P3", "M1"):
    x = out[var]; print(var, {k: (round(x[k], 3) if isinstance(x[k], float) else x[k]) for k in ("n", "d_beta0", "d_beta_half", "d_sr_bias", "d_beta_sr", "verdict", "sr_route")})
