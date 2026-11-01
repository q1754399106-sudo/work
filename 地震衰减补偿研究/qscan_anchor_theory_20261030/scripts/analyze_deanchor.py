"""Analysis of PLAN_DEANCHOR_ZH.md: de-anchored estimates E_fix = (v - k v_a)/(1 - k), E_fix0 (lambda -> 0 intercept).
-> analysis/deanchor.json"""
import json, glob, statistics as st
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; DV = 1 / 0.7 - 1 / 1.4
V = 1 / np.exp(np.linspace(np.log(0.4), np.log(2.5), 15))


def vertex(J):
    J = np.asarray(J, float); j = int(np.argmin(J))
    if j in (0, len(J) - 1):
        return None
    a = np.polyfit(V[j - 1:j + 2], J[j - 1:j + 2], 2)
    return float(np.clip(-a[1] / (2 * a[0]), V[j + 1], V[j - 1])) if a[0] > 0 else float(V[j])


def est(v_half, v0, k, va=1.0):
    return dict(E_cur=v_half, E_fix=(v_half - k * va) / (1 - k), E_fix0=(None if v0 is None else (v0 - k * va) / (1 - k)))


semi = {}
for p in sorted(glob.glob(str(ROOT / "results" / "T_*__a1.0.json"))):
    r = json.loads(Path(p).read_text()); ra = json.loads(Path(p.replace("a1.0", "a0.7")).read_text())
    rb = json.loads(Path(p.replace("a1.0", "a1.4")).read_text())
    va, vb = vertex(ra["diag"]["J"][0]), vertex(rb["diag"]["J"][0]); vh, vl = vertex(r["diag"]["J"][0]), vertex(r["diag"]["J"][1])
    if None in (va, vb, vh):
        semi[r["id"]] = None; continue
    k = (va - vb) / DV; v0 = None if vl is None else vl - (vh - vl) / 7.0
    semi[r["id"]] = dict(u=1 / r["c_true"], kappa=k, v_half=vh, v_l16=vl, v0=v0, **est(vh, v0, k))

field = {}
LAMS = (0.0625, 0.125, 0.25, 0.5)
for p in sorted(glob.glob(str(ROOT / "results11" / "IL_*__w1.4_a1.0.json"))):
    r5 = json.loads(Path(p).read_text()); ln = r5["inline"]
    ra = json.loads(Path(p.replace("a1.0", "a0.7")).read_text()); rb = json.loads(Path(p.replace("a1.0", "a1.4")).read_text())
    k = (ra["v_sel"] - rb["v_sel"]) / DV; pts = []
    for lam in LAMS:
        q = r5 if lam == 0.5 else (json.loads((ROOT / "results23" / f"F_IL_{ln}_lam{lam}.json").read_text())
                                    if (ROOT / "results23" / f"F_IL_{ln}_lam{lam}.json").exists() else None)
        if q is not None and not q["edge"]:
            pts.append((lam, q["v_sel"]))
    v0 = float(np.polyfit(*np.array(pts).T, 1)[1]) if len(pts) >= 3 else None
    field[str(ln)] = dict(kappa=k, c_sr=r5["c_sr"], v_half=r5["v_sel"], points=pts, v0=v0, **est(r5["v_sel"], v0, k))

out = dict(semi=semi, field=field); S = [x for x in semi.values() if x]
summ = {}
for e in ("E_cur", "E_fix", "E_fix0"):
    xs = [x for x in S if x[e] is not None]; u = np.array([x["u"] for x in xs]); y = np.array([x[e] for x in xs])
    slope = float(np.polyfit(u, y, 1)[0]) if len(xs) >= 3 else None
    F = [x[e] for x in field.values() if x[e] is not None]
    summ[e] = dict(semi_n=len(xs), semi_slope=slope, semi_med_abs_err=(float(np.median(np.abs(y - u))) if xs else None),
                   semi_med_err=(float(np.median(y - u)) if xs else None), field_n=len(F),
                   field_c_median=(float(1 / np.median(F)) if F else None), field_v_median=(float(np.median(F)) if F else None),
                   field_v_range=([float(min(F)), float(max(F))] if F else None))
out["summary"] = summ; e0 = summ["E_fix0"]
D1 = None if e0["semi_slope"] is None else bool(abs(e0["semi_slope"] - 1) <= 0.2 and e0["semi_med_abs_err"] <= 0.15)
D2 = None if e0["semi_med_abs_err"] is None else bool(e0["semi_med_abs_err"] < summ["E_cur"]["semi_med_abs_err"])
D3 = None if e0["field_c_median"] is None else bool(1.02 <= e0["field_c_median"] <= 1.29)
out["sr_c_median"] = float(st.median(x["c_sr"] for x in field.values()))
out["verdicts"] = dict(D1=D1, D2=D2, D3=D3)
(ROOT / "analysis" / "deanchor.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for e, x in summ.items():
    print(e, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in x.items()})
print("SR c median %.3f" % out["sr_c_median"]); print(out["verdicts"])
