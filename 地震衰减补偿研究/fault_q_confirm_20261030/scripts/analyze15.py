"""Verdicts and descriptive tables of experiment 15 (PLAN_ZH.md section 5) -> analysis/summary.json."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "test"
T90 = 1.7958848187036696                                   # t(0.95, 11)
METHODS = ["Input", "IQ30", "IQ70", "SCa_SR", "VDB", "QSCa", "QSCa_tuned", "SCa_oracle", "IQ70_oracle"]


def recs(var):
    out = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(RES.glob(f"{var}*.json"))]
    return sorted(out, key=lambda r: (r["structure"], r["c0"], r["r"]))


def tint(d):
    d = np.asarray(d, float); m = d.mean(); h = T90 * d.std(ddof=1) / np.sqrt(len(d))
    return dict(mean=float(m), lo=float(m - h), hi=float(m + h), n=len(d), positive=int((d > 0).sum()))


def paired(R, a, b, key="S"):
    d = [r["scores"][a][key] - r["scores"][b][key] for r in R]
    out = tint(d); out["by_structure"] = {k: float(np.mean([x for x, r in zip(d, R) if r["structure"] == k])) for k in "BC"}
    return out


def verdict_pos(*ints, structures=True):
    if all(i["lo"] > 0 and (not structures or all(v > 0 for v in i["by_structure"].values())) for i in ints):
        return "支持"
    if any(i["hi"] < 0 for i in ints):
        return "相反"
    return "无法判定"


def fit_slope(R, rng, nboot=2000):
    u = np.array([r["theory"]["u_qs"] for r in R]); v = 1 / np.array([r["qscan"] for r in R])
    A, C = np.polyfit(u, v, 1); bs = []
    for _ in range(nboot):
        i = rng.integers(0, len(R), len(R))
        if np.ptp(u[i]) > 0:
            bs.append(np.polyfit(u[i], v[i], 1)[0])
    lo, hi = np.percentile(bs, [5, 95]); return dict(A=float(A), C=float(C), lo=float(lo), hi=float(hi))


def qerr(R):
    out = {}
    for ref in ("c_eff_cum", "c_eff_dif"):
        es = [abs(np.log(r["c_sr"] / r[ref])) for r in R]; eq = [abs(np.log(r["qscan"] / r[ref])) for r in R]
        out[ref] = dict(SR=float(np.mean(es)), QScan=float(np.mean(eq)), diff=tint(np.subtract(es, eq)))
    return out


def group(R, rng):
    s = dict(n=len(R), mean_S={m: float(np.mean([r["scores"][m]["S"] for r in R])) for m in METHODS},
             mean_M1={m: float(np.mean([r["scores"][m]["M1"] for r in R])) for m in METHODS})
    s["H1"] = paired(R, "QSCa_tuned", "QSCa")
    s["H2"] = dict(SCa_SR=paired(R, "SCa_SR", "QSCa"), IQ70=paired(R, "IQ70", "QSCa"))
    s["H3"] = dict(IQ70=paired(R, "IQ70", "Input"), SCa_SR=paired(R, "SCa_SR", "Input"))
    s["H4"] = fit_slope(R, rng); s["H5"] = qerr(R)
    s["IQ70_vs_SCa_SR"] = paired(R, "IQ70", "SCa_SR")
    s["oracle_gap"] = dict(SCa=paired(R, "SCa_oracle", "SCa_SR"), IQ70=paired(R, "IQ70_oracle", "IQ70"))
    s["local_mean"] = {m: {k: float(np.mean([r["local"][m][k] for r in R])) for k in R[0]["local"][m]} for m in METHODS}
    s["qscan_boundary"] = int(sum(r["qscan_boundary"] for r in R))
    s["selections"] = [dict(id=r["id"], c_eff_cum=r["c_eff_cum"], c_eff_dif=r["c_eff_dif"], u_qs=r["theory"]["u_qs"],
                            c_sr=r["c_sr"], qscan=r["qscan"]) for r in R]
    return s


def main():
    rng = np.random.default_rng(20261030); Z, Y = recs("Z"), recs("Y")
    out = dict(Z=group(Z, rng), Y=group(Y, rng) if Y else None)
    z = out["Z"]
    out["verdicts"] = dict(
        H1=verdict_pos(z["H1"]), H2=verdict_pos(z["H2"]["SCa_SR"], z["H2"]["IQ70"]),
        H3=verdict_pos(z["H3"]["IQ70"], z["H3"]["SCa_SR"], structures=False),
        H4=("支持" if abs(z["H4"]["A"] - 0.33) <= 0.2 and z["H4"]["hi"] < 0.75 else "相反" if z["H4"]["lo"] > 0.75 else "无法判定"),
        H5=("支持" if z["H5"]["c_eff_cum"]["diff"]["hi"] < 0 else "相反" if z["H5"]["c_eff_cum"]["diff"]["lo"] > 0 else "无法判定"))
    if Y:
        ym = {r["id"][1:]: r for r in Y}
        out["cost_Z_minus_Y"] = {m: tint([zr["scores"][m]["S"] - ym[zr["id"][1:]]["scores"][m]["S"] for zr in Z]) for m in METHODS}
        out["selection_diff"] = [dict(id=zr["id"][1:], obs=1 / zr["qscan"] - 1 / ym[zr["id"][1:]]["qscan"],
                                      pred=(1 - 0.675) * (zr["theory"]["u_qs"] - ym[zr["id"][1:]]["theory"]["u_qs"]))
                                 for zr in Z]
    (ROOT / "analysis").mkdir(exist_ok=True)
    (ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(out["verdicts"], ensure_ascii=False))


if __name__ == "__main__":
    main()
