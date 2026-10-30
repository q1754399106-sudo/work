"""Verdicts and descriptive tables of experiment 16 (PLAN_ZH.md section 5) -> analysis/summary.json."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "test"
T90 = 1.7958848187036696                                   # t(0.95, 11)
METHODS = ["Input", "IQ70", "VDB", "SR_mu10", "SR_mu1", "SR_mu0", "QSCa", "QSCa_mu0", "QT_mu10", "QT_mu0", "OR_mu10", "OR_mu0"]
STRUCTS = "DEF"


def recs():
    out = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(RES.glob("Z*.json"))]
    return sorted(out, key=lambda r: (r["structure"], r["c0"], r["r"]))


def tint(d):
    d = np.asarray(d, float); m = d.mean(); h = T90 * d.std(ddof=1) / np.sqrt(len(d))
    return dict(mean=float(m), lo=float(m - h), hi=float(m + h), n=len(d), positive=int((d > 0).sum()))


def paired(R, a, b):
    d = [r["scores"][a]["S"] - r["scores"][b]["S"] for r in R]; out = tint(d)
    out["by_structure"] = {k: float(np.mean([x for x, r in zip(d, R) if r["structure"] == k])) for k in STRUCTS}
    out["by_c0"] = {str(c): float(np.mean([x for x, r in zip(d, R) if r["c0"] == c])) for c in (0.75, 1.05)}
    return out


def one_sided(i):
    """支持: lower bound > 0 and every structure mean > 0; 相反: upper bound < 0; else 无法判定."""
    if i["lo"] > 0 and all(v > 0 for v in i["by_structure"].values()):
        return "支持"
    return "相反" if i["hi"] < 0 else "无法判定"


def two_sided(i):
    return "SC-a 更好" if i["lo"] > 0 else "反 Q 更好" if i["hi"] < 0 else "无法判定"


def main():
    R = recs()
    s = dict(n=len(R), mean_S={m: float(np.mean([r["scores"][m]["S"] for r in R])) for m in METHODS},
             mean_M1={m: float(np.mean([r["scores"][m]["M1"] for r in R])) for m in METHODS})
    s["H1"] = paired(R, "SR_mu0", "SR_mu10"); s["H2"] = paired(R, "QT_mu0", "QT_mu10")
    s["H3"] = paired(R, "OR_mu0", "OR_mu10"); s["H4"] = paired(R, "QT_mu0", "QSCa"); s["H5"] = paired(R, "SR_mu0", "IQ70")
    s["verdicts"] = dict(H1=one_sided(s["H1"]), H2=one_sided(s["H2"]), H3=one_sided(s["H3"]), H4=one_sided(s["H4"]),
                         H5=two_sided(s["H5"]))
    s["descriptive"] = dict(SR_mu1_vs_mu10=paired(R, "SR_mu1", "SR_mu10"), SR_mu0_vs_mu1=paired(R, "SR_mu0", "SR_mu1"),
                            QSCa_mu0_vs_QSCa=paired(R, "QSCa_mu0", "QSCa"), QT_mu10_vs_QSCa=paired(R, "QT_mu10", "QSCa"),
                            QT_mu0_vs_IQ70=paired(R, "QT_mu0", "IQ70"), OR_mu0_vs_IQ70=paired(R, "OR_mu0", "IQ70"),
                            VDB_vs_IQ70=paired(R, "VDB", "IQ70"))
    s["local_mean"] = {m: {k: float(np.mean([r["local"][m][k] for r in R])) for k in R[0]["local"][m]} for m in METHODS}
    s["per_dataset"] = [dict(id=r["id"], c_eff_cum=r["c_eff_cum"], c_sr=r["c_sr"], qscan=r["qscan"], boundary=r["qscan_boundary"],
                             fs=r["fs"], fs_last=r["fs_last"], fd=r["fd"], f_true=r["f_true"], level=r["level"],
                             S={m: r["scores"][m]["S"] for m in METHODS}) for r in R]
    s["fs_collapse"] = int(sum(r["fs"] < r["fd"] for r in R))
    (ROOT / "analysis").mkdir(exist_ok=True)
    (ROOT / "analysis" / "summary.json").write_text(json.dumps(s, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(s["verdicts"], ensure_ascii=False))


if __name__ == "__main__":
    main()
