"""Verdicts of the confirmation stage (PLAN_ZH.md section 5) -> analysis/summary.json."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "test"
T90 = 1.7958848187036696                                   # t(0.95, 11)
METHODS = ["Input", "IQ70", "VDB", "SR_mu10", "SR_sel", "QSCa", "QT_mu10", "QT_sel"]
STRUCTS = "GHI"


def tint(d):
    d = np.asarray(d, float); m = d.mean(); h = T90 * d.std(ddof=1) / np.sqrt(len(d))
    return dict(mean=float(m), lo=float(m - h), hi=float(m + h), n=len(d), positive=int((d > 0).sum()))


def paired(R, a, b):
    d = [r["scores"][a]["S"] - r["scores"][b]["S"] for r in R]; out = tint(d)
    out["by_structure"] = {k: float(np.mean([x for x, r in zip(d, R) if r["structure"] == k])) for k in STRUCTS}
    out["by_c0"] = {str(c): float(np.mean([x for x, r in zip(d, R) if r["c0"] == c])) for c in (0.75, 1.05)}
    return out


def verdict(i):
    if i["lo"] > 0 and all(v > 0 for v in i["by_structure"].values()):
        return "支持"
    return "相反" if i["hi"] < 0 else "无法判定"


def main():
    R = sorted((json.loads(p.read_text(encoding="utf-8")) for p in RES.glob("Z*.json")), key=lambda r: (r["structure"], r["c0"], r["r"]))
    choice = R[0]["choice"]
    s = dict(n=len(R), choice=choice, mean_S={m: float(np.mean([r["scores"][m]["S"] for r in R])) for m in METHODS})
    s["C1"] = paired(R, "SR_sel", "SR_mu10"); s["C2"] = paired(R, "QT_sel", "QT_mu10")
    s["C3"] = paired(R, "SR_sel", "IQ70"); s["C4"] = paired(R, "QT_sel", "QSCa")
    s["verdicts"] = {h: ("不适用（选中 10）" if h in ("C1", "C2") and choice["SR" if h == "C1" else "QT"] == 10 else verdict(s[h]))
                     for h in ("C1", "C2", "C3", "C4")}
    s["descriptive"] = dict(QT_mu10_vs_QSCa=paired(R, "QT_mu10", "QSCa"), QT_sel_vs_IQ70=paired(R, "QT_sel", "IQ70"),
                            SR_sel_vs_VDB=paired(R, "SR_sel", "VDB"), SR_mu10_vs_IQ70=paired(R, "SR_mu10", "IQ70"))
    s["per_dataset"] = [dict(id=r["id"], c_eff_cum=r["c_eff_cum"], c_sr=r["c_sr"], qscan=r["qscan"], fs=r["fs"], fd=r["fd"],
                             level=r["level"], S={m: r["scores"][m]["S"] for m in METHODS}) for r in R]
    s["fs_collapse"] = int(sum(r["fs"] < r["fd"] for r in R)); s["qscan_boundary"] = int(sum(r["qscan_boundary"] for r in R))
    (ROOT / "analysis").mkdir(exist_ok=True)
    (ROOT / "analysis" / "summary.json").write_text(json.dumps(s, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(s["verdicts"], ensure_ascii=False))


if __name__ == "__main__":
    main()
