"""Summary of the WR development grid (analysis/dev/*.json) -> analysis/dev_summary.json."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq  # noqa: E402
from wq import np  # noqa: E402

recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((wq.ROOT / "analysis" / "dev").glob("*.json"))]
out = {}
for var in ("sm5", "noRC", "SR"):
    rows = []
    for r in recs:
        for j in range(3):
            e = r["est"][f"n{j}_{var}"]
            rows.append(dict(well=r["well"], c=r["c_true"], err=float(np.log(e["c"] / r["c_true"])),
                             cover=bool(e["ci90"][0] <= r["c_true"] <= e["ci90"][1]) if "ci90" in e else None,
                             width=float(e["ci90"][1] / e["ci90"][0]) if "ci90" in e else None))
    err = np.array([x["err"] for x in rows]); cov = [x["cover"] for x in rows if x["cover"] is not None]
    wid = [x["width"] for x in rows if x["width"] is not None]
    s = dict(n=len(rows), median_err=float(np.median(err)), mean_abs_log_err=float(np.mean(np.abs(err))),
             within25=float(np.mean(np.abs(err) <= np.log(1.25))), coverage=float(np.mean(cov)) if cov else None,
             median_ci_ratio=float(np.median(wid)) if wid else None, by={})
    for w in ("F-12", "F-14"):
        for c in sorted({x["c"] for x in rows}):
            e = np.array([x["err"] for x in rows if x["well"] == w and x["c"] == c])
            s["by"][f"{w}_c{c}"] = dict(median_ratio=float(np.exp(np.median(e))), mean_abs=float(np.mean(np.abs(e))))
    out[var] = s
    print(var, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items() if k != "by"})
    print("   ", {k: round(v["median_ratio"], 2) for k, v in s["by"].items()})
(wq.ROOT / "analysis" / "dev_summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
