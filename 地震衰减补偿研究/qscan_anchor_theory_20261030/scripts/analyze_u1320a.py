"""Score PLAN_U1320A_ZH.md from results18/ -> analysis/u1320a.json."""
import json, glob
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; PRED = {"w2.24": 0.218, "w3.00": 0.884}
R = {}
for f in glob.glob(str(ROOT / "results18" / "*.json")):
    r = json.loads(Path(f).read_text()); R.setdefault((r["line"], r["window"]), {})[r["v_anchor"]] = r
out = {}
for (ln, w), d in sorted(R.items()):
    v = {a: d[a]["v_sel"] for a in d}; ok = all(not d[a]["edge"] and d[a]["ok"] for a in d)
    out[f"{ln}|{w}"] = dict(v=v, kappa=(v[2.0] - v[0.5]) / 1.5, pred=PRED[w], reliable=ok, monotone=bool(v[0.5] < v[1.0] < v[2.0]))
lines = sorted({k.split("|")[0] for k in out})
def crit(w, tol=0.15):
    rs = [out[f"{l}|{w}"] for l in lines if out[f"{l}|{w}"]["reliable"]]
    return "undecided" if not rs else bool(len(rs) == len(lines) and all(abs(r["kappa"] - PRED[w]) <= tol for r in rs))
V = {"U1": crit("w2.24"), "U2": crit("w3.00"),
     "U3": bool(all(out[f"{l}|w3.00"]["kappa"] > out[f"{l}|w2.24"]["kappa"] for l in lines)),
     "U4": sum(r["monotone"] for r in out.values()) >= 3}
out["verdicts"] = V; (ROOT / "analysis" / "u1320a.json").write_text(json.dumps(out, indent=1))
for k, r in out.items():
    print(k, r if k == "verdicts" else {x: (round(y, 3) if isinstance(y, float) else y) for x, y in r.items() if x != "v"})
