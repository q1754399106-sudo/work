"""Verdicts of PLAN_CONFIRM_SHAPING_ZH.md section 5-6 from results/cases/*.json -> analysis/summary.json (and a printed summary).

Unit of resampling: structure (6). Per case, a difference is averaged over the m set first (mid = 0.85, 1, 1.2; wide = 0.6 ... 2.0),
then averaged within each structure (d_s), then over the 6 structures. 90% intervals: t (df 5, t = 2.015) and a structure
bootstrap (10 000 draws, seed 20261119, 5-95% percentiles); decisions use the lower of the two lower bounds and the higher of
the two upper bounds. Failure rate: share of cases with d < -0.02. Worst m: min over the m set of the structure-mean difference.
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; CASES = ROOT / "results" / "cases"; OUT = ROOT / "analysis"
MID = (0.85, 1.0, 1.2); WIDE = (0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); ALL = (0.5,) + WIDE
T90 = {5: 2.015}; FAIL = -0.02; METHODS = ("VDB", "IQ70", "UQW", "UQW_TI", "PH_TI", "HAW", "UQW_TIw")


def load():
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(CASES.glob("*.json"))]
    tab = {}
    for r in recs:
        for row in r["rows"]:
            tab[(r["case"], row["method"], row["m"])] = row
    return recs, tab


def stats(recs, tab, a, b, metric, ms, sel):
    cs = [r for r in recs if sel(r)]
    if not cs:
        return None
    d = {r["case"]: float(np.mean([tab[(r["case"], a, m)][metric] - tab[(r["case"], b, m)][metric] for m in ms])) for r in cs}
    by = {}
    for r in cs:
        by.setdefault(r["structure"], []).append(r["case"])
    S = sorted(by); ds = np.array([np.mean([d[c] for c in by[s]]) for s in S]); n = len(ds); mu = float(ds.mean())
    se = float(ds.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan"); t = T90.get(n - 1, 2.015)
    rng = np.random.default_rng(20261119); boot = np.array([rng.choice(ds, n).mean() for _ in range(10000)])
    blo, bhi = (float(v) for v in np.percentile(boot, [5, 95]))
    per_m = {m: float(np.mean([np.mean([tab[(c, a, m)][metric] - tab[(c, b, m)][metric] for c in by[s]]) for s in S])) for m in ms}
    return dict(mean=mu, median=float(np.median(ds)), worst_structure=float(ds.min()), worst_case=float(min(d.values())),
                fail_rate=float(np.mean([v < FAIL for v in d.values()])), n_cases=len(d), per_structure=dict(zip(S, ds.tolist())),
                t90=[mu - t * se, mu + t * se], boot90=[blo, bhi], lo=min(mu - t * se, blo), hi=max(mu + t * se, bhi),
                per_m=per_m, worst_m=min(per_m.values()), positive_structures=int((ds > 0).sum()))


def main():
    recs, tab = load(); OUT.mkdir(exist_ok=True)
    bad = [r["case"] for r in recs if r["nonfinite"]]
    P = lambda r: r["state"] == "cal" and r["variant"] in ("TV", "ST")  # noqa: E731
    st = lambda *a: stats(recs, tab, *a)  # noqa: E731
    g = {}
    g["G0"] = st("UQW", "VDB", "T2", MID, P); g["G0"]["pass"] = g["G0"]["mean"] < 0 and g["G0"]["hi"] < 0
    for k, ms in (("G1", MID), ("G2", WIDE)):
        s = st("UQW_TI", "VDB", "T2", ms, P); s["pass"] = s["mean"] >= 0 and s["lo"] >= -0.02 and s["fail_rate"] <= 0.25
        if k == "G2":
            s["pass"] = s["pass"] and s["worst_m"] >= -0.03
        g[k] = s
    s = st("UQW_TI", "PH_TI", "T2", WIDE, P); s["pass"] = s["mean"] > 0 and s["lo"] > 0
    s["statement"] = ("UQ 部分有贡献" if s["pass"] else "UQ 部分没有贡献" if s["mean"] <= 0 else "UQ 部分的贡献未得到确认"); g["G3"] = s
    g["G3_mid"] = st("UQW_TI", "PH_TI", "T2", MID, P)
    r = 1 - g["G1"]["mean"] / g["G0"]["mean"] if g["G0"]["mean"] != 0 else float("nan")
    g["G4"] = dict(r=r, applies=g["G0"]["pass"], pass_=bool(g["G0"]["pass"] and r >= 0.6))
    s = st("UQW_TI", "VDB", "S", WIDE, P); s["pass"] = s["mean"] >= -0.02; g["B1"] = s
    s = st("UQW_TI", "UQW", "T2", MID, lambda r: r["state"] == "cal" and r["variant"] == "K0"); s["pass"] = s["mean"] >= -0.03; g["K0"] = s
    ND = lambda r: r["state"] == "nd"  # noqa: E731
    s = st("UQW_TI", "VDB", "T2", WIDE, ND); s["pass"] = s["mean"] > 0 and s["lo"] > 0; g["ND1"] = s
    s = st("UQW_TI", "UQW", "T2", WIDE, ND); s["pass"] = s["mean"] >= -0.05; g["ND2"] = s

    h1 = ("无法判定（新数据未复现差距）" if not g["G0"]["pass"] else "得到支持（半合成）" if g["G1"]["pass"] and g["G4"]["pass_"] else "不支持")
    main_ok = g["G1"]["pass"] and g["G2"]["pass"]; guards = {k: g[k]["pass"] for k in ("K0", "ND1", "ND2")}
    if bad:
        label = "无法判定（运行异常）"
    elif not main_ok:
        label = "已淘汰"
    elif not g["G3"]["pass"]:
        label = "机制研究（整形有效，UQ 无贡献或未确认）"
    elif not g["B1"]["pass"]:
        label = "无法判定（T2 依赖）"
    elif not all(guards.values()):
        label = "方法候选（限定条件：" + "、".join(k for k, v in guards.items() if not v) + " 不过）"
    else:
        label = "方法候选（合成确认通过，实测未确认）"

    desc = {}
    for state in ("cal", "nd"):
        for var in ("K0", "TV", "ST"):
            sel = lambda r, s_=state, v_=var: r["state"] == s_ and r["variant"] == v_  # noqa: E731
            for meth in ("UQW", "UQW_TI", "PH_TI", "HAW", "IQ70", "UQW_TIw"):
                for metric in ("T2", "S", "T2mis"):
                    for tag, ms in (("mid", MID), ("wide", WIDE)):
                        desc[f"{state}/{var}/{meth}-VDB/{metric}/{tag}"] = st(meth, "VDB", metric, ms, sel)
    absval = {}
    for state in ("cal", "nd"):
        for var in ("K0", "TV", "ST"):
            cs = [r["case"] for r in recs if r["state"] == state and r["variant"] == var]
            for meth in METHODS:
                for metric in ("T2", "S", "latcoh", "hf"):
                    absval[f"{state}/{var}/{meth}/{metric}"] = {m: float(np.mean([tab[(c, meth, m)][metric] for c in cs])) for m in ALL}
            for metric in ("T2", "S", "latcoh", "hf"):
                absval[f"{state}/{var}/Input/{metric}"] = float(np.mean([tab[(c, "Input", 0.0)][metric] for c in cs]))
    res = dict(n_cases=len(recs), nonfinite_cases=bad, gates=g, H1=h1, label=label, guards=guards, describe=desc, absolute=absval)
    (OUT / "summary.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    print("cases", len(recs), "nonfinite", len(bad))
    for k in ("G0", "G1", "G2", "G3", "G3_mid", "B1", "K0", "ND1", "ND2"):
        s = g[k]
        print("%-6s mean %+.4f median %+.4f worst_s %+.4f worst_case %+.4f fail %.2f worst_m %+.4f t90 [%+.4f %+.4f] boot90 [%+.4f %+.4f] pos %d/6 pass %s" % (
            k, s["mean"], s["median"], s["worst_structure"], s["worst_case"], s["fail_rate"], s["worst_m"], *s["t90"], *s["boot90"],
            s["positive_structures"], s.get("pass")))
    print("G4 r %.3f applies %s pass %s" % (g["G4"]["r"], g["G4"]["applies"], g["G4"]["pass_"]))
    print("G3:", g["G3"]["statement"]); print("H1:", h1); print("label:", label)


if __name__ == "__main__":
    main()
