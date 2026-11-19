"""Analysis of PLAN_ZH.md: per-scenario CSVs, reproduction check against archived stages, model-stratified summaries,
rules N and L (with the validity conditions) -> analysis/per_scenario_<method>.csv, analysis/summary.json, printed summary.

Model = semi-synthetic structure A-F (4 datasets each) or FD H1-H6 (3 realisations each, per noise level). Each model is
first averaged over its datasets / realisations and the 5 m; mean, median, worst (min) and number positive over the 6.
"""
import sys, json, glob, csv, collections
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; QC = ROOT.parent / "qunknown_comp_20261109" / "results"
MS = (0.5, 0.7, 1.0, 1.4, 2.0); CLASSES = ("semi", "14dB", "08dB")
METHODS = ("V0", "V1", "V2", "V3", "V0s25", "V1s25", "V0d", "W0", "VDB", "IQ70")
TOL_REPRO, TOL_N, TOL_L_DROP, TOL_QUAD = 1e-6, 0.01, 0.03, 0.005
DIFFS = [("V1", "V0"), ("V0", "V2"), ("V3", "V0"), ("V0d", "V0"), ("V1s25", "V0s25"),                     # rules / validity
         ("V0", "VDB"), ("V0", "IQ70"), ("V2", "W0"), ("V2", "VDB"), ("V2", "IQ70"), ("V3", "V2"), ("V0", "V3"),
         ("V0", "V1"), ("V1", "VDB"), ("V1", "IQ70"), ("V0s25", "V1s25"), ("V0", "W0"), ("V3", "VDB"), ("V0s25", "VDB")]


def load():
    rows = []
    for p in sorted(glob.glob(str(ROOT / "results" / "raw" / "*.json"))):
        J = json.loads(Path(p).read_text(encoding="utf-8")); cls = "semi" if J["kind"] == "semi" else J["level"]
        for r in J["records"]:
            rows.append(dict(job=J["id"], cls=cls, model=J["model"], c0=J["c0"], r=J["r"], level=J["level"], m=r["m"],
                             method=r["method"], S=r["S"], T2=r["T2"], finite=r["finite"]))
    return rows


def write_csvs(rows):
    (ROOT / "analysis").mkdir(exist_ok=True)
    for k in METHODS:
        with open(ROOT / "analysis" / f"per_scenario_{k}.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh); w.writerow(["job", "data_class", "model", "c0", "r", "noise_level", "m", "S", "T2", "finite"])
            for q in rows:
                if q["method"] == k:
                    w.writerow([q["job"], q["cls"], q["model"], q["c0"], q["r"], q["level"], q["m"], "%.10f" % q["S"],
                                "" if q["T2"] is None else "%.10f" % q["T2"], q["finite"]])


def reproduce(rows):
    """max |S_new - S_archived| per method (V0 vs stage W UQW_g0 sig 0.447, W0 vs stage W W0, IQ70 / VDB vs stage 0 / 2b)."""
    arch = {}
    for p in glob.glob(str(QC / "stageW" / "*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8"))
        for r in J["records"]:
            if (r["method"] == "UQW_g0" and abs(r["sig"] - 0.447) < 1e-12) or r["method"] == "W0":
                arch[(J["id"], "V0" if r["method"] == "UQW_g0" else "W0", r["m"])] = r["S"]
    for p in glob.glob(str(QC / "stage0" / "*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8"))
        if J["method"] in ("IQ70", "VDB"):
            arch[(J["id"], J["method"], J["m"])] = J["score"]["S"]
    for p in glob.glob(str(QC / "stage2b" / "*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8")); jid = "FD_%s_%s_r%d" % (J["model"], J["level"], J["r"])
        for r in J["records"]:
            if r["method"] in ("IQ70", "VDB"):
                arch[(jid, r["method"], r["m"])] = r["S"]
    out = collections.defaultdict(lambda: dict(n=0, max_abs=0.0))
    for q in rows:
        key = (q["job"], q["method"], q["m"])
        if key in arch:
            o = out[q["method"]]; o["n"] += 1; o["max_abs"] = max(o["max_abs"], abs(q["S"] - arch[key]))
    return dict(out)


def table(rows, metric="S"):
    """{cls: {method: {model: model mean over datasets and m}}} and per-m {cls: {method: [mean over models]}}."""
    acc = collections.defaultdict(list); accm = collections.defaultdict(list)
    for q in rows:
        if q[metric] is None:
            continue
        acc[(q["cls"], q["method"], q["model"])].append(q[metric]); accm[(q["cls"], q["method"], q["model"], q["m"])].append(q[metric])
    T = collections.defaultdict(lambda: collections.defaultdict(dict)); PM = collections.defaultdict(dict)
    for (c, k, mo), v in acc.items():
        T[c][k][mo] = float(np.mean(v))
    for c in T:
        for k in T[c]:
            mods = sorted(T[c][k]); PM[c][k] = [float(np.mean([np.mean(accm[(c, k, mo, m)]) for mo in mods])) for m in MS]
    return T, PM


def stats(v):
    v = np.asarray(v, float); return dict(mean=float(v.mean()), median=float(np.median(v)), worst=float(v.min()), n_pos=int((v > 0).sum()), n=len(v),
                                           per_model=[round(float(x), 5) for x in v])


def diff(T, c, a, b):
    mods = sorted(set(T[c][a]) & set(T[c][b])); return stats([T[c][a][mo] - T[c][b][mo] for mo in mods]), mods


def rules(T):
    R = {}
    for c in CLASSES:
        dN = diff(T, c, "V1", "V0")[0]["mean"]; dN25 = diff(T, c, "V1s25", "V0s25")[0]["mean"]
        dQ = diff(T, c, "V0d", "V0")[0]["mean"]
        def ruleL(ref):
            drop = diff(T, c, ref, "V2")[0]["mean"]; rec = diff(T, c, "V3", ref)[0]["mean"]
            return dict(ref=ref, drop_ref_minus_V2=drop, V3_minus_ref=rec, holds=bool(drop >= TOL_L_DROP and abs(rec) < TOL_N))
        L = ruleL("V0"); Lq = ruleL("V0d") if abs(dQ) >= TOL_QUAD else None
        if Lq is None:
            l_verdict = "左尾是主因" if L["holds"] else "无法判定"
        else:
            l_verdict = ("左尾是主因" if L["holds"] else "无法判定") if L["holds"] == Lq["holds"] else "无法判定（求积敏感）"
        R[c] = dict(N=dict(d=dN, holds=bool(abs(dN) < TOL_N), verdict="噪声项缩放不重要" if abs(dN) < TOL_N else "无法判定"),
                    N_sig025=dict(d=dN25, holds=bool(abs(dN25) < TOL_N)),
                    quad=dict(d=dQ, triggered=bool(abs(dQ) >= TOL_QUAD)), L=L, L_quadref=Lq, L_verdict=l_verdict)
    allN = all(R[c]["N"]["holds"] for c in CLASSES); anyN = any(R[c]["N"]["holds"] for c in CLASSES)
    sigN_same = all(R[c]["N"]["holds"] == R[c]["N_sig025"]["holds"] for c in CLASSES)
    R["overall"] = dict(
        N="噪声项缩放不重要" if allN else ("无法判定（随数据类别而异）" if anyN else "无法判定"),
        N_sigma_note=None if sigN_same else "规则 N 在 σ = 0.25 与 σ = 0.447 的判读不同：结论只在 σ = 0.447 时成立",
        L_classes=[c for c in CLASSES if R[c]["L_verdict"] == "左尾是主因"],
        change_design_wording=any(R[c]["L_verdict"] == "左尾是主因" for c in CLASSES))
    return R


if __name__ == "__main__":
    rows = load(); write_csvs(rows); T, PM = table(rows); T2, PM2 = table([q for q in rows if q["cls"] == "semi"], "T2")
    jobs = sorted({q["job"] for q in rows})
    out = dict(n_jobs=len(jobs), n_rows=len(rows), n_nonfinite=sum(not q["finite"] for q in rows), reproduce=reproduce(rows))
    out["repro_ok"] = all(v["max_abs"] <= TOL_REPRO for v in out["reproduce"].values()) and len(out["reproduce"]) == 4
    out["levels"] = {c: {k: stats(list(T[c][k][mo] for mo in sorted(T[c][k]))) for k in METHODS} for c in CLASSES}
    out["diffs"] = {c: {f"{a}-{b}": diff(T, c, a, b)[0] for a, b in DIFFS} for c in CLASSES}
    out["diffs_T2_semi"] = {f"{a}-{b}": diff(T2, "semi", a, b)[0] for a, b in DIFFS}
    out["levels_T2_semi"] = {k: stats(list(T2["semi"][k][mo] for mo in sorted(T2["semi"][k]))) for k in METHODS}
    out["per_m"] = {c: PM[c] for c in CLASSES}; out["per_m_T2_semi"] = PM2["semi"]
    out["rules"] = rules(T) if out["repro_ok"] else "not evaluated: reproduction check failed"
    (ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print("jobs", out["n_jobs"], "rows", out["n_rows"], "nonfinite", out["n_nonfinite"]); print("reproduce", json.dumps(out["reproduce"]), "ok", out["repro_ok"])
    for c in CLASSES:
        print("==", c)
        for k in METHODS:
            s = out["levels"][c][k]; print("  %-6s mean %.4f median %.4f worst %.4f | by m %s" % (k, s["mean"], s["median"], s["worst"], " ".join("%.3f" % x for x in PM[c][k])))
        for key, s in out["diffs"][c].items():
            print("  %-12s mean %+.4f median %+.4f worst %+.4f  %d/%d  %s" % (key, s["mean"], s["median"], s["worst"], s["n_pos"], s["n"], s["per_model"]))
    print("== T2 semi")
    for key, s in out["diffs_T2_semi"].items():
        print("  %-12s mean %+.4f median %+.4f worst %+.4f  %d/%d" % (key, s["mean"], s["median"], s["worst"], s["n_pos"], s["n"]))
    print(json.dumps(out["rules"], indent=1, ensure_ascii=False))
