"""Verdicts V1-V6 of PLAN_CONFIRM_SC_ZH.md from results/confirmSC (and UQW-S of results/confirmS for FD).
Semi: unit = dataset (12), 90% interval from 2000 bootstrap resamples of datasets (seed 20261198). FD: unit = model (6),
realisations averaged first; mean and number of positive models.  usage: analyze_confirmSC.py [--dev DIR]"""
import sys, json, glob
import numpy as np
from env import ROOT

M5 = ("m0.5", "m0.7", "m1", "m1.4", "m2"); OVER = ("m0.5", "m0.7")
SRC = sys.argv[sys.argv.index("--dev") + 1] if "--dev" in sys.argv else str(ROOT / "results" / "confirmSC")
CSRC = sys.argv[sys.argv.index("--dev") + 2] if "--dev" in sys.argv else str(ROOT / "results" / "confirmS")


def boot(d):
    g = np.random.default_rng(20261198); return [float(x) for x in np.percentile(d[g.integers(0, len(d), (2000, len(d)))].mean(1), [5, 95])]


semi = {}
for p in sorted(glob.glob(f"{SRC}/S6_*.json")):
    did, st = p.split("S6_")[1][:-5].rsplit("_", 1)
    for r in json.load(open(p))["records"]:
        semi.setdefault(st, {}).setdefault(did, {})[r["point"]] = r


def sdiff(st, pts, a="UQWS", b="SCa"):
    ids = sorted(semi[st]); d = np.array([np.mean([semi[st][i][p][a] - semi[st][i][p][b] for p in pts]) for i in ids])
    return dict(mean=float(d.mean()), ci=boot(d), pos=int((d > 0).sum()), n=len(d))


fdsc, fduq = {}, {}
for p in glob.glob(f"{SRC}/FD3_*.json"):
    h, lv = p.split("FD3_")[1].split("_")[:2]
    for r in json.load(open(p))["records"]:
        fdsc.setdefault((lv, h, r["point"]), []).append(r["SC"])
for p in glob.glob(f"{CSRC}/FD3_*.json"):
    h, lv = p.split("FD3_")[1].split("_")[:2]
    for r in json.load(open(p))["records"]:
        if r["method"] == "UQWS" and f"m{r['m']:g}" in M5:
            fduq.setdefault((lv, h, f"m{r['m']:g}"), []).append(r["S"])


def fdiff(lv, pts):
    hs = sorted({k[1] for k in fdsc if k[0] == lv}, key=lambda h: int(h[1:]))
    d = np.array([np.mean([np.mean(fduq[(lv, h, p)]) - np.mean(fdsc[(lv, h, p)]) for p in pts]) for h in hs])
    return dict(mean=float(d.mean()), pos=int((d > 0).sum()), n=len(d), by_model=d.tolist())


R = {}
R["V1"] = dict(test="semi res4, SR route: UQWS - SCa", **sdiff("res4", ("SR",))); R["V1"]["holds"] = R["V1"]["ci"][0] > 0
R["V2"] = dict(test="semi gau4, SR route: UQWS - SCa", **sdiff("gau4", ("SR",))); R["V2"]["holds"] = R["V2"]["ci"][0] > 0
R["V3"] = dict(test="semi res4, 5 m: UQWS - SCa", **sdiff("res4", M5)); R["V3"]["holds"] = R["V3"]["ci"][0] > 0
R["V4"] = dict(test="semi cal, SR route: UQWS - SCa", **sdiff("cal", ("SR",))); R["V4"]["holds"] = R["V4"]["ci"][1] < 0
R["V5"] = dict(test="semi cal, m 1: UQWS - SCa", **sdiff("cal", ("m1",))); R["V5"]["holds"] = R["V5"]["ci"][1] < 0
R["V6_semi"] = dict(test="semi cal, m 0.5-0.7: UQWS - SCa", **sdiff("cal", OVER)); R["V6_semi"]["holds"] = R["V6_semi"]["ci"][0] > 0
for lv in ("30dB", "14dB"):
    k = "V6_fd" + lv; R[k] = dict(test=f"FD {lv}, m 0.5-0.7: UQWS - SC", **fdiff(lv, OVER)); R[k]["holds"] = R[k]["mean"] > 0 and R[k]["pos"] >= 5
desc = {}
for st in semi:
    desc[f"semi_{st}"] = {p: sdiff(st, (p,)) for p in M5 + ("SR",)}; desc[f"semi_{st}"]["5m"] = sdiff(st, M5)
    desc[f"semi_{st}"]["SCa_conv"] = int(sum(r["SCa_conv"] for i in semi[st] for r in semi[st][i].values()))
    desc[f"semi_{st}"]["mean"] = {a: {p: float(np.mean([semi[st][i][p][a] for i in semi[st]])) for p in M5 + ("SR",)} for a in ("UQWS", "SCa")}
for lv in ("30dB", "14dB"):
    if any(k[0] == lv for k in fdsc):
        desc[f"FD_{lv}"] = {p: fdiff(lv, (p,)) for p in M5}; desc[f"FD_{lv}"]["5m"] = fdiff(lv, M5)
out = dict(rules=R, confirmed=all(v["holds"] for v in R.values()), describe=desc)
if "--dev" not in sys.argv:
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "confirmSC.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for k, v in R.items():
    ci = " [%+.3f, %+.3f]" % tuple(v["ci"]) if "ci" in v else ""
    print("%-9s %-38s %+.3f%s %d/%d %s" % (k, v["test"], v["mean"], ci, v["pos"], v["n"], "holds" if v["holds"] else "FAILS"))
print("confirmed:", out["confirmed"])
for k, D in desc.items():
    print(k, "  ".join("%s %+.3f(%d)" % (p, D[p]["mean"], D[p]["pos"]) for p in M5 + (("SR",) if k.startswith("semi") else ())), "| 5m %+.3f" % D["5m"]["mean"])
