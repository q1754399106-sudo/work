"""Summary of results/dev10 (and the FD part): per state, mean by m, mean, 90% bootstrap (datasets) of differences."""
import json, glob, numpy as np
MS = [0.5, 0.7, 1.0, 1.4, 2.0]


def boot(d, n=4000):
    g = np.random.default_rng(1); return np.percentile(d[g.integers(0, len(d), (n, len(d)))].mean(1), [5, 95])


for pre, keys in (("semi", ("S", "T2")), ("fd", ("S",))):
    fs = sorted(glob.glob(f"../results/dev10/{pre}_*.json")); D = {}
    for i, p in enumerate(fs):
        for r in json.load(open(p)):
            if r["method"].startswith("_"):
                continue
            for k in keys:
                D.setdefault((r["state"], k, r["method"]), np.full((len(fs), 5), np.nan))[i, MS.index(r["m"])] = r[k]
    print(pre, "units", len(fs))
    for st in dict.fromkeys(s for s, _, _ in D):
        for k in keys:
            print("== %s %s" % (st, k))
            for meth in ("UQW", "VDBin", "HWin", "UQWS"):
                a = D[(st, k, meth)]; line = "  %-6s " % meth + " ".join("%.3f" % x for x in a.mean(0)) + " | mean %.3f worst-m %.3f" % (a.mean(), a.mean(0).min())
                for ref in ("UQW", "VDBin"):
                    if meth != ref and meth != "VDBin":
                        d = (a - D[(st, k, ref)]).mean(1); lo, hi = boot(d); line += "  vs %s %+.3f [%+.3f,%+.3f] %d/%d" % (ref, d.mean(), lo, hi, (d > 0).sum(), len(d))
                if meth == "UQWS":
                    d = a - D[(st, k, "UQW")]; line += "\n         UQWS-UQW by m " + " ".join("%+.3f" % x for x in d.mean(0))
                print(line)
