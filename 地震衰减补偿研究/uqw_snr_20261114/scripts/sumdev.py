"""Summary of a results directory: sumdev.py DIR PREFIX [reference method] -> mean over datasets by m, S and T2."""
import sys, json, glob
import numpy as np
d, pre = sys.argv[1], sys.argv[2]; ref = sys.argv[3] if len(sys.argv) > 3 else None; files = sorted(glob.glob(f"{d}/{pre}*.json")); R = {}
for i, p in enumerate(files):
    for r in json.load(open(p)):
        for k in ("S", "T2"):
            if k in r:
                R.setdefault((k, r.get("state", ""), r["method"]), {}).setdefault(i, {})[r["m"]] = r[k]
print(pre, "files", len(files))
for (k, st, meth), v in R.items():
    a = np.array([[v[i][m] for m in sorted(v[i])] for i in sorted(v)]); line = "%-3s %-5s %-10s " % (k, st, meth) + " ".join("%.3f" % x for x in a.mean(0)) + " | mean %.3f" % a.mean()
    if ref and (k, st, ref) in R and meth != ref:
        b = np.array([[R[(k, st, ref)][i][m] for m in sorted(v[i])] for i in sorted(v)]); line += "  vs %s %+.3f pos %d/%d" % (ref, (a - b).mean(), ((a - b).mean(1) > 0).sum(), len(a))
    print(line)
