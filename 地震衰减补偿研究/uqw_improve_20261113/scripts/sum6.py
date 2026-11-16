"""Summary of dev6.py: means over the development sets by Q multiplier, differences to VDB and UQ-W."""
import json, glob, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; K = ("UQW", "VDB", "H", "HA", "HAW"); fs = sorted(glob.glob(str(ROOT / "results" / "dev6" / "*.json"))); print("sets", len(fs))
for st in ("cal", "nd"):
    for var in "OF":
        for kk in ("T2", "S"):
            v = {k: np.array([[q[kk] for q in json.load(open(f)) if q["state"] == st and q["var"] == var and q["method"] == k] for f in fs]) for k in K}
            for k in K:
                print(st, var, kk, "%-4s" % k, " ".join("%.3f" % x for x in v[k].mean(0)), "| mean %.3f" % v[k].mean())
            for a in ("HA", "HAW"):
                for b in ("VDB", "UQW"):
                    d = v[a] - v[b]; print("   %s-%s" % (a, b), " ".join("%+.3f" % x for x in d.mean(0)), "| mean %+.3f pos %d/%d" % (d.mean(), (d.mean(1) > 0).sum(), len(fs)))
