"""Tables from results/kirch.json: mean over noise realisations, target zone and whole image."""
import json
from common import ROOT, MS
r = json.load(open(ROOT / "results" / "kirch.json")); names = ("NONE", "PH", "LP", "STAB", "WIEN", "UQ", "UQ25")
print("effective Q %.1f, max path t* %.4f s" % (r["tq"], r["amax"]))
for tag in ("none", "30dB", "14dB"):
    reals = [k for k in r[tag] if k.startswith("real")]
    print("\n== noise %s (%d realisation(s)); tuned LP cutoff %s Hz, STAB gain limit %s dB" % (
        tag, len(reals), r[tag]["tune"]["LP"]["best"], r[tag]["tune"]["STAB"]["best"]))
    for zone in ("target", "all"):
        for met in ("E", "corr"):
            print("  %s, %s:   m = %s | mean | worst" % (zone, met, "  ".join("%5g" % m for m in MS)))
            for n in names:
                v = [sum(r[tag][k]["%s_m%g" % (n, m)][zone][met] for k in reals) / len(reals) for m in MS]
                print("   %-5s        %s | %6.3f | %6.3f" % (n, " ".join("%6.3f" % x for x in v), sum(v) / len(v), min(v)))
