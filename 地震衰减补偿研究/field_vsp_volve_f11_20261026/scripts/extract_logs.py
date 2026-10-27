"""Extract AC (us/ft) and DEN (g/cm3) against MD (m) from raw/composite/WLC_COMPOSITE_1.DLIS (PLAN_ZH.md "T2").

Reads with dlisio 1.0.4 installed in the ignored cache/deps. Missing values (absent or the DLIS null -999.25)
become NaN. Writes inputs/logs_md.csv; no reflectivity or tie is computed here.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "cache" / "deps"))
import numpy as np  # noqa: E402
from dlisio import dlis  # noqa: E402

with dlis.load(str(ROOT / "raw" / "composite" / "WLC_COMPOSITE_1.DLIS")) as files:
    (f,) = files
    (fr,) = [fr for fr in f.frames if {"AC", "DEN"} <= {c.name for c in fr.channels}]
    units = {c.name: c.units for c in fr.channels}
    assert (units["DEPTH"], units["AC"], units["DEN"]) == ("mm", "us/ft", "g/cm3"), units
    c = fr.curves()
md = c["DEPTH"] / 1000.0
ac = c["AC"].astype(float)
den = c["DEN"].astype(float)
for x in (ac, den):
    x[(x == -999.25) | ~np.isfinite(x) | (x <= 0)] = np.nan
o = np.argsort(md)
out = ROOT / "inputs"
out.mkdir(exist_ok=True)
np.savetxt(out / "logs_md.csv", np.c_[md[o], ac[o], den[o]], delimiter=",", fmt="%.4f",
           header="MD_m,AC_us_per_ft,DEN_g_cm3", comments="")
ok = np.isfinite(ac) & np.isfinite(den)
print(f"{len(md)} samples, MD {md.min():.1f}-{md.max():.1f} m; AC and DEN both valid: {ok.sum()} "
      f"(MD {md[ok].min():.1f}-{md[ok].max():.1f} m)")
