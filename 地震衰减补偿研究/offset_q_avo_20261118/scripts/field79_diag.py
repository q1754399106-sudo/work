"""Descriptive diagnostic after the field run of PLAN_FIELD_I (not a verdict): is the delivered PZ data spectrally whitened?
1. near-offset (|x| <= 350 m) spectra per window: centroid 8-60 Hz and log-power slope; an attenuating earth without
   whitening makes the centroid fall with time (synthetic of the self-check for comparison).
2. mean 0.5 ln(P_x / P_ref) per offset class for the 2.0 s window (field vs synthetic sign of the slope).
-> results/field79_diag.json"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import field79 as F  # noqa: E402

H, h = F.read_headers(); sel = np.abs(h["cmp"]) < F.CMP_LIM; zw = float(np.median(h["wd"][sel])); R = F.Rays(1.0, zw)
tabs = {tc: R.table(tc) for tc in F.CENT}; k = np.where(sel)[0]; G = F.traces(H, k); off = h["off"][k]
ib = (F.F >= F.BAND[0]) & (F.F <= F.BAND[1]); out = dict(near={}, ratio_2s={})
kr = off <= F.REF; xr = float(np.sqrt(np.mean(off[kr] ** 2)))
for tc in F.CENT:
    tr = float(tabs[tc](np.array([xr]))[0][0]); P = F.power(G[kr], tr - F.WL / 2, tr + F.WL / 2)
    if P is None:
        continue
    c = float((F.F[ib] * P[ib]).sum() / P[ib].sum()); s = float(np.polyfit(F.F[ib], 0.5 * np.log(P[ib]), 1)[0])
    out["near"]["%.1f" % tc] = dict(centroid_hz=c, half_log_power_slope_per_hz=s, power_60_over_peak=float(P[ib][-1] / P[ib].max()))
    print("near tc %.1f centroid %.1f Hz  slope %+.4f /Hz  P(60)/peak %.3f" % (tc, c, s, P[ib][-1] / P[ib].max()), flush=True)
tc = 2.0; tr = float(tabs[tc](np.array([xr]))[0][0]); P0 = F.power(G[kr], tr - F.WL / 2, tr + F.WL / 2)
for o in F.OFFS:
    kk = (off >= o - F.HW) & (off < o + F.HW); xo = float(np.sqrt(np.mean(off[kk] ** 2))); tx = float(tabs[tc](np.array([xo]))[0][0])
    if not np.isfinite(tx) or tx - F.WL / 2 < F.MUTE(o + F.HW) + F.MUTE_PAD:
        continue
    Px = F.power(G[kk], tx - F.WL / 2, tx + F.WL / 2)
    if Px is None:
        continue
    lr = 0.5 * np.log(Px[ib] / P0[ib]); sl = float(np.polyfit(F.F[ib], lr, 1)[0])
    out["ratio_2s"]["%d" % o] = dict(slope_per_hz=sl, at_10_30_50=[float(np.interp(f, F.F[ib], lr - lr[0])) for f in (10, 30, 50)])
    print("2.0 s offset %4d  slope %+.4f /Hz  (0.5 ln ratio rel. 8 Hz at 10/30/50 Hz: %s)" % (o, sl, np.round(out["ratio_2s"]["%d" % o]["at_10_30_50"], 2)), flush=True)
(F.ROOT / "results" / "field79_diag.json").write_text(json.dumps(out, indent=1))
