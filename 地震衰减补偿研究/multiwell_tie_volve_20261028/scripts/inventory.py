"""Well inventory for a multi-well tie (descriptive, before any plan): for every downloaded composite log,
the MD interval where AC and DEN are both valid, its TWT interval and position from the Volve well picks
(MD/TWT/Easting/Northing), and the distance to the nearest MIG_RAW trace. Writes inventory/wells.csv and
inventory/logs/<well>.csv (MD m, AC us/ft, DEN g/cm3).

The MIG_RAW cube is the local copy of experiment 6 (checked against its SHA-256 in experiment 7/8).
"""
import sys, csv, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "field_vsp_volve_f11_20261026" / "cache" / "deps"))   # dlisio 1.0.4
import numpy as np  # noqa: E402
from dlisio import dlis, lis  # noqa: E402

CUBE = ROOT.parent / "field_vsp_volve_f11_20261026" / "cache" / "seismic" / "MIG_RAW.segy"
OUT = ROOT / "inventory"
UNIT = {"mm": 1e-3, "0.1 in": 0.00254, "M": 1.0, "m": 1.0}


def read_dlis(p):
    with dlis.load(str(p)) as fs:
        for f in fs:
            for fr in f.frames:
                names = {c.name: c for c in fr.channels}
                if {"AC", "DEN"} <= set(names):
                    cv = fr.curves(); idx = fr.index
                    return cv[idx] * UNIT[names[idx].units], cv["AC"].astype(float), cv["DEN"].astype(float)
    return None


def read_lis(p):
    with lis.load(str(p)) as fs:
        for f in fs:
            for spec in f.data_format_specs():
                m = [s.mnemonic.strip() for s in spec.specs]
                if {"AC", "DEN"} <= set(m):
                    cv = lis.curves(f, spec); names = cv.dtype.names
                    pick = lambda k: names[m.index(k)]  # noqa: E731
                    return cv[pick("DEPT")] * UNIT[spec.specs[0].units.strip()], cv[pick("AC")].astype(float), cv[pick("DEN")].astype(float)
    return None


def picks():
    """{well name without 'NO ': [(surface, MD, TWT s, E, N), ...]} from the fixed-width pick file."""
    out, spans = {}, None
    for line in (ROOT / "raw" / "picks" / "Well_picks_Volve_v1.dat").read_text(encoding="latin-1").splitlines():
        if line.strip().startswith("-----"):
            spans = [(m.start(), m.end()) for m in re.finditer(r"-+", line)]; continue
        if spans is None or not line.strip() or line.startswith("Well NO") or line.lstrip().startswith("Well name"):
            continue
        cols = [line[a:b + 1].strip() if i < len(spans) - 1 else line[a:].strip() for i, (a, b) in enumerate(spans)]
        try:
            md, twt, e, n = float(cols[4]), float(cols[7]), float(cols[10]), float(cols[11])
        except (ValueError, IndexError):
            continue
        out.setdefault(cols[0].replace("NO ", ""), []).append((cols[1], md, twt / 1000, e, n))
    return out


def cube_xy():
    d = np.memmap(CUBE, dtype=np.uint8, mode="r"); ns = 1126; tb = 240 + 4 * ns; n = (len(d) - 3600) // tb
    h = np.asarray(d[3600:3600 + n * tb]).reshape(n, tb)[:, :240]
    sc = h[:, 70:72].copy().view(">i2").ravel().astype(float); sc = np.where(sc < 0, -1 / sc, np.where(sc == 0, 1, sc))
    x = h[:, 180:184].copy().view(">i4").ravel() * sc; y = h[:, 184:188].copy().view(">i4").ravel() * sc
    il = h[:, 188:192].copy().view(">i4").ravel(); xl = h[:, 192:196].copy().view(">i4").ravel()
    return x, y, il, xl


def main():
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    P = picks(); x, y, il, xl = cube_xy(); rows = []
    for d in sorted((ROOT / "raw" / "composite").iterdir()):
        f = next(d.iterdir()); well = d.name.replace("15_9-", "15/9-").replace("_", " ")
        r = read_dlis(f) if f.suffix.upper() == ".DLIS" else read_lis(f)
        row = dict(well=well, file=f.name, ac_and_den=r is not None)
        if r is not None:
            md, ac, den = r
            for v in (ac, den):
                v[(v == -999.25) | ~np.isfinite(v) | (v <= 0)] = np.nan
            o = np.argsort(md); md, ac, den = md[o], ac[o], den[o]; ok = np.isfinite(ac) & np.isfinite(den)
            np.savetxt(OUT / "logs" / f"{d.name}.csv", np.c_[md, ac, den], delimiter=",", fmt="%.4f",
                       header="MD_m,AC_us_per_ft,DEN_g_cm3", comments="")
            row.update(valid_samples=int(ok.sum()))
            if ok.any():
                row.update(md_top=float(md[ok].min()), md_base=float(md[ok].max()))
                pk = sorted(P.get(well, []), key=lambda t: t[1])
                row["picks_with_twt"] = len(pk)
                if len(pk) >= 2:
                    pm = np.array([t[1] for t in pk]); pt = np.array([t[2] for t in pk])
                    pe = np.array([t[3] for t in pk]); pn = np.array([t[4] for t in pk])
                    lo, hi = max(row["md_top"], pm.min()), min(row["md_base"], pm.max())
                    row.update(md_covered_by_picks=[float(lo), float(hi)] if hi > lo else None)
                    if hi > lo:
                        t0, t1 = np.interp([lo, hi], pm, pt); mid = (lo + hi) / 2
                        e, n = np.interp(mid, pm, pe), np.interp(mid, pm, pn)
                        j = int(np.argmin((x - e) ** 2 + (y - n) ** 2))
                        e1, n1 = np.interp([lo, hi], pm, pe), np.interp([lo, hi], pm, pn)
                        row.update(twt_top=float(t0), twt_base=float(t1), twt_len=float(t1 - t0),
                                   mid_E=float(e), mid_N=float(n), nearest_inline=int(il[j]), nearest_xline=int(xl[j]),
                                   nearest_trace_m=float(np.hypot(x[j] - e, y[j] - n)),
                                   lateral_drift_m=float(np.hypot(e1[1] - e1[0], n1[1] - n1[0])))
        rows.append(row)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
    with open(OUT / "wells.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["well", "file", "ac_and_den", "valid_samples", "md_top", "md_base",
                                           "picks_with_twt", "md_covered_by_picks", "twt_top", "twt_base", "twt_len",
                                           "mid_E", "mid_N", "nearest_inline", "nearest_xline", "nearest_trace_m",
                                           "lateral_drift_m"])
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
