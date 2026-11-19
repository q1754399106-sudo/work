"""Stage 4 well input (PLAN_STAGE4_ZH.md section 4.2): reflectivity in depth (10 m, TVDSS) and the well
position along the line, for F-5 (DLIS) and F-14 (LIS). Output cache/fld/well_<name>.bin/.json.
usage (PYTHONPATH must contain dlisio): fld_well.py <traj_dir> <log_dir> <outdir>"""
import sys, os, re, glob, json, math
import numpy as np
from dlisio import dlis, lis

traj_dir, log_dir, out = sys.argv[1:4]
A = math.radians(14); X0, Y0 = 435050.02, 6478563.53
il = lambda x, y: -(x - X0) * math.cos(A) + (y - Y0) * math.sin(A)
cl = lambda x, y: (x - X0) * math.sin(A) + (y - Y0) * math.cos(A)


def traj(name):
    f = [g for g in glob.glob(os.path.join(traj_dir, '*.txt')) if re.search(r'_%s_ACTUAL' % re.escape(name), g)][0]
    txt = open(f, encoding='latin-1').read(); kb = float(re.search(r'KB-MSL:\s*([\d.]+)', txt).group(1)); r = []
    for ln in txt.splitlines():
        p = ln.split()
        if len(p) >= 8 and re.match(r'^\d', p[0]):
            try: r.append([float(p[0]), float(p[3]) - kb, float(p[6]), float(p[7])])
            except ValueError: pass
    return np.array(r), kb                                                       # MD (RKB), TVDSS, E, N


def curves_dlis(fn):
    with dlis.load(fn) as fs:
        fr = fs[0].frames[0]; c = fr.curves(); d = c['DEPTH'] * 0.00254          # 0.1 in -> m
        return d, c['AC'].astype(float), c['DEN'].astype(float)


def curves_lis(fn):
    with lis.load(fn) as fs:
        f = fs[0]; spec = f.data_format_specs()[0]; a = lis.curves(f, spec)
        names = [n.strip() for n in a.dtype.names]; a = {n: a[o] for n, o in zip(names, a.dtype.names)}
        return a['DEPT'].astype(float), a['AC'].astype(float), a['DEN'].astype(float)


for name, fn, rd in (('F-5', 'F-5_WLC_COMPOSITE_1.DLIS', curves_dlis), ('F-14', 'F-14_WLC_COMPOSITE_1.LIS', curves_lis)):
    md, ac, den = rd(os.path.join(log_dir, fn)); T, kb = traj(name)
    ok = (ac > 30) & (ac < 200) & (den > 1.5) & (den < 3.2) & np.isfinite(ac) & np.isfinite(den)
    md, ac, den = md[ok], ac[ok], den[ok]
    tvd = np.interp(md, T[:, 0], T[:, 1]); e = np.interp(md, T[:, 0], T[:, 2]); n = np.interp(md, T[:, 0], T[:, 3])
    vp = 0.3048 / (ac * 1e-6); imp = vp * den
    z = np.arange(np.ceil(tvd.min() / 10) * 10, np.floor(tvd.max() / 10) * 10 + 1, 10.0)
    zi = np.array([imp[(tvd >= a - 5) & (tvd < a + 5)].mean() if np.any((tvd >= a - 5) & (tvd < a + 5)) else np.nan for a in z])
    good = np.isfinite(zi); z, zi = z[good], zi[good]
    r = np.zeros_like(zi); r[1:] = (zi[1:] - zi[:-1]) / (zi[1:] + zi[:-1])
    s = np.interp(z, tvd, [il(a, b) for a, b in zip(e, n)]); c = np.interp(z, tvd, [cl(a, b) for a, b in zip(e, n)])
    np.c_[z, s, c, r].astype(np.float32).tofile(os.path.join(out, f'well_{name}.bin'))
    json.dump(dict(name=name, n=int(z.size), z=[float(z.min()), float(z.max())], s=[float(s.min()), float(s.max())],
                   cross=[float(c.min()), float(c.max())], vp=[float(vp.min()), float(vp.max())]),
              open(os.path.join(out, f'well_{name}.json'), 'w'), indent=1)
    print(name, 'TVDSS %.0f-%.0f m, line coordinate %.0f..%.0f m, crossline %.0f..%.0f m, %d samples' %
          (z.min(), z.max(), s.min(), s.max(), c.min(), c.max(), z.size))
