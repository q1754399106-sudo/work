"""Distance of each Volve wellbore from the seq-79 line (cable 5081) in the reservoir depth range.
Line coordinates: origin at the Volve F site, inline along 284 deg, crossline along 14 deg; the
reflection points of this line lie at crossline ~ +115 m (mean of source 60-118 m and receiver 130-156 m).
usage: well_line_distance.py <traj_dir>"""
import sys, glob, math, os, re
d = sys.argv[1]; a = math.radians(14); x0, y0 = 435050.02, 6478563.53; CL = 115.0
cl = lambda x, y: (x - x0) * math.sin(a) + (y - y0) * math.cos(a)
il = lambda x, y: -(x - x0) * math.cos(a) + (y - y0) * math.sin(a)
rows = []
for f in sorted(glob.glob(os.path.join(d, '*.txt'))):
    txt = open(f, encoding='latin-1').read(); kb = float(re.search(r'KB-MSL:\s*([\d.]+)', txt).group(1))
    name = re.search(r'WELLBORE NAME:\s*(.+)', txt).group(1).strip(); pts = []
    for ln in txt.splitlines():
        p = ln.split()
        if len(p) >= 8 and re.match(r'^\d', p[0]):
            try: pts.append((float(p[3]) - kb, float(p[6]), float(p[7])))
            except ValueError: pass
    res = [(t, cl(x, y) - CL, il(x, y)) for t, x, y in pts if 2600 <= t <= 3300]
    if not res:
        rows.append((name, None, None, None, max(t for t, _, _ in pts))); continue
    dmin = min(abs(r[1]) for r in res); dmed = sorted(abs(r[1]) for r in res)[len(res) // 2]
    rows.append((name, dmin, dmed, (min(r[2] for r in res), max(r[2] for r in res)), max(t for t, _, _ in pts)))
print(f"{'wellbore':10s} {'min|dCL| m':>11s} {'median':>8s}  inline range (m) in TVDSS 2600-3300   max TVDSS")
for n, a1, a2, r, tmax in sorted(rows, key=lambda r: (r[1] is None, r[1] or 0)):
    if a1 is None: print(f"{n:10s} {'-':>11s} {'-':>8s}  (does not reach 2600 m)  {tmax:8.0f}")
    else: print(f"{n:10s} {a1:11.0f} {a2:8.0f}  {r[0]:7.0f} .. {r[1]:7.0f}   {tmax:8.0f}")
