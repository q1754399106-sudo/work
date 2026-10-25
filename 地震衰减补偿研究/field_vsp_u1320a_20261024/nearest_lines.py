"""Find NAMSS 2D lines closest to IODP Hole U1320A from SEG-P1 navigation (NAD27 vs WGS84 shift ignored, ~tens of m)."""
import json, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WELL = (27 + 18.0809 / 60, -(94 + 23.2537 / 60))


def dms(s, neg):
    d = int(s[:-6]); m = int(s[-6:-4]); sec = int(s[-4:]) / 100
    v = d + m / 60 + sec / 3600
    return -v if neg else v


def dist_km(a, b):
    p1, p2 = math.radians(a[0]), math.radians(b[0]); dl = math.radians(b[1] - a[1])
    return 6371 * math.acos(min(1, math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)))


rows = []
for survey in ['b-12-88-gm', 'b-13-88-gm']:
    best = {}
    for line in (ROOT / 'raw' / 'namss' / f'{survey}.segp1').read_text(errors='replace').splitlines():
        if line.startswith('H') or len(line) < 45:
            continue
        name = line[1:20].strip(); sp = line[20:26].strip()
        lat = line[26:35].strip(); lon = line[35:45].strip()
        try:
            pos = (dms(lat[:-1], lat.endswith('S')), dms(lon[:-1], lon.endswith('W')))
        except ValueError:
            continue
        d = dist_km(WELL, pos)
        if name not in best or d < best[name][0]:
            best[name] = (d, sp, pos)
    for name, (d, sp, pos) in best.items():
        rows.append(dict(survey=survey, line=name, distance_km=round(d, 3), nearest_sp=sp, lat=round(pos[0], 5), lon=round(pos[1], 5)))
rows.sort(key=lambda r: r['distance_km'])
(ROOT / 'nearest_lines.json').write_text(json.dumps(rows[:20], indent=2), encoding='utf-8')
for r in rows[:12]:
    print(r)
