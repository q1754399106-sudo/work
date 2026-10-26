"""Coverage probe: does the ST0202 MIG_RAW cube cover the 15/9-F-11 T2 VSP receivers? (headers only)

Receivers: all NI levels and the RS levels at MD <= 2228.8 m (segment RS of PLAN_ZH.md).

Reads trace-header positions of the cube restored into cache/seismic/MIG_RAW.segy and of the raw Z shots;
writes probe/cube_coverage.json. No amplitudes are read.
"""
import json
import struct
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CUBE = ROOT / 'cache' / 'seismic' / 'MIG_RAW.segy'


def cube_positions(path):
    with open(path, 'rb') as f:
        f.seek(3200)
        bh = f.read(400)
        ns = struct.unpack('>H', bh[20:22])[0]
        tb = 240 + 4 * ns
        n = (path.stat().st_size - 3600) // tb
        rows = np.empty((n, 4))
        for j in range(n):
            f.seek(3600 + j * tb)
            h = f.read(240)
            sc = struct.unpack('>h', h[70:72])[0]
            s = -1 / sc if sc < 0 else (sc or 1)
            rows[j] = (struct.unpack('>i', h[188:192])[0], struct.unpack('>i', h[192:196])[0],
                       struct.unpack('>i', h[180:184])[0] * s, struct.unpack('>i', h[184:188])[0] * s)
    return rows  # inline, xline, x, y


def main():
    geo = json.loads((ROOT / 'probe' / 'vsp_geometry.json').read_text(encoding='utf-8'))
    c = cube_positions(CUBE)
    print(f"cube: {len(c)} traces, inline {c[:, 0].min():.0f}-{c[:, 0].max():.0f}, xline {c[:, 1].min():.0f}-{c[:, 1].max():.0f}, "
          f"x {c[:, 2].min():.1f}-{c[:, 2].max():.1f}, y {c[:, 3].min():.1f}-{c[:, 3].max():.1f}")
    out = []
    levels = [dict(l, survey='RS') for l in geo['RS'] if l['md'] <= 2228.8] + [dict(l, survey='NI') for l in geo['NI']]
    for l in levels:
        d = np.hypot(c[:, 2] - l['rec_x'], c[:, 3] - l['rec_y'])
        j = int(np.argmin(d))
        out.append(dict(survey=l['survey'], trace=l['trace'], md=l['md'], tvdss=l['tvdss'], nearest_inline=int(c[j, 0]),
                        nearest_xline=int(c[j, 1]), distance_m=round(float(d[j]), 1)))
    dist = np.array([o['distance_m'] for o in out])
    print(f"receiver-to-nearest-bin distance: median {np.median(dist):.1f} max {dist.max():.1f} m")
    for o in out[:: max(1, len(out) // 6)] + [out[-1]]:
        print(o)
    (ROOT / 'probe' / 'cube_coverage.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
