"""Geometry probe of the ST0202R08 raw PSDM time stack (read-only on raw/); writes probe/cube_geometry.json.

Reads binary and trace headers only (inline, crossline, bin coordinates), fits a straight line to
the VSP receiver positions of levels 1-191 and reports the nearest bins along it. No trace
amplitudes are used.
"""
import json
import struct
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CUBE = ROOT / 'raw' / 'seismic' / 'ST0202R08_PZ_PSDM_RAW_STACK_PP_TIME.MIG_RAW.POST_STACK.3D.JS-017534.segy'


def main():
    with open(CUBE, 'rb') as fh:
        head = fh.read(3600)
        dt = struct.unpack('>H', head[3216:3218])[0]
        ns = struct.unpack('>H', head[3220:3222])[0]
        fmt = struct.unpack('>H', head[3224:3226])[0]
        tb = 240 + 4 * ns
        ntr = (CUBE.stat().st_size - 3600) // tb
        th = np.empty((ntr, 240), np.uint8)
        for j in range(ntr):
            fh.seek(3600 + j * tb)
            th[j] = np.frombuffer(fh.read(240), np.uint8)
    be = lambda o, t: th[:, o:o + np.dtype(t).itemsize].copy().view(np.dtype(t).newbyteorder('>')).ravel()
    il, xl = be(188, 'i4'), be(192, 'i4')
    sc = be(70, 'i2').astype(float)
    sc = np.where(sc < 0, -1 / sc, np.where(sc == 0, 1, sc))
    x, y = be(180, 'i4') * sc, be(184, 'i4') * sc
    print(f'dt {dt} us, ns {ns}, format {fmt}, traces {ntr}')
    print(f'inline {il.min()}-{il.max()} (step {np.min(np.diff(np.unique(il)))}), '
          f'xline {xl.min()}-{xl.max()} (step {np.min(np.diff(np.unique(xl)))}), xy scalar {np.unique(sc)}')
    lv = json.loads((ROOT / 'probe' / 'vsp_geometry.json').read_text(encoding='utf-8'))['levels'][:191]
    rx = np.array([l['rec_x'] for l in lv]); ry = np.array([l['rec_y'] for l in lv])
    fb = np.array([l['first_break_s'] for l in lv])
    c = np.array([rx.mean(), ry.mean()])
    u = np.linalg.svd(np.c_[rx - c[0], ry - c[1]], full_matrices=False)[2][0]
    if u[0] > 0:
        u = -u  # point from the shallow (east) end toward the bottom hole
    s = (rx - c[0]) * u[0] + (ry - c[1]) * u[1]
    perp = -(rx - c[0]) * u[1] + (ry - c[1]) * u[0]
    print(f'well-path line: centre {c.round(1)}, azimuth {np.degrees(np.arctan2(u[0], u[1])) % 360:.1f} deg, '
          f'receiver along-line {s.min():.0f}..{s.max():.0f} m, max perpendicular {np.abs(perp).max():.1f} m')
    near = []
    for k in [0, 60, 120, 190]:
        d = np.hypot(x - rx[k], y - ry[k]); j = int(np.argmin(d))
        near.append(dict(level=k + 1, twt=round(2 * fb[k], 4), nearest_inline=int(il[j]), nearest_xline=int(xl[j]),
                         distance_m=round(float(d[j]), 1)))
        print(near[-1])
    (ROOT / 'probe').mkdir(exist_ok=True)
    (ROOT / 'probe' / 'cube_geometry.json').write_text(json.dumps(dict(
        dt_us=int(dt), ns=int(ns), format=int(fmt), traces=int(ntr), inline=[int(il.min()), int(il.max())],
        xline=[int(xl.min()), int(xl.max())], line_centre=c.tolist(), line_unit=u.tolist(),
        receiver_along=[float(s.min()), float(s.max())], nearest=near), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
