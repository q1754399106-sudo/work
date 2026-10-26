"""Geometry probe of the 15/9-F-11 T2 VSPs (NI and RS; read-only on raw/); writes probe/vsp_geometry.json.

Per stacked level of each survey: MD, TVDSS (header bytes 53-56), receiver and source UTM, straight-ray
incidence angle, contractor first-break pick (bytes 233-236, ms), and the field records of the downhole Z
shots and of the surface reference hydrophone that share the level. No spectra, Q or tie is computed here.
"""
import json
import struct
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'raw' / 'vsp_velocity'
KB_MSL = 54.9  # KB elevation above MSL (textual header C07)


def headers(path):
    d = path.read_bytes()
    ns = struct.unpack('>H', d[3220:3222])[0]
    tb = 240 + 4 * ns
    ntr = (len(d) - 3600) // tb
    get = lambda o, f: np.array([struct.unpack(f, d[3600 + j * tb + o:3600 + j * tb + o + struct.calcsize(f)])[0]
                                 for j in range(ntr)])
    sc = get(70, '>h').astype(float)  # bytes 71-72: coordinate scalar
    sc = np.where(sc < 0, -1 / sc, np.where(sc == 0, 1, sc))
    zs = get(68, '>h').astype(float)  # bytes 69-70: elevation/depth scalar
    zs = np.where(zs < 0, -1 / zs, np.where(zs == 0, 1, zs))
    return dict(field_rec=get(8, '>i'), md=get(40, '>i') * zs, tvdss=get(52, '>i') * zs, offset=get(36, '>i') * zs,
                src_x=get(72, '>i') * sc, src_y=get(76, '>i') * sc, rec_x=get(80, '>i') * sc, rec_y=get(84, '>i') * sc,
                fb_ms=get(232, '>i'), ns=ns, dt=struct.unpack('>H', d[3216:3218])[0] / 1e6, ntr=ntr)


def survey(tag):
    st = headers(RAW / f'VSP{tag}_COMPUTED_1.SEGY')
    z = headers(RAW / f'VSP{tag}_RAW_2.SEGY')
    nfh = headers(RAW / f'VSP{tag}_RAW_1.SEGY')
    nfh_recs = set(nfh['field_rec'].tolist())
    levels = []
    for k in range(st['ntr']):
        hit = np.abs(z['md'] - st['md'][k]) < 0.5
        recs = sorted(set(z['field_rec'][hit].tolist()))
        off = float(np.hypot(st['src_x'][k] - st['rec_x'][k], st['src_y'][k] - st['rec_y'][k]))
        levels.append(dict(trace=k + 1, md=round(float(st['md'][k]), 3), tvdss=round(float(st['tvdss'][k]), 3),
                           rec_x=float(st['rec_x'][k]), rec_y=float(st['rec_y'][k]),
                           src_x=float(st['src_x'][k]), src_y=float(st['src_y'][k]), src_rec_horizontal_m=off,
                           incidence_deg=float(np.degrees(np.arctan2(off, st['tvdss'][k]))),
                           fb_ms=int(st['fb_ms'][k]), raw_z_traces=int(hit.sum()), records=recs,
                           nfh_records=[r for r in recs if r in nfh_recs]))
    md = np.array([l['md'] for l in levels]); tvd = np.array([l['tvdss'] for l in levels])
    fb = np.array([l['fb_ms'] for l in levels]); ang = np.array([l['incidence_deg'] for l in levels])
    n_raw = [l['raw_z_traces'] for l in levels]; n_nfh = [len(l['nfh_records']) for l in levels]
    print(f"{tag}: {st['ntr']} levels, MD {md.min():.1f}-{md.max():.1f}, TVDSS {tvd.min():.1f}-{tvd.max():.1f}, "
          f"first break {fb.min()}-{fb.max()} ms; raw Z {z['ntr']} traces, NFH {nfh['ntr']} (dt {nfh['dt']})")
    print(f"  raw-Z traces per level min {min(n_raw)} median {int(np.median(n_raw))}; NFH records per level "
          f"min {min(n_nfh)} median {int(np.median(n_nfh))}; levels without NFH {sum(n == 0 for n in n_nfh)}")
    for cut in (2228.8, 4530.1):
        m = md <= cut
        print(f"  MD <= {cut}: {m.sum()} levels, incidence max {ang[m].max():.1f} deg, "
              f"one-way span {(fb[m].max() - fb[m].min()) / 1000:.3f} s")
    return levels


def main():
    out = ROOT / 'probe'
    out.mkdir(exist_ok=True)
    geo = dict(kb_msl=KB_MSL, NI=survey('NI'), RS=survey('RS'))
    (out / 'vsp_geometry.json').write_text(json.dumps(geo, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
