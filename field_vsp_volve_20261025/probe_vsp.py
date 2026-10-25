"""Geometry probe of the 15/9-F-15A VSP (read-only on raw/); writes probe/vsp_geometry.json.

Per stacked level: MD, TVD (MSL), receiver UTM, first-break time and the field records of the
near-field hydrophone that share the level's records. No spectra, Q or tie is computed here.
"""
import json
import struct
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'raw'
DF_MSL = 54.9  # drill floor elevation above MSL (textual header C06)


def headers(path):
    d = path.read_bytes()
    ns = struct.unpack('>H', d[3220:3222])[0]
    tb = 240 + 4 * ns
    ntr = (len(d) - 3600) // tb
    get = lambda o, f: np.array([struct.unpack(f, d[3600 + j * tb + o:3600 + j * tb + o + struct.calcsize(f)])[0]
                                 for j in range(ntr)])
    sc = get(68, '>h').astype(float)  # bytes 69-70: scalar for bytes 37-68 (divide if negative), varies per trace
    sc = np.where(sc < 0, -1 / sc, np.where(sc == 0, 1, sc))
    return dict(field_rec=get(8, '>i'), channel=get(12, '>i'), md=get(36, '>i') * sc, tvd_rkb=get(40, '>i') * sc,
                src_x=get(72, '>i') / 10, src_y=get(76, '>i') / 10, rec_x=get(80, '>i') / 10, rec_y=get(84, '>i') / 10,
                fb_us=get(228, '>i'), ns=ns, dt=struct.unpack('>H', d[3216:3218])[0] / 1e6)


def main():
    st = headers(RAW / 'vsp_computed' / 'VSPNI_COMPUTED_1.SEGY')
    z = headers(RAW / 'vsp_raw' / 'VSPNI_RAW_2.SEGY')
    nfh = headers(RAW / 'vsp_raw' / 'VSPNI_RAW_1.SEGY')
    tvd_msl = st['tvd_rkb'] - DF_MSL
    one_way = st['fb_us'] / 1e6
    print(f"stack: {len(tvd_msl)} levels, MD {st['md'].min():.1f}-{st['md'].max():.1f} m, "
          f"TVD MSL {tvd_msl.min():.1f}-{tvd_msl.max():.1f} m, first break {one_way.min():.4f}-{one_way.max():.4f} s")
    print(f"receiver x {st['rec_x'].min():.1f}-{st['rec_x'].max():.1f}, y {st['rec_y'].min():.1f}-{st['rec_y'].max():.1f}")
    print(f"raw Z: {len(z['md'])} traces, records {z['field_rec'].min()}-{z['field_rec'].max()}, "
          f"unique records {len(np.unique(z['field_rec']))}, channels {sorted(set(z['channel'].tolist()))}")
    print(f"NFH: {len(nfh['field_rec'])} traces, unique records {len(np.unique(nfh['field_rec']))}, dt {nfh['dt']}, ns {nfh['ns']}")
    # Map each stacked level to raw-Z traces at the same MD (within 0.5 m) and to their NFH records.
    levels = []
    for k in range(len(tvd_msl)):
        hit = np.abs(z['md'] - st['md'][k]) < 0.5
        recs = sorted(set(z['field_rec'][hit].tolist()))
        levels.append(dict(level=k + 1, md=round(float(st['md'][k]), 3), tvd_msl=round(float(tvd_msl[k]), 3),
                           rec_x=float(st['rec_x'][k]), rec_y=float(st['rec_y'][k]),
                           src_x=float(st['src_x'][k]), src_y=float(st['src_y'][k]),
                           first_break_s=float(one_way[k]), raw_z_traces=int(hit.sum()), records=recs,
                           nfh_records=[r for r in recs if r in set(nfh['field_rec'].tolist())]))
    n_raw = [l['raw_z_traces'] for l in levels]
    n_nfh = [len(l['nfh_records']) for l in levels]
    print(f"raw-Z traces per level: min {min(n_raw)} median {int(np.median(n_raw))} max {max(n_raw)}; "
          f"levels with no raw match {sum(n == 0 for n in n_raw)}")
    print(f"NFH records per level: min {min(n_nfh)} median {int(np.median(n_nfh))} max {max(n_nfh)}; "
          f"levels without NFH {sum(n == 0 for n in n_nfh)}")
    off = np.hypot(st['src_x'] - st['rec_x'], st['src_y'] - st['rec_y'])
    print(f"source-receiver horizontal offset at stack level: max {off.max():.1f} m")
    out = ROOT / 'probe'
    out.mkdir(exist_ok=True)
    (out / 'vsp_geometry.json').write_text(json.dumps(dict(df_msl=DF_MSL, levels=levels), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
