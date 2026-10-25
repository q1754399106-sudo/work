"""Shared readers for the U1320A field experiment (read-only on raw/)."""
import struct
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'raw'
SEAFLOOR_MBRF = 1478.0          # logger's seafloor depth, LDEO hole summary
WELL_TRACE = {'4508-32': 2566, '3677-32': 3826, '3685-32': 2236, '4500-32': 2568}   # 1-based, well_traces.json


def ibm2float(u):
    u = u.astype(np.uint32); sign = np.where(u >> 31, -1.0, 1.0)
    exp = ((u >> 24) & 0x7F).astype(np.int32) - 64; frac = (u & 0x00FFFFFF) / float(1 << 24)
    return sign * frac * np.power(16.0, exp)


def read_segy(path, first=None, count=None):
    """Return (data[nt, ntr], dt, headers dict of int arrays). first is 0-based."""
    d = Path(path).read_bytes(); dt = struct.unpack('>H', d[3216:3218])[0] / 1e6; ns = struct.unpack('>H', d[3220:3222])[0]
    fmt = struct.unpack('>H', d[3224:3226])[0]; assert fmt == 1, 'only IBM float handled'
    tb = 240 + ns * 4; ntr = (len(d) - 3600) // tb
    first = 0 if first is None else max(0, first); last = ntr if count is None else min(ntr, first + count)
    X = np.empty((ns, last - first)); sp = np.empty(last - first, int); seq = sp.copy()
    for k, i in enumerate(range(first, last)):
        o = 3600 + i * tb
        seq[k], sp[k] = struct.unpack('>i', d[o:o + 4])[0], struct.unpack('>i', d[o + 16:o + 20])[0]
        X[:, k] = ibm2float(np.frombuffer(d[o + 240:o + tb], dtype='>u4'))
    return X, dt, dict(trace_seq=seq, shotpoint=sp, first_index=first)


def read_log(name):
    """LDEO standard ASCII log: returns dict column -> array (depth m WMSF first)."""
    lines = (RAW / 'ldeo' / 'standard' / f'308-U1320A_{name}.dat').read_text(errors='replace').splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith('DEPTH'))
    cols = lines[i].split('\t'); vals = np.array([[float(x) for x in l.split('\t')[:len(cols)]] for l in lines[i + 2:] if l.strip()])
    vals[vals <= -999] = np.nan
    return {c.strip(): vals[:, k] for k, c in enumerate(cols)}


def vsp_levels():
    """Valid VSP stack levels: (stack trace index 0-based, depth mbsf, one-way transit s), duplicates kept."""
    txt = (RAW / 'ldeo' / 'vsp_ascii' / '308-U1320A_stacks-summary.txt').read_text(errors='replace').splitlines()
    rows = []
    for l in txt:
        p = l.split()
        if len(p) == 6 and p[0].isdigit() and p[1].isdigit():
            rows.append((int(p[0]), float(p[3]) - SEAFLOOR_MBRF, float(p[5])))
    return rows
