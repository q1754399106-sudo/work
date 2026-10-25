"""Data probe for the Volve 15/9-F-15A VSP files (read-only on raw/).

Prints the textual header, binary-header essentials and selected trace-header ranges
of every SEG-Y file. No Q estimate, compensation or tie is computed here.
"""
import struct
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'raw'

# Trace-header fields (1-based byte, format) printed as min / max over the file.
FIELDS = {
    'seq_line': (1, '>i'), 'field_rec': (9, '>i'), 'channel': (13, '>i'), 'cdp': (21, '>i'),
    'offset': (37, '>i'), 'rec_elev': (41, '>i'), 'src_depth': (49, '>i'), 'scalar_elev': (69, '>h'),
    'scalar_xy': (71, '>h'), 'src_x': (73, '>i'), 'src_y': (77, '>i'), 'rec_x': (81, '>i'), 'rec_y': (85, '>i'),
    'uphole_t': (95, '>h'), 'delay_ms': (109, '>h'), 'ns': (115, '>H'), 'dt_us': (117, '>H'),
    'cdp_x': (181, '>i'), 'cdp_y': (185, '>i'), 'inline': (189, '>i'), 'xline': (193, '>i'),
    'b197': (197, '>i'), 'b201': (201, '>i'), 'b205': (205, '>i'), 'b209': (209, '>i'), 'b213': (213, '>i'),
    'b217': (217, '>i'), 'b221': (221, '>i'), 'b225': (225, '>i'), 'b229': (229, '>i'), 'b233': (233, '>i'),
    'b237': (237, '>i'),
}


def ebcdic(b):
    try:
        return b.decode('cp500')
    except Exception:
        return b.decode('latin-1')


def probe(path):
    d = path.read_bytes()
    txt = ebcdic(d[:3200])
    if sum(c.isprintable() for c in txt) < 0.8 * len(txt):
        txt = d[:3200].decode('ascii', errors='replace')
    dt = struct.unpack('>H', d[3216:3218])[0]
    ns = struct.unpack('>H', d[3220:3222])[0]
    fmt = struct.unpack('>H', d[3224:3226])[0]
    size = {1: 4, 2: 4, 3: 2, 5: 4, 8: 1}[fmt]
    tb = 240 + ns * size
    ntr = (len(d) - 3600) / tb
    print(f'\n######## {path.relative_to(ROOT)}  bytes={len(d)}  dt_us={dt} ns={ns} fmt={fmt} ntr={ntr}')
    for i in range(0, 3200, 80):
        line = txt[i:i + 80].rstrip()
        if line.strip():
            print('  |', line)
    ntr = int(ntr)
    vals = {k: np.array([struct.unpack(f, d[3600 + j * tb + b - 1:3600 + j * tb + b - 1 + struct.calcsize(f)])[0]
                         for j in range(ntr)]) for k, (b, f) in FIELDS.items()}
    for k, v in vals.items():
        if v.min() != 0 or v.max() != 0:
            print(f'  {k:12s} min={v.min()} max={v.max()} first={v[:4].tolist()} last={v[-2:].tolist()}')


if __name__ == '__main__':
    files = sorted(RAW.rglob('*.SEGY'), key=lambda p: (p.parent.name, int(p.stem.split('_')[-1])))
    for p in files:
        if len(sys.argv) > 1 and not any(a in p.name for a in sys.argv[1:]):
            continue
        probe(p)
