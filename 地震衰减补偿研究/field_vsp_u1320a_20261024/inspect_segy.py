"""Read-only SEG-Y header inspection (textual, binary and selected trace headers).

Usage: py -3 -I inspect_segy.py <file.segy> [...]
Downloaded files are treated as data; nothing in them is executed.
"""
import json, struct, sys
from pathlib import Path

FORMATS = {1: ('IBM float', 4), 2: ('int32', 4), 3: ('int16', 2), 5: ('IEEE float', 4), 8: ('int8', 1)}


def ebcdic(b):
    try:
        return b.decode('cp037')
    except Exception:
        return b.decode('latin-1', 'replace')


def inspect(path):
    data = Path(path).read_bytes()
    text = ebcdic(data[:3200])
    if sum(c.isprintable() for c in text) < 2000:
        text = data[:3200].decode('ascii', 'replace')
    binh = data[3200:3600]
    dt_us, ns, fmt = struct.unpack('>H', binh[16:18])[0], struct.unpack('>H', binh[20:22])[0], struct.unpack('>H', binh[24:26])[0]
    name, size = FORMATS.get(fmt, ('unknown', 4))
    tr_bytes = 240 + ns * size
    ntr = (len(data) - 3600) // tr_bytes
    def th(i):
        h = data[3600 + i * tr_bytes: 3600 + i * tr_bytes + 240]
        g = lambda o, f='>i', n=4: struct.unpack(f, h[o:o + n])[0]
        return dict(trace_seq=g(0), field_record=g(8), cdp=g(20), offset=g(36), scalar=struct.unpack('>h', h[70:72])[0],
                    src_x=g(72), src_y=g(76), cdp_x=g(180), cdp_y=g(184), shotpoint=g(196),
                    ns=struct.unpack('>H', h[114:116])[0], dt_us=struct.unpack('>H', h[116:118])[0])
    return dict(file=Path(path).name, bytes=len(data), sample_interval_us=dt_us, samples=ns, format=fmt, format_name=name,
                traces=ntr, remainder_bytes=(len(data) - 3600) % tr_bytes, record_length_s=ns * dt_us / 1e6,
                first_trace=th(0), last_trace=th(ntr - 1) if ntr else None,
                text_header=[text[i:i + 80].rstrip() for i in range(0, 3200, 80)])


if __name__ == '__main__':
    out = [inspect(p) for p in sys.argv[1:]]
    for o in out:
        print(json.dumps({k: v for k, v in o.items() if k != 'text_header'}, indent=1))
        print('\n'.join(l for l in o['text_header'] if l.strip())[:2500])
