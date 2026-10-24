"""Locate the trace nearest U1320A on each downloaded line (shotpoint in trace-header bytes 17-20)
and report seafloor time and bandwidth of the well-side trace. Read-only on raw/."""
import json, struct
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
NEAREST = {r['line']: r for r in json.loads((ROOT / 'nearest_lines.json').read_text(encoding='utf-8'))}


def ibm2float(u):
    u = u.astype(np.uint32); sign = np.where(u >> 31, -1.0, 1.0)
    exp = ((u >> 24) & 0x7F).astype(np.int32) - 64; frac = (u & 0x00FFFFFF) / float(1 << 24)
    return sign * frac * np.power(16.0, exp)


out = []
for f in sorted((ROOT / 'raw' / 'namss' / 'segy').glob('*.segy')):
    d = f.read_bytes(); ns = struct.unpack('>H', d[3220:3222])[0]; dt = struct.unpack('>H', d[3216:3218])[0] / 1e6
    tb = 240 + ns * 4; ntr = (len(d) - 3600) // tb
    sps = np.array([struct.unpack('>i', d[3600 + i * tb + 16:3600 + i * tb + 20])[0] for i in range(ntr)])
    target = int(NEAREST[f.stem]['nearest_sp']); idx = int(np.argmin(abs(sps - target)))
    traces = []
    for i in range(max(0, idx - 10), min(ntr, idx + 11)):
        raw = np.frombuffer(d[3600 + i * tb + 240:3600 + (i + 1) * tb], dtype='>u4')
        traces.append(ibm2float(raw))
    A = np.array(traces); x = A[len(A) // 2]; t = np.arange(ns) * dt
    env = np.abs(x); first = int(np.argmax(env > 0.2 * env.max()))
    win = A[:, first:first + int(1.0 / dt)]; spec = np.abs(np.fft.rfft(win, axis=1)).mean(0); fr = np.fft.rfftfreq(win.shape[1], dt)
    p = spec ** 2; c = np.cumsum(p) / p.sum()
    out.append(dict(line=f.stem, distance_km=NEAREST[f.stem]['distance_km'], traces=ntr, sample_interval_s=dt, samples=ns,
                    well_trace_index_1based=idx + 1, well_trace_sp=int(sps[idx]), sp_range=[int(sps.min()), int(sps.max())],
                    first_strong_arrival_s=round(float(t[first]), 3),
                    band_5_95pct_Hz=[round(float(fr[np.searchsorted(c, .05)]), 1), round(float(fr[np.searchsorted(c, .95)]), 1)],
                    peak_Hz=round(float(fr[np.argmax(spec)]), 1)))
(ROOT / 'well_traces.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
for o in out:
    print(o)
