"""Data probe before pre-registration: VSP trace order, log ranges, water-column muting. No method is run."""
import numpy as np
from fieldlib import RAW, read_segy, read_log, vsp_levels, WELL_TRACE

V, dt, h = read_segy(RAW / 'ldeo' / 'vsp_segy' / '308-U1320A_stacks-traces.segy')
print('VSP stacks', V.shape, dt, 'trace_seq', h['trace_seq'].tolist())
lv = vsp_levels(); print('levels (seq, mbsf, transit):', [(a, round(b, 1), round(c, 4)) for a, b, c in lv])
for k in range(V.shape[1]):
    x = V[:, k]; i = int(np.argmax(np.abs(x) > .3 * np.abs(x).max()))
    print(f'  stack {k:2d} seq {h["trace_seq"][k]:2d} first strong sample t={i*dt:.3f}s peak={np.abs(x).max():.3g}')
for name, col in [('dsi2', 'VCO'), ('hldsm', 'RHOM')]:
    L = read_log(name); z = L['DEPTH_WMSF']; v = L[col]; ok = np.isfinite(v)
    print(name, col, 'depth', z[ok].min(), z[ok].max(), 'median', np.nanmedian(v), 'n', ok.sum())
X, dts, hs = read_segy(RAW / 'namss' / 'segy' / '4508-32.segy', WELL_TRACE['4508-32'] - 50, 99)
t = np.arange(X.shape[0]) * dts; sf = np.argmax(np.abs(X) > .3 * np.abs(X).max(0), axis=0)
print('seismic window', X.shape, 'SP', hs['shotpoint'][[0, 49, 98]], 'seafloor t range', t[sf].min(), t[sf].max())
wc = X[(t > 1.4) & (t < 1.85)]; deep = X[(t > 12) & (t < 14.5)]
print('water column MAD median', np.median(np.median(np.abs(wc - np.median(wc, 0)), 0)), 'zeros frac', np.mean(wc == 0))
print('deep MAD median', np.median(np.median(np.abs(deep - np.median(deep, 0)), 0)), 'signal rms 2.0-3.5s', np.sqrt(np.mean(X[(t > 2.0) & (t < 3.5)] ** 2)))
