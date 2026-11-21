"""PLAN_FIELDNOISE_ZH.md: lateral and shot-to-shot noise coherence on Volve seq 79 / cable 5081 (read only).
usage: field_noise.py <segy> <out.json>"""
import sys, json
import numpy as np

fn, out = sys.argv[1:3]
DT = 0.004; I0 = 63; NW = 100                      # 0.252-0.648 s
OFFMIN = 2600.0; NFFT = 128; DMAX = 20
BANDS = {'A_4_25': (4, 25), 'B_25_60': (25, 60)}

raw = np.memmap(fn, dtype=np.uint8, mode='r')
ns = int(np.frombuffer(raw[3220:3222].tobytes(), '>i2')[0]); tsz = 240 + 4 * ns; n = (raw.size - 3600) // tsz
H = np.ndarray((n, tsz), np.uint8, raw[3600:3600 + n * tsz])


def hdr(o, t='>i4'):
    w = 4 if t == '>i4' else 2
    return np.frombuffer(np.ascontiguousarray(H[:, o:o + w]).tobytes(), t).astype(np.int64)


def ibm2f(b):
    b = b.astype(np.uint32); s = np.where(b >> 31, -1.0, 1.0); e = ((b >> 24) & 0x7F).astype(np.int32) - 64
    m = (b & 0x00FFFFFF).astype(np.float64) / float(1 << 24)
    return s * m * np.power(16.0, e)


fldr = hdr(8); off = hdr(36).astype(float); gx = hdr(80); gy = hdr(84)
shots = np.unique(fldr); ns_ = shots.size
# receiver index from the receiver coordinate (fixed OBC cable)
rkey = gx * 10_000_000 + gy
recs, rid = np.unique(rkey, return_inverse=True)
# order receivers along the cable by the first principal direction
gxy = np.array([[k // 10_000_000, k % 10_000_000] for k in recs], float)
u = np.linalg.svd(gxy - gxy.mean(0), full_matrices=False)[2][0]
order = np.argsort((gxy - gxy.mean(0)) @ u); pos = np.empty_like(order); pos[order] = np.arange(order.size)
rj = pos[rid]; nr = recs.size
sj = np.searchsorted(shots, fldr)
print('shots', ns_, 'receivers', nr, 'traces', n)

# read the window, taper, FFT
w = np.hanning(NW)
tr = ibm2f(np.frombuffer(np.ascontiguousarray(H[:, 240 + 4 * I0:240 + 4 * (I0 + NW)]).tobytes(), '>u4').reshape(n, NW))
X = np.fft.rfft((tr - tr.mean(1, keepdims=True)) * w, NFFT, axis=1)
f = np.fft.rfftfreq(NFFT, DT)
G = np.zeros((ns_, nr, f.size), complex); OFF = np.full((ns_, nr), np.nan)
G[sj, rj] = X; OFF[sj, rj] = off
valid = np.abs(OFF) >= OFFMIN
side = np.sign(np.nan_to_num(OFF))

# QC: rms of last vs first 25 samples of the window
q_first = np.sqrt((tr[:, :25] ** 2).mean(1)); q_last = np.sqrt((tr[:, -25:] ** 2).mean(1))
qc = {}
for lo in range(2600, 5001, 400):
    k = (np.abs(off) >= lo) & (np.abs(off) < lo + 400)
    if k.any():
        qc[f'{lo}-{lo+400}'] = dict(n=int(k.sum()), ratio_median=float(np.median(q_last[k] / q_first[k])))

rng = np.random.default_rng(20261121)
EPS = rng.choice([-1.0, 1.0], size=(ns_, nr))[:, :, None]


def lateral(Gx):
    C = []
    for d in range(1, DMAX + 1):
        a, b = Gx[:, :-d], Gx[:, d:]
        m = valid[:, :-d] & valid[:, d:] & (side[:, :-d] == side[:, d:])
        a = a * m[:, :, None]; b = b * m[:, :, None]
        C.append(abs((a * b.conj()).sum((0, 1))) / np.sqrt((abs(a) ** 2).sum((0, 1)) * (abs(b) ** 2).sum((0, 1))))
    return np.array(C)                                   # (DMAX, nf)


def crossshot(Gx, k):
    a, b = Gx[:-k], Gx[k:]
    m = valid[:-k] & valid[k:]
    a = a * m[:, :, None]; b = b * m[:, :, None]
    return abs((a * b.conj()).sum((0, 1))) / np.sqrt((abs(a) ** 2).sum((0, 1)) * (abs(b) ** 2).sum((0, 1)))


C = lateral(G); C0 = lateral(G * EPS)
tmax = 1 + 2 * np.clip(C - C0, 0, None).sum(0)
Gam = {k: crossshot(G, k) for k in (1, 2)}; Gam0 = {k: crossshot(G * EPS, k) for k in (1, 2)}

res = dict(n_shots=int(ns_), n_receivers=int(nr), window_s=[I0 * DT, (I0 + NW - 1) * DT], offmin=OFFMIN,
           n_valid_traces=int(valid.sum()), qc=qc, f=f.tolist(), bands={})
for name, (lo, hi) in BANDS.items():
    k = (f >= lo) & (f <= hi)
    res['bands'][name] = dict(
        tmax_median=float(np.median(tmax[k])),
        C_by_lag_median=[float(np.median(C[d, k])) for d in range(DMAX)],
        Cnull_by_lag_median=[float(np.median(C0[d, k])) for d in range(DMAX)],
        Gamma1_minus_null_median=float(np.median(Gam[1][k] - Gam0[1][k])),
        Gamma2_minus_null_median=float(np.median(Gam[2][k] - Gam0[2][k])),
        Gamma1_median=float(np.median(Gam[1][k])), Gamma1_null_median=float(np.median(Gam0[1][k])))
res['by_f'] = dict(tmax=tmax.tolist(), C1=C[0].tolist(), C1null=C0[0].tolist(),
                   Gamma1=Gam[1].tolist(), Gamma1null=Gam0[1].tolist())
a = res['bands']['A_4_25']
res['L1'] = 'negligible' if a['tmax_median'] <= 1.5 else ('per-trace flip underestimates >=2x' if a['tmax_median'] >= 2 else 'undecided')
g = a['Gamma1_minus_null_median']
res['S1'] = 'holds' if g <= 0.1 else ('fails' if g >= 0.3 else 'undecided')
json.dump(res, open(out, 'w', encoding='utf-8'), indent=1)
print('QC', qc)
for nm, b in res['bands'].items():
    print(nm, 'tmax %.2f' % b['tmax_median'], 'C lag1-5', [round(x, 3) for x in b['C_by_lag_median'][:5]],
          'null', round(b['Cnull_by_lag_median'][0], 3), '| Gamma1 %.3f null %.3f, G1-null %.3f, G2-null %.3f' %
          (b['Gamma1_median'], b['Gamma1_null_median'], b['Gamma1_minus_null_median'], b['Gamma2_minus_null_median']))
print('L1:', res['L1'], ' S1:', res['S1'])
