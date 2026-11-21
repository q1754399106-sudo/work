"""Step-0 check, part 2: kinematically inconsistent but shot-coherent energy (the class that a per-shot
random sign counts as signal, DERIVATION_ZH.md 4.4), and acquisition-gap migration artefacts.

Same model as t_hfocus.py. Classes
  mv : a multiple-like event in every shot: same zero-offset time as the primary at the image point,
       moveout velocity v2 = a * v (a = 0.9, 0.8); deterministic, identical in every shot
  gap: no signal at the image point; a flat reflector 200 m deeper; receivers 1000-1400 m missing
       (platform-gap analogue) -> the image at the point is pure migration artefact
For mv, e = log2(est/true) with the event scaled so that true N(0)/S(0) = q.
For gap, report phi = (p - pf) / p, the fraction of the artefact an estimator calls signal (ideal <= 0).
"""
import json
import sys
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else 't_hfocus2.json'
xr_all = np.arange(300) * 8.0
xs = 200 + np.arange(51) * 40.0
zs = 8.0
pts = [(x, z) for x in np.arange(300, 2101, 150) for z in np.arange(820, 1551, 120)]
freqs = np.arange(10, 51, 4.0)
HFAR = [1.5, 2.0, 2.5]
QS = [0.3, 1.0, 3.0]


def green(xa, x, z, f, v):
    r = np.hypot(xa - x, z - zs)
    return np.exp(-2j * np.pi * f * r / v) / np.sqrt(r)


def kern(x, z, h, f, v, xr):
    return np.conj(green(xs, x - h, z, f, v)[:, None] * green(xr, x + h, z, f, v)[None, :])


def refl(f, v, z_r, xr):
    T = np.hypot(xr[None, :] - xs[:, None], 2 * (z_r - zs)) / v
    return np.exp(-2j * np.pi * f * T) / np.sqrt(T * v)


def lags(f, v):
    lam = v / f
    return [s * a * lam for a in HFAR for s in (1, -1)]


def mv_point(x, z, f, v, a):
    xr = xr_all
    d = refl(f, v, z, xr)                               # primary through the point
    t0 = 2 * (z - zs) / v
    v2 = a * v
    z2 = zs + v2 * t0 / 2                               # same zero-offset time, slower moveout
    dm = refl(f, v2, z2, xr)
    g0 = kern(x, z, 0.0, f, v, xr)
    S = abs((g0 * d).sum()) ** 2
    N = abs((g0 * dm).sum()) ** 2
    B1n = (abs(g0 * dm) ** 2).sum()
    B300n = (abs((g0 * dm).sum(axis=1)) ** 2).sum()
    L1 = (abs(g0 * d) ** 2).sum()
    L300 = (abs((g0 * d).sum(axis=1)) ** 2).sum()
    far_s, far_n = [], []
    for h in lags(f, v):
        g = kern(x, z, h, f, v, xr)
        far_s.append(abs((g * d).sum()) ** 2)
        far_n.append(abs((g * dm).sum()) ** 2)
    return dict(S=S, N=N, B1n=B1n, B300n=B300n, L1=L1, L300=L300, Hs=np.mean(far_s), Hn=np.mean(far_n))


def gap_point(x, z, f, v):
    keep = (xr_all < 1000) | (xr_all > 1400)
    xr = xr_all[keep]
    d = refl(f, v, z + 200.0, xr)                       # reflector below the point
    g0 = kern(x, z, 0.0, f, v, xr)
    s = g0 * d
    p = abs(s.sum()) ** 2
    pf1 = (abs(s) ** 2).sum()
    pf300 = (abs(s.sum(axis=1)) ** 2).sum()
    pfh = np.mean([abs((kern(x, z, h, f, v, xr) * d).sum()) ** 2 for h in lags(f, v)])
    # reference: the same point with full aperture (artefact without the gap)
    d_full = refl(f, v, z + 200.0, xr_all)
    p_full = abs((kern(x, z, 0.0, f, v, xr_all) * d_full).sum()) ** 2
    return dict(p=p, B1=(p - pf1) / p, B300=(p - pf300) / p, H=(p - pfh) / p, gain=p / max(p_full, 1e-300))


def e_mv(r, q):
    k = q * r['S'] / r['N']
    Nt = k * r['N']
    p = r['S'] + Nt
    out = {}
    for name, pf in [('B1', k * r['B1n'] + r['L1']), ('B300', k * r['B300n'] + r['L300']),
                     ('H', k * r['Hn'] + r['Hs'])]:
        den = p - pf
        out[name] = float(np.log2(pf / den / q)) if den > 0 else 30.0   # 30 = 'all noise' (saturated)
    return out


res = {}
for v in [2300.0, 2700.0]:
    for a in [0.9, 0.8]:
        R = {f: [mv_point(x, z, f, v, a) for (x, z) in pts] for f in freqs}
        key = f'mv_v{int(v)}_a{a:g}'
        tab = {}
        for q in QS:
            for est in ['B1', 'B300', 'H']:
                per_f = [np.median([e_mv(r, q)[est] for r in R[f]]) for f in freqs]
                tab[f'q{q:g}_{est}'] = float(np.median(per_f))
                tab[f'q{q:g}_{est}_byf'] = [float(t) for t in per_f]
        coh = float(np.median([np.median([r['N'] / r['B300n'] for r in R[f]]) for f in freqs]))
        hfoc = float(np.median([np.median([r['Hn'] / r['N'] for r in R[f]]) for f in freqs]))
        res[key] = dict(e=tab, shot_coherence=coh, far_over_zero=hfoc)
        print(f'--- {key}: event true/per-shot-sum {coh:.2f}; event far-lag/zero-lag {hfoc:.2f}')
        for q in QS:
            print(f'    q {q:g}: ' + '  '.join(f'{est} {tab[f"q{q:g}_{est}"]:+.2f}' for est in ['B1', 'B300', 'H']))
    G = {f: [gap_point(x, z, f, v) for (x, z) in pts] for f in freqs}
    key = f'gap_v{int(v)}'
    tab = {est: float(np.median([np.median([g[est] for g in G[f]]) for f in freqs])) for est in ['B1', 'B300', 'H']}
    tab['frac_called_noise'] = {est: float(np.mean([g[est] <= 0 for f in freqs for g in G[f]]))
                                for est in ['B1', 'B300', 'H']}
    tab['gap_gain_median'] = float(np.median([g['gain'] for f in freqs for g in G[f]]))
    res[key] = tab
    print(f'--- {key}: artefact fraction called signal (median phi; ideal <= 0): '
          + '  '.join(f'{est} {tab[est]:+.2f}' for est in ['B1', 'B300', 'H'])
          + '; share of points called all-noise ' + str({k: round(x, 2) for k, x in tab['frac_called_noise'].items()})
          + f'; artefact energy vs full aperture x{tab["gap_gain_median"]:.1f}')
json.dump(res, open(OUT, 'w'), indent=1)
