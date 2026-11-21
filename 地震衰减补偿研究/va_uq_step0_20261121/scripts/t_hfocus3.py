"""Step-0 check, part 3: the estimator as it would be coded in est_noise2.m / est_block.m.

Per centre frequency fc in {15, 25, 35, 45} Hz the real parts of the single-frequency images are summed
over fc +- 5 Hz (1 Hz steps) and squared (box smoothing omitted: one point at a time). The H estimator
uses ONE lag per band, h = kappa * v / (fc - 5), kappa in {1.5, 2, 2.5}, both signs (fixed lag inside a
band keeps the cross-frequency structure of the noise image identical to h = 0).

  A  signal leakage and noise stationarity (m, l) after band summation
  B  multiple-like event (moveout velocity 0.8 v and 0.9 v), estimator bias e (q = N/S = 1)
  C  pure migration artefact: point 300 m above a flat reflector, full aperture and with a
     receiver gap 1000-1400 m; phi = (p - pf)/p is the share an estimator calls signal (ideal <= 0)
"""
import json
import sys
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else 't_hfocus3.json'
xr_all = np.arange(300) * 8.0
xs = 200 + np.arange(51) * 40.0
zs = 8.0
SIG, DMAX = 5.0, 30
Cl = np.exp(-np.arange(DMAX + 1) ** 2 / (4 * SIG ** 2))
pts = [(x, z) for x in np.arange(300, 2101, 150) for z in np.arange(820, 1551, 120)]
FC = [15.0, 25.0, 35.0, 45.0]
KAP = [1.5, 2.0, 2.5]


def green(xa, x, z, f, v):
    r = np.hypot(xa - x, z - zs)
    return np.exp(-2j * np.pi * f * r / v) / np.sqrt(r)


def kern(x, z, h, f, v, xr):
    return np.conj(green(xs, x - h, z, f, v)[:, None] * green(xr, x + h, z, f, v)[None, :])


def refl(f, v, z_r, xr):
    T = np.hypot(xr[None, :] - xs[:, None], 2 * (z_r - zs)) / v
    return np.exp(-2j * np.pi * f * T) / np.sqrt(T * v)


def quad_l(g):
    tot = (abs(g) ** 2).sum()
    for d in range(1, DMAX + 1):
        tot += 2 * Cl[d] * (g[:, d:] * g[:, :-d].conj()).sum().real
    return tot


def band(fc):
    return np.arange(fc - 5, fc + 5.01, 1.0)


def lagset(fc, v):
    lam = v / (fc - 5)
    return [s * k * lam for k in KAP for s in (1, -1)]


def part_ab(x, z, fc, v):
    """signal through the point, noise m and l, multiple-like events a = 0.8, 0.9."""
    t0 = 2 * (z - zs) / v
    rec = {}
    for h in [0.0] + lagset(fc, v):
        sS = 0.0; nm = 0.0; nl = 0.0; mv = {0.8: 0.0, 0.9: 0.0}
        sshot = 0.0; mvshot = {0.8: 0.0, 0.9: 0.0}; mvtr = {0.8: 0.0, 0.9: 0.0}; str_ = 0.0
        for f in band(fc):
            g = kern(x, z, h, f, v, xr_all)
            s = g * refl(f, v, z, xr_all)
            sS = sS + s.sum().real
            nm += 0.5 * (abs(g) ** 2).sum()
            nl += 0.5 * quad_l(g)
            if h == 0.0:
                sshot = sshot + s.sum(axis=1).real
                str_ = str_ + s.real
            for a in mv:
                dm = refl(f, a * v, zs + a * v * t0 / 2, xr_all)
                gm = g * dm
                mv[a] = mv[a] + gm.sum().real
                if h == 0.0:
                    mvshot[a] = mvshot[a] + gm.sum(axis=1).real
                    mvtr[a] = mvtr[a] + gm.real
        r = dict(S=sS ** 2, m=nm, l=nl, mv8=mv[0.8] ** 2, mv9=mv[0.9] ** 2)
        if h == 0.0:
            r.update(L300=(sshot ** 2).sum(), L1=(str_ ** 2).sum(),
                     B300mv8=(mvshot[0.8] ** 2).sum(), B300mv9=(mvshot[0.9] ** 2).sum(),
                     B1mv8=(mvtr[0.8] ** 2).sum(), B1mv9=(mvtr[0.9] ** 2).sum())
        rec[h] = r
    return rec


def part_c(x, z, fc, v, gap):
    xr = xr_all[(xr_all < 1000) | (xr_all > 1400)] if gap else xr_all
    out = {}
    for h in [0.0] + lagset(fc, v):
        acc = 0.0; ashot = 0.0; atr = 0.0
        for f in band(fc):
            s = kern(x, z, h, f, v, xr) * refl(f, v, z + 300.0, xr)
            acc = acc + s.sum().real
            if h == 0.0:
                ashot = ashot + s.sum(axis=1).real
                atr = atr + s.real
        out[h] = acc ** 2
        if h == 0.0:
            out['B300'] = (ashot ** 2).sum(); out['B1'] = (atr ** 2).sum()
    p = out[0.0]
    far = np.mean([out[h] for h in lagset(fc, v)])
    return dict(p=p, B1=(p - out['B1']) / p, B300=(p - out['B300']) / p, H=(p - far) / p)


def e_of(p, pf, q):
    den = p - pf
    return float(np.log2(pf / den / q)) if den > 0 else 30.0


res = {}
for v in [2300.0, 2700.0]:
    A = {fc: [part_ab(x, z, fc, v) for (x, z) in pts] for fc in FC}
    prof = {}
    for kind in ['S', 'm', 'l', 'mv8', 'mv9']:
        prof[kind] = float(np.median([np.median([np.mean([r[h][kind] for h in lagset(fc, v)]) / r[0.0][kind]
                                                  for r in A[fc]]) for fc in FC]))
    eb = {}
    for kind in ['m', 'l', 'mv8', 'mv9']:
        for q in [0.3, 1.0, 3.0]:
            for est in ['B1', 'B300', 'H']:
                per = []
                for fc in FC:
                    vals = []
                    for r in A[fc]:
                        r0 = r[0.0]; S = r0['S']; k = q * S / r0[kind]; Nt = q * S; p = S + Nt
                        far = [r[h] for h in lagset(fc, v)]
                        if est == 'H':
                            pf = np.mean([k * a[kind] + a['S'] for a in far])
                        elif est == 'B300':
                            pf = (Nt if kind in ('m', 'l') else k * r0['B300' + kind]) + r0['L300']
                        else:
                            nb = {'m': r0['m'], 'l': r0['m']}.get(kind)
                            pf = (k * nb if nb is not None else k * r0['B1' + kind]) + r0['L1']
                        vals.append(e_of(p, pf, q))
                    per.append(np.median(vals))
                eb[f'{kind}_q{q:g}_{est}'] = float(np.median(per))
    C = {}
    for gap in [False, True]:
        G = [part_c(x, z, fc, v, gap) for fc in FC for (x, z) in pts]
        C['gap' if gap else 'full'] = {est: dict(phi_median=float(np.median([g[est] for g in G])),
                                                 share_called_noise=float(np.mean([g[est] <= 0 for g in G])))
                                       for est in ['B1', 'B300', 'H']}
    res[f'v{int(v)}'] = dict(profile_far_over_zero=prof, e=eb, artefact=C)
    print(f'--- v {v:.0f}: far-lag / zero-lag (band-summed):', {k: round(x, 3) for k, x in prof.items()})
    for kind in ['m', 'l', 'mv8', 'mv9']:
        for q in [0.3, 1.0, 3.0]:
            print(f'    {kind:4s} q {q:g}: ' + '  '.join(f'{est} {eb[f"{kind}_q{q:g}_{est}"]:+.2f}'
                                                     for est in ['B1', 'B300', 'H']))
    for k2, d in C.items():
        print(f'    artefact ({k2}): ' + '  '.join(
            f'{est} phi {d[est]["phi_median"]:+.2f} / called noise {d[est]["share_called_noise"]:.2f}'
            for est in ['B1', 'B300', 'H']))
json.dump(res, open(OUT, 'w'), indent=1)
