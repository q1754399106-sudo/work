"""Step-0 geometric check for a focusing (subsurface-offset) noise estimator in UQ-RTM.

Same far-field model and acquisition as direction_theory_20261121/scripts/t_blockflip.py
(constant v, 2-D Green's functions, 300 receivers x 8 m, 51 shots x 40 m, depth 8 m, W_deep points,
10-50 Hz), extended with a horizontal subsurface offset h:

  g_sj(h) = conj( G(x_s -> (x-h, z)) * G(x_rj -> (x+h, z)) )       extended imaging kernel
  signal  I_s(h) = sum_sj g_sj(h) d_sj,  d = flat reflector through the image point
  noise   N(h)   = E|sum_sj g_sj(h) n_sj|^2

Noise classes
  m : independent (receiver and shot)
  l : Gaussian lateral correlation, sigma = 5 receivers (the l30 recipe), independent between shots
  cr: one linear event per shot (apparent velocity c), full lateral coherence, random sign per shot
  c : the same linear event in every shot (shot-coherent, deterministic) -> the case B = 300 cannot see

Estimators of 1/snr = N(0)/S(0), all in the form pf / (p - pf), p = S(0) + N(0):
  B1   : one random sign per trace      pf = sum_sj |g n|^2-type energy + leak L1
  B300 : one random sign per shot       pf = sum_s |sum_j g n|^2-type energy + leak L300
  H    : far subsurface offset          pf = mean_{|h| in {1.5, 2, 2.5} lambda} [ N(h) + |I_s(h)|^2 ]
         (lag rule written before any bias number was computed)
e = log2(est / true): median over points, then median over frequencies (as t_blockflip.py).
"""
import json
import sys
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else 't_hfocus.json'
xr = np.arange(300) * 8.0
xs = 200 + np.arange(51) * 40.0
zs = 8.0
SIG = 5.0
DMAX = 30
pts = [(x, z) for x in np.arange(300, 2101, 150) for z in np.arange(820, 1551, 120)]
freqs = np.arange(10, 51, 4.0)
HFAR = [1.5, 2.0, 2.5]                      # |h| / lambda used by the H estimator (fixed in advance)
HPROF = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]   # profile, descriptive
Cl = np.exp(-np.arange(DMAX + 1) ** 2 / (4 * SIG ** 2))
QS = [0.3, 1.0, 3.0]


def green(xa, x, z, f, v):
    r = np.hypot(xa - x, z - zs)
    return np.exp(-2j * np.pi * f * r / v) / np.sqrt(r)


def quad_l(g):
    """sum_s sum_jk g_sj g_sk* C(j-k), Gaussian lateral correlation."""
    tot = (abs(g) ** 2).sum()
    for d in range(1, DMAX + 1):
        tot += 2 * Cl[d] * (g[:, d:] * g[:, :-d].conj()).sum().real
    return tot


def point(x, z, f, v, c_noise, vmig=None):
    vm = v if vmig is None else vmig
    lam = vm / f
    T = np.hypot(xr[None, :] - xs[:, None], 2 * (z - zs)) / v
    d = np.exp(-2j * np.pi * f * T) / np.sqrt(T * v)                 # flat reflector through (x, z)
    off = abs(xr[None, :] - xs[:, None])
    nc = np.exp(-2j * np.pi * f * off / c_noise) / np.sqrt(np.maximum(off, 25.0))   # linear event
    res = {}
    for hl in sorted(set(HPROF + HFAR)):
        for sgn in ([1] if hl == 0 else [1, -1]):
            h = sgn * hl * lam
            g = np.conj(green(xs, x - h, z, f, vm)[:, None] * green(xr, x + h, z, f, vm)[None, :])
            s = g * d
            Is = s.sum()
            gc = g * nc
            rec = dict(S=abs(Is) ** 2,
                       m=(abs(g) ** 2).sum(),
                       l=quad_l(g),
                       cr=(abs(gc.sum(axis=1)) ** 2).sum(),       # random sign per shot: incoherent over s
                       c=abs(gc.sum()) ** 2)                      # identical in every shot: coherent over s
            if hl == 0:
                rec['L1'] = (abs(s) ** 2).sum()
                rec['L300'] = (abs(s.sum(axis=1)) ** 2).sum()
                rec['B1_c'] = (abs(gc) ** 2).sum()                # per-trace energy of the coherent event
                rec['B300_c'] = (abs(gc.sum(axis=1)) ** 2).sum()  # per-shot energy of the coherent event
                rec['B300_l'] = rec['l']
                rec['B1_l'] = rec['m']
            res[(hl, sgn)] = rec
    return res


def estimates(R, kind, q):
    """Return dict estimator -> log2(est/true) for one point/frequency. Noise scaled so N(0)/S(0) = q."""
    r0 = R[(0.0, 1)]
    S = r0['S']
    k = q * S / r0[kind]                    # scale noise so that true 1/snr = q
    Nt = k * r0[kind]
    p = S + Nt
    out = {}
    # B1: per-trace random signs
    if kind == 'm':
        nb1 = k * r0['m']
    elif kind == 'l':
        nb1 = k * r0['B1_l']
    else:                                   # c, cr: per-trace flipping keeps only the trace energies
        nb1 = k * r0['B1_c']
    # B300: per-shot random signs
    if kind in ('m', 'l', 'cr'):
        nb300 = Nt                          # shot-independent noise: unbiased noise part
    else:
        nb300 = k * r0['B300_c']
    for name, nb, leak in [('B1', nb1, r0['L1']), ('B300', nb300, r0['L300'])]:
        pf = nb + leak
        out[name] = np.log2(max(pf / max(p - pf, 1e-30 * S), 1e-30) / q)
    far = [R[(hl, sg)] for hl in HFAR for sg in (1, -1)]
    pf = np.mean([k * r[kind] + r['S'] for r in far])
    out['H'] = np.log2(max(pf / max(p - pf, 1e-30 * S), 1e-30) / q)
    return out


def run(v, c_noise, vmig=None):
    allR = {f: [point(x, z, f, v, c_noise, vmig) for (x, z) in pts] for f in freqs}
    tab = {}
    for kind in ['m', 'l', 'cr', 'c']:
        for q in QS:
            for est in ['B1', 'B300', 'H']:
                per_f = [np.median([estimates(R, kind, q)[est] for R in allR[f]]) for f in freqs]
                tab[f'{kind}_q{q:g}_{est}'] = float(np.median(per_f))
                tab[f'{kind}_q{q:g}_{est}_byf'] = [float(a) for a in per_f]
    prof = {}
    for kind in ['S', 'm', 'l', 'cr', 'c']:
        for hl in HPROF[1:]:
            vals = []
            for f in freqs:
                vals.append(np.median([0.5 * (R[(hl, 1)][kind] + R[(hl, -1)][kind]) / R[(0.0, 1)][kind]
                                       for R in allR[f]]))
            prof[f'{kind}_{hl:g}'] = float(np.median(vals))
    # shot-coherence factor of the coherent event: true / per-shot sum (what B300 sees)
    coh = float(np.median([np.median([R[(0.0, 1)]['c'] / R[(0.0, 1)]['B300_c'] for R in allR[f]]) for f in freqs]))
    return dict(e=tab, profile=prof, shot_coherence_c=coh)


if __name__ == '__main__':
    res = {}
    for v in [2300.0, 2700.0]:
        for cn in [1500.0, 1000.0]:
            key = f'v{int(v)}_c{int(cn)}'
            res[key] = run(v, cn)
            r = res[key]
            print(f'--- {key}: coherent-event shot coherence (true/B300 view) {r["shot_coherence_c"]:.2f}')
            print('    profile ratio X(h)/X(0), median:')
            for kind in ['S', 'm', 'l', 'cr', 'c']:
                print('     ', kind, ' '.join(f'{hl:g}:{r["profile"][f"{kind}_{hl:g}"]:.3f}' for hl in HPROF[1:]))
            for kind in ['m', 'l', 'cr', 'c']:
                for q in QS:
                    print(f'    {kind:2s} q {q:g}: ' + '  '.join(
                        f'{est} {r["e"][f"{kind}_q{q:g}_{est}"]:+.2f}' for est in ['B1', 'B300', 'H']))
    # velocity error: migrate with 3 % fast / slow velocity, independent noise, q = 1 (descriptive)
    for vm in [0.97, 1.03]:
        key = f'v2500_vmig{vm:g}'
        r = run(2500.0, 1500.0, vmig=2500.0 * vm)
        res[key] = r
        print(f'--- {key}: m q1 H {r["e"]["m_q1_H"]:+.2f}  B300 {r["e"]["m_q1_B300"]:+.2f};'
              f'  S(1.5..2.5 lam)/S(0) {np.mean([r["profile"][f"S_{h:g}"] for h in HFAR]):.3f}')
    json.dump(res, open(OUT, 'w'), indent=1)
