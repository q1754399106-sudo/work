"""Theory prediction for uq_rtm_20261117 PLAN_BLOCKFLIP_ZH.md step 1, written before any step-1 number
was seen (blockflip runs on D2/D4 were at 47/105 frequencies when this was written).

Model (DERIVATION_ZH.md section 4): constant velocity v, 2-D far-field Green's functions, geometry of
rtm_setup.m (300 receivers at 8 m, 51 shots at 40 m, depth 8 m), cross-correlation imaging at a point x.
  noise image of shot s:   sum_j g_sj n_sj,  g_sj = conj(G(x_s,x) G(x_rj,x))
  truth  N_t = sum_s sum_jk g_sj g_sk* C(j-k)
  block flip (B receivers share one sign): N_B = sum_s sum_{j,k same block} g_sj g_sk* C(j-k)
  signal (flat reflector through x): s_sj = g_sj d_sj, stacked S = |sum_s sum_j s_sj|^2,
  leak  L_B = sum_s sum_blocks |sum_{j in block} s_sj|^2
  l30: C(D) = exp(-D^2/(4 sigma^2)), sigma = 5 receivers (Gaussian kernel convolved with itself);
  m30/m40: C = delta.
Estimator of est_noise2/est_block: 1/snr_est = pf/(p - pf), p = S + N_t, pf = N_B + L_B; truth N_t/S.
q = N_1/S (independent-noise level per receiver) is a free parameter: the stage-3 m30 estimate at 40 Hz
was 0.86-1.23, so q ~ 1 for m30 and ~10 for m40 (10 dB more noise).
e = log2(est/true) per point and frequency, median over points, then median over 10-50 Hz.
"""
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
xr = np.arange(300) * 8.0
xs = 200 + np.arange(51) * 40.0
zs = 8.0
sig = 5.0
BS = [1, 8, 32, 300]
pts = [(x, z) for x in np.arange(300, 2101, 150) for z in np.arange(820, 1551, 120)]
freqs = np.arange(10, 51, 4.0)
D = np.arange(-300, 301)
Cl = np.exp(-D ** 2 / (4 * sig ** 2))


def green(xa, x, z, f, v):
    r = np.hypot(xa - x, z - zs)
    return np.exp(-2j * np.pi * f * r / v) / np.sqrt(r), r


def blocks(B, n=300):
    return np.arange(n) // B


DMAX = 30   # C(30) = exp(-9) for sigma = 5: truncation error < 1e-3


def quad_corr_all(g, lab):
    """sum over shots of sum_jk g_j g_k* C(j-k) (|j-k| <= DMAX), restricted to lab_j == lab_k if lab given.
    g: (ns, nr)."""
    tot = (abs(g) ** 2).sum()                     # D = 0, C = 1
    for d in range(1, DMAX + 1):
        prod = g[:, d:] * g[:, :-d].conj()        # pairs (j, j-d)
        if lab is not None:
            prod = prod * (lab[d:] == lab[:-d])[None, :]
        tot += 2 * Cl[d + 300] * prod.sum().real
    return tot


def run(v):
    out = {}
    labs = {B: (None if B == 300 else blocks(B)) for B in BS}
    for f in freqs:
        rows = []
        for (x, z) in pts:
            Gr, rr = green(xr, x, z, f, v)
            Gs, rs = green(xs, x, z, f, v)
            # flat reflector at depth z: two-way data, image phase = tau_s + tau_r - T
            T = np.hypot(xr[None, :] - xs[:, None], 2 * (z - zs)) / v
            d = np.exp(-2j * np.pi * f * T) / np.sqrt(T * v)
            g = np.conj(Gs[:, None] * Gr[None, :])          # (ns, nr)
            s = g * d
            S = abs(s.sum()) ** 2
            N1 = (abs(g) ** 2).sum()
            Nt_l = quad_corr_all(g, None)
            rec = dict(S=S, N1=N1, Nt_l=Nt_l)
            for B in BS:
                lab = labs[B]
                if B == 300:
                    NB = Nt_l
                    LB = (abs(s.sum(axis=1)) ** 2).sum()
                else:
                    NB = quad_corr_all(g, lab)
                    nb = lab.max() + 1
                    LB = sum((abs(np.bincount(lab, weights=s[k].real, minlength=nb) +
                                  1j * np.bincount(lab, weights=s[k].imag, minlength=nb)) ** 2).sum()
                             for k in range(len(xs)))
                rec[f'NB{B}'] = NB; rec[f'LB{B}'] = LB
            rows.append(rec)
        out[float(f)] = rows
    return out


def e_of(rows, B, q, kind):
    es = []
    for r in rows:
        S = r['S']; n1 = q * S                      # q = N1/S defines the noise scale
        k = n1 / r['N1']
        lam = r[f'LB{B}'] / S                       # leak relative to stacked signal
        if kind == 'l30':
            Nt = k * r['Nt_l']; NB = k * r[f'NB{B}']
        else:
            Nt = n1; NB = n1
        pf = NB + lam * S
        est = pf / max(S + Nt - pf, 1e-30 * S)
        es.append(np.log2(est / (Nt / S)))
    return float(np.median(es))


res = {}
for v in [2300.0, 2700.0]:
    R = run(v)
    tab = {}
    for kind, qs in [('m30', [0.3, 1.0, 3.0]), ('m40', [3.0, 10.0, 30.0]), ('l30', [0.3, 1.0, 3.0])]:
        for q in qs:
            for B in BS:
                ef = [e_of(R[f], B, q, kind) for f in R]
                tab[f'{kind}_q{q:g}_B{B}'] = float(np.median(ef))
    # geometric ratios (median over points and f)
    geo = {}
    for B in BS:
        geo[f'rB{B}'] = float(np.median([np.median([r[f'NB{B}'] / r['N1'] for r in R[f]]) for f in R]))
        geo[f'leak{B}'] = float(np.median([np.median([r[f'LB{B}'] / r['S'] for r in R[f]]) for f in R]))
    geo['t'] = float(np.median([np.median([r['Nt_l'] / r['N1'] for r in R[f]]) for f in R]))
    geo['t_by_f'] = {f: float(np.median([r['Nt_l'] / r['N1'] for r in R[f]])) for f in R}
    res[f'v{int(v)}'] = dict(e=tab, geo=geo)
    print(f'--- v {v:.0f}: t = N_true/N_1 (l30) {geo["t"]:.2f};  by f',
          {int(k): round(x, 2) for k, x in geo['t_by_f'].items()})
    print('    r_B = N_B/N_1:', {B: round(geo[f"rB{B}"], 2) for B in BS},
          '  leak L_B/S:', {B: round(geo[f"leak{B}"], 4) for B in BS})
    for key in sorted(tab):
        pass
    for kind in ['m30', 'm40', 'l30']:
        for q in sorted({float(k.split('_q')[1].split('_')[0]) for k in tab if k.startswith(kind)}):
            print(f'    {kind} q {q:g}: ' + '  '.join(f'B{B} {tab[f"{kind}_q{q:g}_B{B}"]:+.2f}' for B in BS))
os.makedirs(os.path.join(ROOT, 'analysis'), exist_ok=True)
json.dump(res, open(os.path.join(ROOT, 'analysis', 't_blockflip.json'), 'w', encoding='utf-8'), indent=1)
