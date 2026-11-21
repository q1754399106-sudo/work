"""Theory checks for DERIVATION_ZH.md, parts 1-2 (no new migration runs).

Part 1: Bayes risk of the expected (UQ) operator for one frequency,
    R = 1 - |E H|^2 S / (E|H|^2 S + n),  H = exp(-u z),  z = alpha - i beta.
  Claim (exact for Gaussian u ~ N(ubar, s^2)):
    R = 1/(1+rho) + (1 - exp(-s^2 |z|^2)) rho/(1+rho),   rho = E|H|^2 S / n.
  Checked against Monte Carlo and against the lognormal prior of the code (41-point Gauss-Hermite).

Part 2: break-even calibration accuracy from the cached stage-3 scores (uq_rtm_20261117/analysis,
  confirmation models C1-C4, deep window, E1). A calibrated Q leaves a residual ln-error
  e ~ N(b, sc^2); its score is approximated by W1 at m = exp(e) (sigma -> 0 limit of a posterior UQ),
  interpolated linearly in ln m. Compared with UQW1 (sigma 0.447) averaged over the m grid and its worst m.
"""
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RTM = os.path.join(os.path.dirname(ROOT), 'uq_rtm_20261117', 'analysis')
OUT = os.path.join(ROOT, 'analysis')
os.makedirs(OUT, exist_ok=True)
res = {}

# ---------------- part 1 ----------------
rng = np.random.default_rng(20261121)
rows = []
for (alpha, beta, s, rho0) in [(1.26, 1.11, .447, 3.), (2.0, 1.2, .447, 1.), (2.0, 1.2, .1, 1.),
                               (3.45, 1.5, .447, 10.), (3.45, 1.5, .1, 10.), (1.0, 0.0, .3, .1)]:
    z = alpha - 1j * beta
    u = 1 + s * rng.standard_normal(2_000_000)           # Gaussian u, ubar = 1
    H = np.exp(-u * z)
    S, n = 1.0, None
    EH = H.mean(); EH2 = (abs(H) ** 2).mean()
    n = EH2 * S / rho0                                     # noise so that rho = rho0
    R_mc = 1 - abs(EH) ** 2 * S / (EH2 * S + n)
    c = np.exp(-s ** 2 * abs(z) ** 2)
    R_th = 1 / (1 + rho0) + (1 - c) * rho0 / (1 + rho0)
    # lognormal u = exp(s e) (code prior), same formula with s -> s (first order)
    e, w = np.polynomial.hermite_e.hermegauss(41); w = w / w.sum()
    ul = np.exp(s * e); Hl = np.exp(-np.outer(ul, [z]).ravel())
    EHl = (w * Hl).sum(); EH2l = (w * abs(Hl) ** 2).sum(); nl = EH2l * S / rho0
    R_ln = 1 - abs(EHl) ** 2 * S / (EH2l * S + nl)
    # oracle (known u) Bayes risk averaged over the same Gaussian prior
    rho_u = abs(H) ** 2 * S / n
    R_or = (1 / (1 + rho_u)).mean()
    rows.append(dict(alpha=alpha, beta=beta, s=s, rho=rho0, R_mc=R_mc, R_theory=R_th,
                     R_lognormal=R_ln, R_oracle=R_or, c=c))
res['part1'] = rows

# ---------------- part 2 ----------------
def interp_lnm(m, y, x):
    return np.interp(np.log(x), np.log(m), y)

grid = dict(sc=[0.05, 0.10, 0.15, 0.20, 0.30, 0.447], b=[-0.15, -0.10, 0.0, 0.10, 0.15])
ee, ww = np.polynomial.hermite_e.hermegauss(31); ww = ww / ww.sum()
part2 = {}
for mod in ['C1', 'C2', 'C3', 'C4']:
    d = json.load(open(os.path.join(RTM, f'stage2_{mod}.json'), encoding='utf-8'))
    for ds in ['m30', 'm40', 'l30']:
        T = d[ds]['deep']['E1']; m = np.array(T['m']); W1 = np.array(T['W1']); U1 = np.array(T['UQW1'])
        Bp = np.maximum.reduce([np.array(T[k]) for k in ['R0', 'PH', 'IQ30', 'IQ70', 'LP', 'W0', 'W1']])
        cell = dict(UQW1_mean=U1.mean(), UQW1_worst=U1.min(), Bplus_mean=Bp.mean(), Bplus_worst=Bp.min(),
                    W1_at1=float(W1[m == 1][0]), UQW1_at1=float(U1[m == 1][0]), cal={})
        for b in grid['b']:
            for sc in grid['sc']:
                x = np.exp(np.clip(b + sc * ee, np.log(m[0]), np.log(m[-1])))   # clip to the grid range
                cell['cal'][f'b{b:+.2f}_s{sc:.3f}'] = float((ww * interp_lnm(m, W1, x)).sum())
        part2[f'{mod}_{ds}'] = cell
res['part2'] = part2

# summary: break-even sc (b = 0) where E[W1] = UQW1 mean over the grid
summ = {}
for ds in ['m30', 'm40', 'l30']:
    keys = [k for k in part2 if k.endswith(ds)]
    summ[ds] = {
        'UQW1_mean': np.mean([part2[k]['UQW1_mean'] for k in keys]),
        'UQW1_worst': np.mean([part2[k]['UQW1_worst'] for k in keys]),
        'Bplus_mean': np.mean([part2[k]['Bplus_mean'] for k in keys]),
        'W1_at1': np.mean([part2[k]['W1_at1'] for k in keys]),
        'UQW1_at1': np.mean([part2[k]['UQW1_at1'] for k in keys]),
        'cal': {c: np.mean([part2[k]['cal'][c] for k in keys]) for c in part2[keys[0]]['cal']},
    }
res['summary'] = summ
json.dump(res, open(os.path.join(OUT, 't_bayes.json'), 'w', encoding='utf-8'), indent=1, default=float)

for r in rows:
    print('P1 a%.2f b%.2f s%.3f rho%.1f  MC %.4f  theory %.4f  lognormal %.4f  oracle %.4f  c %.3f' %
          (r['alpha'], r['beta'], r['s'], r['rho'], r['R_mc'], r['R_theory'], r['R_lognormal'], r['R_oracle'], r['c']))
for ds, s in summ.items():
    print(f"\nP2 {ds}: UQW1 mean {s['UQW1_mean']:.3f} worst {s['UQW1_worst']:.3f} | B+ mean {s['Bplus_mean']:.3f} | "
          f"W1@1 {s['W1_at1']:.3f} UQW1@1 {s['UQW1_at1']:.3f}")
    for b in grid['b']:
        print('   b %+.2f: ' % b + '  '.join('s%.2f %.3f' % (sc, s['cal'][f'b{b:+.2f}_s{sc:.3f}']) for sc in grid['sc']))
