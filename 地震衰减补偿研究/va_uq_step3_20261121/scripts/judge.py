"""Formal synthetic test of the residual Delta t* metric (PLAN_ZH.md section 5), written before any
formal number. Reads analysis/formal_<model>.json (metric.py, default windows) for D4, D1, D3 and
analysis/op_pred.json; writes analysis/formal_judge.json and prints the gates.
usage: python -I judge.py <analysis dir>"""
import sys, os, json
import numpy as np

A = sys.argv[1]
OP = json.load(open(os.path.join(A, 'op_pred.json')))['synthetic']
BAND = 'b8_40'
MAIN = {'D4': ['VA', 'UQW1_25'], 'D1': ['UQW1'], 'D3': ['UQW1']}
SECOND = {'D4': ['VAfix', 'UQW1_100'], 'D1': [], 'D3': []}
CTRL = {'D4': 'PH25', 'D1': 'PH', 'D3': 'PH'}
OPKEY = {'VA': 'VA', 'UQW1_25': 'UQW1_25', 'VAfix': 'VAfix', 'UQW1_100': 'UQW1_100', 'UQW1': 'UQW1_100'}


def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


J = {'band': BAND, 'models': {}}
allrel, all30, allctrl, allzero = True, True, True, True
for mdl in ('D4', 'D1', 'D3'):
    fn = os.path.join(A, 'formal_%s.json' % mdl)
    if not os.path.isfile(fn):
        J['models'][mdl] = 'missing'
        allrel = all30 = allctrl = allzero = False
        continue
    R = json.load(open(fn))
    T = R['arms']['main']
    m = np.array(R['m'])
    r0, r0se = T['R0'][BAND]['dt_lin_s'][0], T['R0'][BAND]['se_lin_s'][0]
    G = {'R0_dt_s': r0, 'R0_se_s': r0se, 'dt_vert_s': R['dtstar_true_vert_s'],
         'S1_R0_positive': bool(r0 - 1.645 * r0se > 0)}
    i1 = int(np.argmin(abs(m - 1)))
    opm = OP[mdl]['m']
    for q in MAIN[mdl] + SECOND[mdl]:
        res, se = np.array(T[q][BAND]['dt_lin_s']), np.array(T[q][BAND]['se_lin_s'])
        rho = spearman(m, res)
        a_obs = (r0 - res) / r0
        fpred = OP[mdl]['band_8_40'][OPKEY[q]]['applied_fraction']
        # op_pred uses m = 0.5 0.7 1 1.4 2; compare at m = 1
        g = dict(res_s=res.tolist(), se_s=se.tolist(), spearman=rho,
                 S2_monotone=bool(rho >= 0.999) if m.size == 5 else bool(rho >= 0.95),
                 S3_direction=bool(res[0] < res[-1] and np.all(res < r0)),
                 a_obs_m1=float(a_obs[i1]), a_pred_m1=float(fpred[opm.index(1.0)]),
                 S4_operator_descriptive=bool(abs(a_obs[i1] - fpred[opm.index(1.0)]) <= 0.2),
                 S7_near_zero_m1=bool(abs(res[i1]) <= 0.25 * r0), res_over_R0_m1=float(res[i1] / r0),
                 S5_precision=bool(se[i1] <= 0.107 * r0), se_over_R0_m1=float(se[i1] / r0),
                 main=q in MAIN[mdl])
        G[q] = g
        if q in MAIN[mdl]:
            allrel &= G['S1_R0_positive'] and g['S2_monotone'] and g['S3_direction']
            all30 &= g['S5_precision']
            allzero &= g['S7_near_zero_m1']
    c = np.array(T[CTRL[mdl]][BAND]['dt_lin_s'])
    G['S6_phase_control'] = bool(np.all(abs(c - r0) <= 0.15 * r0))
    G['S6_max_dev_over_R0'] = float(np.max(abs(c - r0)) / r0)
    allctrl &= G['S6_phase_control']
    FL = R['arms'].get('fieldlike_b6_23_single_well', {})
    r0full = T['R0']['b6_23']['dt_lin_s'][0]
    fl = {}
    for tie in ('exact', 'rho063'):
        v = [FL['x%d_%s' % (x0, tie)]['R0']['dt_lin_s'][0] for x0 in (600, 1200, 1800)]
        s = [FL['x%d_%s' % (x0, tie)]['R0']['se_lin_s'][0] for x0 in (600, 1200, 1800)]
        fl[tie] = dict(R0_dt_s=v, R0_se_s=s, FL1_se_large=bool(max(s) >= 0.3 * r0full),
                       FL2_spread_large=bool(max(v) - min(v) >= 0.5 * r0full))
    G['fieldlike'] = dict(R0_full_b6_23_s=r0full, **fl)
    J['models'][mdl] = G
J['verdict'] = dict(relative_valid=bool(allrel and allctrl), relative_without_control=bool(allrel),
                    near_zero_at_m1=bool(allzero), resolves_30pct=bool(all30),
                    field_run_allowed=bool(allrel and allctrl and all30))
json.dump(J, open(os.path.join(A, 'formal_judge.json'), 'w'), indent=1)
for mdl, G in J['models'].items():
    if G == 'missing':
        print(mdl, 'missing')
        continue
    print('%s R0 %.5f (se %.5f, vertical %.5f) S1 %s  S6 control %s (max dev %.2f R0)' %
          (mdl, G['R0_dt_s'], G['R0_se_s'], G['dt_vert_s'], G['S1_R0_positive'], G['S6_phase_control'], G['S6_max_dev_over_R0']))
    for q in MAIN[mdl] + SECOND[mdl]:
        g = G[q]
        print('  %-9s%s res %s  S2 %s S3 %s S7 %s (res/R0 %.2f) S5 %s (se/R0 %.2f) | S4 desc %s (a %.2f vs %.2f)' %
              (q, '*' if g['main'] else ' ', ' '.join('%+.5f' % v for v in g['res_s']), g['S2_monotone'], g['S3_direction'],
               g['S7_near_zero_m1'], g['res_over_R0_m1'], g['S5_precision'], g['se_over_R0_m1'],
               g['S4_operator_descriptive'], g['a_obs_m1'], g['a_pred_m1']))
    print('  field-like', {k: (v['FL1_se_large'], v['FL2_spread_large']) for k, v in G['fieldlike'].items() if isinstance(v, dict)})
print('verdict', J['verdict'])
