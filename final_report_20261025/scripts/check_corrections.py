"""Reproduce the numbers behind ERRATA_ZH.md from archived results (read via `git show HEAD:...`).

1. CFS coefficient: centroids and variance from the *power* spectrum need 2*pi, not pi
   (Quan & Harris 1997 state the formula for the amplitude spectrum). Ideal Gaussian test,
   then the corrected Q_CFS of experiments 4 and 5 (Q_CFS is linear in the coefficient).
2. fd2d scores are normalised squared errors ||X-T||^2/||T||^2; convert the reported
   relative changes x to relative-error changes sqrt(1+x)-1 and list the noise-free
   penalty at Q0.
Usage: python check_corrections.py <repo root>  -> final_report_20261025/analysis/corrections.json
"""
import csv, io, json, math, subprocess, sys
from pathlib import Path
import numpy as np

root = Path(sys.argv[1])
show = lambda p: subprocess.run(['git', '-C', str(root), 'show', f'HEAD:{p}'], capture_output=True, check=True).stdout.decode('utf-8')
out = {}

# 1. CFS on an ideal Gaussian amplitude spectrum, scored exactly as vsp_q.m (power-spectrum moments).
f = np.linspace(0, 300, 30001); f0, s, Q = 50., 15., 100.; taus = np.array([0., .2, .4, .6])
A = np.exp(-(f - f0) ** 2 / (2 * s * s)); fc = []
for t in taus:
    p = (A * np.exp(-np.pi * f * t / Q)) ** 2; fc.append((f * p).sum() / p.sum())
p = A ** 2; var = ((f - fc[0]) ** 2 * p).sum() / p.sum(); C = np.polyfit(taus, fc[0] - np.array(fc), 1)[0]
out['cfs_gaussian_test'] = dict(true_Q=Q, Q_with_pi=np.pi * var / C, Q_with_2pi=2 * np.pi * var / C)
for name, path, key in [('U1320A', 'field_vsp_u1320a_20261024/analysis/vsp_q.json', 'Q_SR'),
                        ('Volve', 'field_vsp_volve_20261025/analysis/vsp_q.json', 'Q_int')]:
    j = json.loads(show(path)); q = 2 * j['Q_CFS']; r = q / j[key]
    out[f'cfs_{name}'] = dict(archived_Q_CFS=j['Q_CFS'], corrected_Q_CFS=q, reference=key, reference_Q=j[key],
                             archived_ratio=j['Q_CFS'] / j[key], corrected_ratio=r, ratio_in_0p5_2=0.5 <= r <= 2,
                             archived_G1_pass=j['gate']['G1_pass'], median_r2=j['gate']['median_r2'])

# 2. fd2d: squared-error changes -> relative-error changes.
rel = lambda x: math.sqrt(1 + x) - 1
S = list(csv.DictReader(io.StringIO(show('fd2d_20261025/analysis/fd_summary.csv'))))
C2 = list(csv.DictReader(io.StringIO(show('fd2d_20261025/analysis/fd_curves.csv'))))
rows = []
for r in S:
    d = dict(model=r['model'], scene=r['scene'], offset_steps=int(r['offset_steps']),
             q_qscan=float(r['q_qscan_mode']), q_best=float(r['q_best']), q_best_noise_free=float(r['q_best_noise_free']))
    for k in ['rel_B_vs_A', 'rel_B', 'rel_S1', 'rel_S2', 'rel_opt']:
        d[k + '_sq'] = float(r[k]); d[k + '_relerr'] = rel(float(r[k]))
    rows.append(d)
out['fd2d_rows'] = rows
p = [r for r in rows if r['scene'] == 'Q100_14dB']; u = [r for r in rows if r['scene'] == 'Q070_14dB']
rng = lambda v: [min(v), max(v)]
out['fd2d_ranges'] = dict(
    E3_B_vs_A_sq=rng([r['rel_B_vs_A_sq'] for r in p]), E3_B_vs_A_relerr=rng([r['rel_B_vs_A_relerr'] for r in p]),
    E4_Q070_B_vs_D_sq=rng([r['rel_B_sq'] for r in u]), E4_Q070_B_vs_D_relerr=rng([r['rel_B_relerr'] for r in u]),
    S1_worst_sq=max(r['rel_S1_sq'] for r in rows), S1_worst_relerr=max(r['rel_S1_relerr'] for r in rows),
    S2_worst_sq=max(r['rel_S2_sq'] for r in rows), S2_worst_relerr=max(r['rel_S2_relerr'] for r in rows),
    S1_median_sq=float(np.median([r['rel_S1_sq'] for r in rows])), S2_median_sq=float(np.median([r['rel_S2_sq'] for r in rows])))
out['fd2d_ranges']['S1_median_relerr'] = rel(out['fd2d_ranges']['S1_median_sq'])
out['fd2d_ranges']['S2_median_relerr'] = rel(out['fd2d_ranges']['S2_median_sq'])
out['fd2d_D_better_than_B_Q0'] = {r['model']: 1 - 1 / (1 + r['rel_B_sq']) for r in p if r['rel_B_sq'] > 0}
pen = {}
for m in sorted({c['model'] for c in C2}):
    c = [x for x in C2 if x['model'] == m and x['scene'] == 'Q100_14dB']
    nf = [float(x['noise_free_mse']) for x in c]; q = [float(x['q_ratio_to_Q0']) for x in c]
    i0 = min(range(len(q)), key=lambda i: abs(q[i] - 1)); mn = min(nf)
    pen[m] = dict(q_best=q[nf.index(mn)], sq_penalty_at_Q0=nf[i0] / mn - 1, relerr_penalty_at_Q0=math.sqrt(nf[i0] / mn) - 1)
out['fd2d_noise_free_penalty_at_Q0'] = pen

# Squared-error changes quoted in REPORT_ZH.md from other experiments, converted at the range ends
# (sqrt(1+x)-1 is monotone, so per-item ranges convert exactly; means do not and are left as squared).
quoted = {'combo_stage1_Q0_combo_worse': [0.21, 0.50, -0.065], 'combo_stage1_Q070_combo_better': [-0.63, -0.79],
          'combo_stage2_opt_vs_D': [-0.24, -0.33], 'qscan_offset_Q070_B_worse': [0.47, 3.57],
          'qscan_offset_shift_worst': [0.34, 0.90], 'qscan_abstain_A2_worst': [0.141], 'critcombo_A1A3_worst': [0.39, 0.62]}
out['report_squared_to_relerr'] = {k: [dict(sq=x, relerr=rel(x)) for x in v] for k, v in quoted.items()}

dst = root / 'final_report_20261025' / 'analysis'; dst.mkdir(exist_ok=True)
(dst / 'corrections.json').write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding='utf-8')
print(json.dumps({k: v for k, v in out.items() if k != 'fd2d_rows'}, indent=1, ensure_ascii=False))
