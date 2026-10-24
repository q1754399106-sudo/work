"""Score the sealed confirmation predictions on U1-U4 against 12-replicate expected truth (PLAN_ZH.md)."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
A, F, PRED = ROOT / 'analysis', ROOT / 'figures', ROOT / 'predictions'
plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'], 'axes.unicode_minus': False,
                     'font.size': 9.5, 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
BLUE, VERM, GREEN, PURPLE, BLACK, GREY, ORANGE = '#0072B2', '#D55E00', '#009E73', '#CC79A7', '#000000', '#777777', '#E69F00'
COLS = {'B0': 'B0_rel_delta_mse', 'B0M': 'B0M_rel_delta_mse', 'B0C': 'B0C_rel_delta_mse'}


def rho(x, y):
    x, y = pd.Series(np.asarray(x)).rank(), pd.Series(np.asarray(y)).rank()
    return float(np.corrcoef(x, y)[0, 1])


def save(fig, name, caption):
    for ext in ['png', 'pdf']:
        fig.savefig(F / f'{name}.{ext}', dpi=200)
    (F / f'{name}.txt').write_text(caption, encoding='utf-8')
    plt.close(fig)


def main():
    seal = json.loads((PRED / 'SEALED_SHA256.json').read_text(encoding='utf-8'))
    scoring = {'scripts/analyze_confirm.py', 'scripts/run_stage.py'}
    for rel, digest in seal['hashes'].items():
        if rel not in scoring:
            assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest, f'sealed file changed: {rel}'
    sel = json.loads((A / 'dev_selection.json').read_text(encoding='utf-8'))
    var, tau = sel['selected']['variant'], sel['selected']['threshold']
    taus = {v: (s['threshold'] if s['separable'] else None) for v, s in sel['variants'].items()}
    d = pd.concat([pd.read_csv(f) for f in sorted(PRED.glob('confirm_w*.csv'))]).sort_values('id')
    t = pd.concat([pd.read_csv(f) for f in sorted(A.glob('truth_U*.csv'))])
    t['delta'] = t.mse_qscan - t.mse_base
    exp = t.groupby(['model', 'condition']).agg(R=('delta', 'size'), mse_base=('mse_base', 'mean'), delta=('delta', 'mean'),
                                                se=('delta', lambda s: s.std(ddof=1) / np.sqrt(len(s))),
                                                worse_fraction=('delta', lambda s: float((s > 0).mean()))).reset_index()
    exp['rel_delta'] = exp.delta / exp.mse_base
    exp['decisive'] = exp.rel_delta.abs() >= 0.05
    exp['boundary'] = exp.decisive & (exp.rel_delta.abs() < 0.20)
    exp.to_csv(A / 'confirm_expected_truth.csv', index=False)
    real = t[t.replicate <= 3].rename(columns={'replicate': 'replicate', 'mse_base': 'real_mse_base', 'mse_qscan': 'real_mse_qscan'})
    d = d.merge(exp[['model', 'condition', 'rel_delta', 'decisive', 'boundary']], on=['model', 'condition'], validate='many_to_one')
    d = d.merge(real[['model', 'condition', 'replicate', 'real_mse_base', 'real_mse_qscan']], on=['model', 'condition', 'replicate'],
                validate='one_to_one')
    d['truth_worse'] = d.rel_delta > 0
    calls = {'selected': d[COLS[var]] > tau, 'B0 (tau=0)': d.B0_rel_delta_mse > 0, 'B0M (tau=0)': d.B0M_rel_delta_mse > 0,
             'B0C (tau=0)': d.B0C_rel_delta_mse > 0, 'V': d.B0_delta_var > 0, 'G': d.gate_gain <= 0.01}
    for v, tv in taus.items():
        if tv is not None:
            calls[f'{v} (dev tau={tv:.3f})'] = d[COLS[v]] > tv
    dec = d.decisive.astype(bool)
    acc = {k: dict(decisive=float((c[dec] == d.truth_worse[dec]).mean()) if dec.any() else None,
                   boundary=float((c[d.boundary] == d.truth_worse[d.boundary]).mean()) if d.boundary.any() else None,
                   all=float((c == d.truth_worse).mean())) for k, c in calls.items()}
    C1 = bool(dec.any() and acc['selected']['decisive'] >= 0.80)
    C2 = bool(dec.any() and all(acc['selected']['decisive'] >= acc[k]['decisive'] for k in ['B0 (tau=0)', 'V', 'G']))
    cond = d.groupby(['model', 'condition']).agg(rel_delta=('rel_delta', 'first'), decisive=('decisive', 'first'),
                                                 **{v: (c, 'mean') for v, c in COLS.items()},
                                                 qscan_factor=('qscan_factor', 'median')).reset_index()
    cond.to_csv(A / 'confirm_conditions.csv', index=False)
    out = dict(selected_variant=var, threshold=tau, decisive_gathers=int(dec.sum()), boundary_gathers=int(d.boundary.sum()),
               C1_accuracy_at_least_0_80=C1, C2_not_worse_than_B0_V_G=C2, accuracy=acc,
               spearman_conditions={v: rho(cond[v], cond.rel_delta) for v in COLS})
    (A / 'confirm_decisions.json').write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding='utf-8')
    for k, c in calls.items():
        d[f'call_{k}'] = c
    d.to_csv(A / 'confirm_scored.csv', index=False)

    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.6), layout='constrained')
    for v, col, mk in [('B0', GREY, 'o'), ('B0M', BLUE, 's'), ('B0C', VERM, '^')]:
        axs[0].scatter(d.rel_delta, d[COLS[v]], color=col, marker=mk, s=26, alpha=.8, label=v)
        c = cond
        axs[1].plot(c.rel_delta, c[v], ls='none', marker=mk, color=col, ms=6, label=v)
    for ax in axs:
        ax.axhline(0, color=GREY, lw=.8)
        ax.axvline(0, color=GREY, lw=.8)
        ax.axvspan(-.05, .05, color=GREY, alpha=.12, lw=0)
        ax.set_xlabel('期望 ΔMSE / MSE(FullCalRefit)（12 个实现）')
    if tau:
        for ax in axs:
            ax.axhline(tau, color={'B0': GREY, 'B0M': BLUE, 'B0C': VERM}[var], lw=1, ls='--')
    axs[0].set_yscale('symlog', linthresh=.3)
    axs[0].set(ylabel='自助世界中的 ΔMSE* / MSE*（symlog）', title='U1–U4 全部 48 个道集')
    axs[1].set(xlim=(-.3, .3), ylim=(-.4, .5), title=f'符号边界附近（场景均值）；虚线为冻结阈值 τ={tau:.3f}' if tau else '符号边界附近（场景均值）')
    axs[0].legend(fontsize=8)
    save(fig, 'confirm_b0_variants', '新模型 U1–U4 上三种世界真值的 B0 与期望 ΔMSE。灰带为不明确场景；虚线为开发阶段冻结的阈值。')
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
