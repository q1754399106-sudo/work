"""Scoring of the sealed truth-free predictions: decisions E1-E4 of PLAN_ZH.md, tables and figures."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
A, F, PRED = ROOT / 'analysis', ROOT / 'figures', ROOT / 'predictions'
PREV = ROOT / 'cache' / 'inputs' / 'qscan_bias_variance_20261024' / 'analysis'
plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'], 'axes.unicode_minus': False,
                     'font.size': 9.5, 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
BLUE, VERM, GREEN, PURPLE, BLACK, GREY, ORANGE = '#0072B2', '#D55E00', '#009E73', '#CC79A7', '#000000', '#777777', '#E69F00'
# Each candidate's "QScan worse" call, as registered.
CALLS = {'B0': lambda d: d.B0_rel_delta_mse > 0, 'V': lambda d: d.B0_delta_var > 0, 'G': lambda d: d.gate_gain <= 0.01,
         'CV': lambda d: d.cv_rel_delta > 0, 'RS': lambda d: d.rs_delta > 0}
NAMES = {'B0': 'B0 所用Q正确世界自助', 'V': 'V 只看方差', 'G': 'G 第三轮门限', 'CV': 'CV 留出预测', 'RS': 'RS 残差能量'}


def rho(x, y):
    x, y = pd.Series(np.asarray(x)).rank(), pd.Series(np.asarray(y)).rank()
    return float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 and x.std() > 0 and y.std() > 0 else float('nan')


def save(fig, name, caption):
    for ext in ['png', 'pdf']:
        fig.savefig(F / f'{name}.{ext}', dpi=200)
    (F / f'{name}.txt').write_text(caption, encoding='utf-8')
    plt.close(fig)


def expected_truth():
    """Expected relative Delta MSE (QScan - FullCalRefit) per Q0 scene condition."""
    rows = []
    dec = pd.read_csv(PREV / 'decomposition.csv')
    for (m, c), g in dec.groupby(['model', 'condition']):
        b, s = g[g.estimator == 'c=1.0000'].iloc[0], g[g.estimator == 'QScan'].iloc[0]
        rows.append(dict(set='bae', model=m, condition=c, R=int(b.R), mse_base=b.mse_rel, delta=s.mse_rel - b.mse_rel, se=np.nan))
    snr = pd.read_csv(PREV / 'snr_probe.csv')
    for (m, db), g in snr.groupby(['model', 'snr_db']):
        b, s = g[g.estimator == 'c=1.0000'].iloc[0], g[g.estimator == 'QScan'].iloc[0]
        rows.append(dict(set='bae', model=m, condition=f'SNR{db}', R=int(b.R), mse_base=b.mse_rel, delta=s.mse_rel - b.mse_rel, se=np.nan))
    for f in sorted(A.glob('truth_replicates_*.csv')):
        t = pd.read_csv(f)
        d = t.mse_qscan - t.mse_base
        rows.append(dict(set=t.set[0], model=t.model[0], condition=t.condition[0], R=len(t), mse_base=t.mse_base.mean(),
                         delta=d.mean(), se=d.std(ddof=1) / np.sqrt(len(t))))
    e = pd.DataFrame(rows)
    e['rel_delta'] = e.delta / e.mse_base
    e['decisive'] = e.rel_delta.abs() >= 0.05
    return e


def main():
    # Predictions, plan, inputs and prediction-side code must be unchanged since sealing; the
    # scoring-side scripts may be corrected afterwards (any change is logged in logs/).
    seal = json.loads((PRED / 'SEALED_SHA256.json').read_text(encoding='utf-8'))
    scoring = {'scripts/analyze.py', 'scripts/score_diagnostics.m', 'scripts/run_truth_replicates.m', 'scripts/run_stage.py'}
    for rel, digest in seal['hashes'].items():
        if rel in scoring:
            continue
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest, f'sealed file changed: {rel}'
    man = pd.read_csv(PRED / 'manifest.csv')
    diag = pd.concat([pd.read_csv(f) for f in sorted(PRED.glob('diagnostics_w*.csv'))])
    real = pd.read_csv(A / 'realized_scores.csv')
    d = man.merge(diag, on='id', validate='one_to_one').merge(real, on='id', validate='one_to_one')
    d['realized_rel_delta'] = (d.mse_qscan - d.mse_base) / d.mse_base
    exp = expected_truth()
    exp.to_csv(A / 'expected_truth.csv', index=False)
    d = d.merge(exp[['set', 'model', 'condition', 'rel_delta', 'decisive', 'delta', 'mse_base']].rename(
        columns={'rel_delta': 'expected_rel_delta', 'delta': 'expected_delta', 'mse_base': 'expected_mse_base'}),
        on=['set', 'model', 'condition'], how='left')
    q0 = d.q_err == 0
    assert d.loc[q0, 'expected_rel_delta'].notna().all(), 'expected truth missing for a Q0 gather'
    d['decisive'] = d.decisive.fillna(False).astype(bool)
    d['truth_worse'] = np.where(q0, d.expected_rel_delta > 0, d.realized_rel_delta > 0)
    d['label_source'] = np.where(q0, 'expected', 'realized single')
    for k, f in CALLS.items():
        d[f'{k}_worse'] = f(d)
        d[f'{k}_correct'] = d[f'{k}_worse'] == d.truth_worse
    d.to_csv(A / 'scored_gathers.csv', index=False)

    # E1: development F1/Moderate (expected worse) and F2/Moderate (expected better).
    f1 = d[(d.set == 'bae') & (d.model == 'F1') & (d.condition == 'Moderate') & q0]
    f2 = d[(d.set == 'bae') & (d.model == 'F2') & (d.condition == 'Moderate') & q0]
    E1 = bool(len(f1) == 3 and len(f2) == 3 and f1.B0_worse.all() and (~f2.B0_worse).all())
    # E2: confirmation Q0, decisive scenes.
    c0 = d[(d.group == 'confirm') & q0 & d.decisive.astype(bool)]
    acc = {k: float(c0[f'{k}_correct'].mean()) if len(c0) else float('nan') for k in CALLS}
    E2 = 'not testable' if len(c0) == 0 else bool(acc['B0'] >= .75 and acc['B0'] > acc['G'] and acc['B0'] > acc['V'])
    summary = []
    for name, sub in [('dev Q0 (all)', d[(d.group == 'dev') & q0]), ('dev Q0 decisive', d[(d.group == 'dev') & q0 & d.decisive.astype(bool)]),
                      ('confirm Q0 (all)', d[(d.group == 'confirm') & q0]), ('confirm Q0 decisive', c0),
                      ('dev Q-error', d[(d.group == 'dev') & ~q0]), ('confirm Q-error', d[(d.group == 'confirm') & ~q0]),
                      ('Q-30 all', d[d.q_err == -30]), ('Q+30 all', d[d.q_err == 30])]:
        row = dict(subset=name, gathers=len(sub), truth_worse=int(sub.truth_worse.sum()))
        for k in CALLS:
            row[f'{k}_accuracy'] = float(sub[f'{k}_correct'].mean()) if len(sub) else float('nan')
            row[f'{k}_calls_worse'] = int(sub[f'{k}_worse'].sum())
        summary.append(row)
    pd.DataFrame(summary).to_csv(A / 'accuracy_summary.csv', index=False)
    # Per-condition view for Q0 (B0 averaged over the gathers of a condition).
    cond = d[q0].groupby(['group', 'set', 'model', 'condition']).agg(
        gathers=('id', 'size'), expected_rel_delta=('expected_rel_delta', 'first'), decisive=('decisive', 'first'),
        B0_mean=('B0_rel_delta_mse', 'mean'), B0_min=('B0_rel_delta_mse', 'min'), B0_max=('B0_rel_delta_mse', 'max'),
        B0_abs_mean=('B0_delta_bias2', 'mean'), gate_mean=('gate_gain', 'mean'), cv_mean=('cv_rel_delta', 'mean'),
        rs_mean=('rs_delta', 'mean'), qscan_factor=('qscan_factor', 'median'),
        realized_mean=('realized_rel_delta', 'mean')).reset_index()
    cond.to_csv(A / 'q0_conditions.csv', index=False)
    decisions = dict(E1_dev_F1M_F2M_discriminated=E1, E2_confirmation=E2, confirmation_decisive_gathers=len(c0),
                     confirmation_accuracy=acc,
                     spearman_B0_vs_expected_dev=rho(cond[cond.group == 'dev'].B0_mean, cond[cond.group == 'dev'].expected_rel_delta),
                     spearman_B0_vs_expected_confirm=rho(cond[cond.group == 'confirm'].B0_mean, cond[cond.group == 'confirm'].expected_rel_delta),
                     spearman_B0_vs_expected_all=rho(cond.B0_mean, cond.expected_rel_delta))
    (A / 'decisions.json').write_text(json.dumps(decisions, indent=2, ensure_ascii=False), encoding='utf-8')

    # Figure 1: B0 versus expected Delta MSE at Q0 (full range, and the sign boundary).
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.8), layout='constrained')
    z = d[q0]
    style = {('dev', 'bae'): (BLUE, 'o', '开发 F1/F2'), ('confirm', 'qc'): (VERM, '^', '确认 QC-F3/F4'), ('confirm', 'r3'): (GREEN, 's', '确认 R3-F3/F4')}
    for (grp, st), (col, mk, lab) in style.items():
        s = z[(z.group == grp) & (z.set == st)]
        axs[0].scatter(s.expected_rel_delta, s.B0_rel_delta_mse, color=col, marker=mk, s=30, alpha=.85, label=lab, edgecolor='white', lw=.4)
        c = cond[(cond.group == grp) & (cond.set == st)]
        axs[1].errorbar(c.expected_rel_delta, c.B0_mean, yerr=[c.B0_mean - c.B0_min, c.B0_max - c.B0_mean], fmt=mk, color=col, ms=6,
                        capsize=3, label=lab)
    for _, r in cond.iterrows():
        if -.15 <= r.expected_rel_delta <= .25 and -.4 <= r.B0_mean <= .4:
            axs[1].annotate(f'{r.model}/{r.condition[:3]}' if r.set == 'bae' else f'{r.set.upper()}-{r.model}/{r.condition[:3]}',
                            (r.expected_rel_delta, r.B0_mean), xytext=(5, -3), textcoords='offset points', fontsize=7.5, color=BLACK)
    for ax in axs:
        ax.axhline(0, color=GREY, lw=.8)
        ax.axvline(0, color=GREY, lw=.8)
        ax.axvspan(-0.05, 0.05, color=GREY, alpha=.12, lw=0)
        ax.set_xlabel('期望 ΔMSE / MSE(FullCalRefit)（真值，重复实现）')
    axs[0].set_yscale('symlog', linthresh=.3)
    axs[0].set(ylabel='B0：自助世界中的 ΔMSE* / MSE*（symlog）', title='全部 Q0 道集（右上、左下为判对）')
    axs[1].set(xlim=(-.15, .25), ylim=(-.4, .4), title='符号边界附近：场景均值，竖线为该场景道集的最小–最大')
    axs[0].legend(fontsize=8, loc='upper left')
    save(fig, 'b0_vs_expected', '横轴为该场景期望 ΔMSE（开发来自上一轮 100/40 个实现，确认来自 12 个新实现），纵轴为只用观测算出的 B0。'
         '左图每点一个道集，纵轴 symlog（线性阈 0.3）；右图放大符号边界。灰带 |横轴|<5% 为不明确场景，不计入 E2。')

    # Figure 2: sign accuracy per candidate and subset.
    S = pd.DataFrame(summary).set_index('subset')
    subsets = ['dev Q0 decisive', 'confirm Q0 decisive', 'Q-30 all', 'Q+30 all']
    fig, ax = plt.subplots(figsize=(8.4, 3.8), layout='constrained')
    w = .16
    for i, (k, col, hatch) in enumerate(zip(CALLS, [BLUE, GREY, ORANGE, GREEN, PURPLE], ['', '..', '//', 'xx', '\\\\'])):
        vals = [S.loc[s, f'{k}_accuracy'] if s in S.index else np.nan for s in subsets]
        ax.bar(np.arange(len(subsets)) + (i - 2) * w, vals, w, color=col, hatch=hatch, edgecolor=BLACK, lw=.4, label=NAMES[k])
    ax.axhline(.5, color=GREY, lw=.8, ls=':')
    zh = {'dev Q0 decisive': '开发 Q0（明确场景）', 'confirm Q0 decisive': '确认 Q0（明确场景）', 'Q-30 all': 'Q 误差 −30%', 'Q+30 all': 'Q 误差 +30%'}
    ax.set_xticks(range(len(subsets)), [f'{zh[s]}\n(n={int(S.loc[s, "gathers"])})' for s in subsets], fontsize=8)
    ax.set(ylim=(0, 1.05), ylabel='判断 QScan 更差/更好的符号正确率')
    ax.legend(fontsize=7.5, ncols=3, loc='upper center', bbox_to_anchor=(.5, 1.25))
    save(fig, 'sign_accuracy', '各候选判断“QScan 是否比 FullCalRefit 更差”的符号正确率。Q0 子集用期望 ΔMSE 作标签，只计明确场景；'
         'Q±30% 子集用单个归档实现的实际 ΔMSE。虚线为 0.5。')

    # Figure 3: candidate values on the Q0 gathers, ordered by expected Delta MSE.
    z = d[q0].sort_values('expected_rel_delta').reset_index(drop=True)
    fig, axs = plt.subplots(4, 1, figsize=(9.0, 7.6), layout='constrained', sharex=True)
    xs = np.arange(len(z))
    cols = np.where(z.truth_worse, VERM, BLUE)
    for ax, (field, lab, thr) in zip(axs, [('B0_rel_delta_mse', 'B0  ΔMSE*/MSE*', 0), ('gate_gain', 'G 目标相对下降', .01),
                                           ('cv_rel_delta', 'CV  Δ留出误差', 0), ('rs_delta', 'RS  Δ残差能量', 0)]):
        ax.scatter(xs, z[field], c=cols, s=22, marker='o')
        ax.axhline(thr, color=GREY, lw=.8, ls='--')
        ax.set_ylabel(lab, fontsize=8)
    axs[-1].set_xticks(xs, [f'{r.set}-{r.model}/{r.condition[:3]}' for r in z.itertuples()], rotation=90, fontsize=6)
    axs[0].set_title('Q0 道集按期望 ΔMSE 由小到大排列；橙=真实更差，蓝=真实更好；虚线为各候选的判断阈值', fontsize=9)
    save(fig, 'candidate_values', '每列一个 Q0 道集。B0 越过 0、CV/RS 越过 0、G 越过 0.01 的方向即各候选的判断。')

    # Figure 4: Q-error gathers.
    z = d[~q0]
    fig, ax = plt.subplots(figsize=(6.8, 4.4), layout='constrained')
    for q, mk in [(-30, 'v'), (30, '^')]:
        s = z[z.q_err == q]
        ax.scatter(s.realized_rel_delta, s.B0_rel_delta_mse, marker=mk, s=40, c=np.where(s.group == 'dev', BLUE, VERM),
                   label=f'Q 误差 {q:+d}%（蓝=开发，橙=确认）', edgecolor='white', lw=.4)
    ax.axhline(0, color=GREY, lw=.8)
    ax.axvline(0, color=GREY, lw=.8)
    ax.set(xlabel='实际 ΔMSE / MSE(FullCalRefit)（单个归档实现）', ylabel='B0  ΔMSE*/MSE*')
    ax.legend(fontsize=8)
    ax.set_title('Q 有误时 B0 的假设不成立：左上为误报（实际 QScan 更好，B0 判更差）', fontsize=9)
    save(fig, 'q_error_cases', 'Q 误差 ±30% 的归档观测。B0 假设所用 Q 正确，因此在 Q 有误、QScan 实际纠正误差时会误报。')
    print(json.dumps(decisions, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
