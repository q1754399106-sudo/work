"""Decisions D1-D5 of PLAN_ZH.md, window-level index checks and figures. Reads stage outputs only."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines
import matplotlib.patches

ROOT = Path(__file__).resolve().parents[1]
A, F = ROOT / 'analysis', ROOT / 'figures'
DIAG = ROOT.parent / 'theory_predict_diagnostic_20261024'
plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'], 'axes.unicode_minus': False,
                     'font.size': 9.5, 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
BLUE, VERM, GREEN, PURPLE, BLACK, GREY = '#0072B2', '#D55E00', '#009E73', '#CC79A7', '#000000', '#777777'
SCENES = [('F1', 'Moderate'), ('F2', 'Moderate'), ('F1', 'Strong'), ('F2', 'Strong')]
BASE, SCAN = 'c=1.0000', 'QScan'
NAME = {BASE: 'FullCalRefit (c=1)', SCAN: 'QScan'}


def rho(x, y):
    x, y = pd.Series(np.asarray(x)).rank(), pd.Series(np.asarray(y)).rank()
    return float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else float('nan')


def auc(score, label):
    score, label = np.asarray(score, float), np.asarray(label, bool)
    ok = np.isfinite(score)
    score, label = score[ok], label[ok]
    n1, n0 = label.sum(), (~label).sum()
    if n1 == 0 or n0 == 0:
        return float('nan')
    r = pd.Series(score).rank().to_numpy()
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def save(fig, name, caption):
    for ext in ['png', 'pdf']:
        fig.savefig(F / f'{name}.{ext}', dpi=200)
    (F / f'{name}.txt').write_text(caption, encoding='utf-8')
    plt.close(fig)


def main():
    dec = pd.read_csv(A / 'decomposition.csv')
    pair = pd.read_csv(A / 'qscan_vs_base.csv')
    sel = pd.read_csv(A / 'qscan_selection.csv')
    win = pd.read_csv(A / 'window_moments.csv')
    prof = pd.read_csv(A / 'time_profile.csv')
    spec = pd.read_csv(A / 'spectral_transfer.csv')
    slope = pd.read_csv(A / 'residual_attenuation.csv')
    rep = pd.read_csv(A / 'replicates.csv')
    cobj = pd.read_csv(A / 'clean_qscan_objective.csv')
    robj = pd.read_csv(A / 'replicate_qscan_objective.csv')
    risk = pd.read_csv(DIAG / 'predictions' / 'risk.csv')
    scored = pd.read_csv(DIAG / 'analysis' / 'scored_windows.csv')
    qwin = pd.read_csv(DIAG / 'analysis' / 'qscan_windows.csv')
    qdiag = pd.read_csv(DIAG / 'analysis' / 'qscan_diagnosis.csv')
    key = ['model', 'condition', 'trace', 'window']
    w = win.merge(risk[key + ['live', 'U_r', 'U_Q_f_phase', 'U_cal', 'U_scan']], on=key, how='left', validate='one_to_one')
    w = w.merge(scored[key + ['base_mse', 'scan_mse']].rename(columns={'base_mse': 'arch_base_mse', 'scan_mse': 'arch_scan_mse'}),
                on=key, how='left', validate='one_to_one')
    w = w.merge(qwin[key + ['predicted_log_risk_ratio']], on=key, how='left', validate='one_to_one')
    assert w[['U_r', 'U_scan', 'arch_base_mse']].notna().all().all() and len(w) == 4 * 35 * 26
    for m in ['base', 'qscan']:
        w[f'std_{m}'] = np.sqrt(w[f'var_{m}'])
        w[f'absbias_{m}'] = np.sqrt(w[f'bias2_{m}'].clip(lower=0))
        w[f'rmse_{m}'] = np.sqrt(w[f'mse_{m}'])
    w['expected_degraded'] = w.mean_delta_mse > 0
    w['archived_degraded'] = w.arch_scan_mse > w.arch_base_mse

    decisions, checks, agree = {}, [], []
    for model, cond in SCENES:
        p = pair[(pair.model == model) & (pair.condition == cond)].iloc[0]
        d = dec[(dec.model == model) & (dec.condition == cond)]
        r = rep[(rep.model == model) & (rep.condition == cond)]
        diff = r.rel_QScan - r.rel_c1p0000
        same = max((diff > 0).mean(), (diff < 0).mean())
        fixed = d[d.estimator != SCAN]
        b_opt = float(fixed.loc[fixed.bias2_rel.idxmin(), 'factor'])
        m_opt = float(fixed.loc[fixed.mse_rel.idxmin(), 'factor'])
        c_opt = float(fixed.loc[fixed.noise_free_rel_error.idxmin(), 'factor'])
        g = w[(w.model == model) & (w.condition == cond)]
        vb, vs = d[d.estimator == BASE].var_rel.item(), d[d.estimator == SCAN].var_rel.item()
        pred = qdiag[(qdiag.model == model) & (qdiag.condition == cond)].predicted_log_risk_ratio.item()
        row = dict(R=int(p.R), mse_base=p.mse_base, mse_qscan=p.mse_qscan, delta_mse=p.delta_mse, delta_bias2=p.delta_bias2,
                   delta_var=p.delta_var, delta_mse_over_base=p.delta_mse_over_base,
                   delta_bias2_over_base=p.delta_bias2_over_base, delta_var_over_base=p.delta_var_over_base,
                   fraction_degraded=p.fraction_degraded, same_sign_fraction=same,
                   D1_systematic=bool(same >= .9), median_ratio=p.median_ratio, p05_ratio=p.p05_ratio, p95_ratio=p.p95_ratio,
                   archived_ratio=p.archived_ratio, fraction_selected_1p2803=p.fraction_selected_1p2803,
                   noise_free_selected_factor=p.noise_free_selected_factor,
                   variance_down_bias_up_more=bool(p.delta_var < 0 and p.delta_bias2 > abs(p.delta_var)),
                   bias_optimal_factor=b_opt, mse_optimal_factor=m_opt, noise_free_optimal_factor=c_opt,
                   archived_predicted_log_risk_ratio=pred, empirical_log_std_ratio=0.5 * float(np.log(vs / vb)),
                   index_sign_matches_variance=bool(np.sign(pred) == np.sign(np.log(vs / vb))),
                   index_sign_matches_mse=bool(np.sign(pred) == np.sign(p.delta_mse)),
                   rho_Ur_std_base=rho(g.U_r, g.std_base), rho_Ur_absbias_base=rho(g.U_r, g.absbias_base),
                   rho_Uscan_std_qscan=rho(g.U_scan, g.std_qscan), rho_Uscan_absbias_qscan=rho(g.U_scan, g.absbias_qscan))
        row['D5_index_tracks_variance'] = bool(row['rho_Ur_std_base'] > row['rho_Ur_absbias_base']
                                               and row['rho_Uscan_std_qscan'] > row['rho_Uscan_absbias_qscan']
                                               and row['index_sign_matches_variance'])
        decisions[f'{model}/{cond}'] = row
        for est, u in [('base', 'U_r'), ('base', 'U_cal'), ('qscan', 'U_scan')]:
            for target in ['std', 'absbias', 'rmse']:
                checks.append(dict(model=model, condition=cond, estimator=est, index=u, target=target,
                                   statistic='spearman', value=rho(g[u], g[f'{target}_{est}'])))
            ratio = g[f'std_{est}'] / g[u]
            for q in [.25, .5, .75]:
                checks.append(dict(model=model, condition=cond, estimator=est, index=u, target='std/index',
                                   statistic=f'quantile_{q:.2f}', value=float(ratio.quantile(q))))
        bias_share = (g.bias2_base.clip(lower=0).sum() / g.mse_base.sum(), g.bias2_qscan.clip(lower=0).sum() / g.mse_qscan.sum())
        agree.append(dict(model=model, condition=cond, windows=len(g),
                          expected_degraded_fraction=g.expected_degraded.mean(),
                          archived_degraded_fraction=g.archived_degraded.mean(),
                          label_agreement=(g.expected_degraded == g.archived_degraded).mean(),
                          auroc_index_vs_archived_label=auc(g.predicted_log_risk_ratio, g.archived_degraded),
                          auroc_index_vs_expected_label=auc(g.predicted_log_risk_ratio, g.expected_degraded),
                          spearman_index_vs_expected_delta=rho(g.predicted_log_risk_ratio, g.mean_delta_mse),
                          spearman_index_vs_log_var_ratio=rho(g.predicted_log_risk_ratio, np.log(g.var_qscan / g.var_base)),
                          window_bias_share_base=bias_share[0], window_bias_share_qscan=bias_share[1]))
    moderate = [decisions['F1/Moderate'], decisions['F2/Moderate']]
    decisions['D2_F1_Moderate_variance_down_bias_up_more'] = decisions['F1/Moderate']['variance_down_bias_up_more']
    decisions['D3_F2_Moderate_improvement_from_bias_reduction'] = bool(decisions['F2/Moderate']['delta_bias2'] < 0)
    decisions['D4_opposite_effect_from_different_bias_optimal_compensation'] = bool(
        moderate[0]['bias_optimal_factor'] <= 1.149 and moderate[1]['bias_optimal_factor'] >= 1.2803)
    decisions['D5_index_tracks_variance_both_moderate'] = bool(all(m['D5_index_tracks_variance'] for m in moderate))
    (A / 'decisions.json').write_text(json.dumps(decisions, indent=2), encoding='utf-8')
    pd.DataFrame(checks).to_csv(A / 'window_index_check.csv', index=False)
    pd.DataFrame(agree).to_csv(A / 'window_label_agreement.csv', index=False)
    w.to_csv(A / 'window_joined.csv', index=False)

    # Figure 1: bias^2 + variance stacks (bars start at zero).
    fig, ax = plt.subplots(figsize=(8.6, 4.2), layout='constrained')
    for i, (model, cond) in enumerate(SCENES):
        d = dec[(dec.model == model) & (dec.condition == cond)].set_index('estimator')
        for k, est in enumerate([BASE, SCAN]):
            x = i * 3 + k
            b2, va = d.loc[est, 'bias2_rel'], d.loc[est, 'var_rel']
            col = BLUE if est == BASE else VERM
            ax.bar(x, b2, color=col, edgecolor=BLACK, linewidth=.6)
            ax.bar(x, va, bottom=b2, color='white', edgecolor=col, hatch='///', linewidth=.9)
            ax.plot(x, d.loc[est, 'noise_free_rel_error'] ** 2, marker='D', color=BLACK, ms=5, ls='none')
            ax.plot(x, d.loc[est, 'archived_rel_error'] ** 2, marker='o', mfc='none', mec=PURPLE, mew=1.4, ms=7, ls='none')
            ax.text(x, b2 + va + .004, f'{b2 + va:.3f}', ha='center', va='bottom', fontsize=7.5)
        ax.text(i * 3 + .5, -.085, f'{model}/{cond}', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8.5)
    ax.set_xticks([i * 3 + k for i in range(4) for k in range(2)], ['FullCalRefit', 'QScan'] * 4, fontsize=7.5)
    ax.set_ylabel(r'相对均方误差  $E\|\hat{x}-x^*\|^2/\|x^*\|^2$')
    ax.set_ylim(0, None)
    proxies = [matplotlib.patches.Patch(facecolor=GREY, edgecolor=BLACK, label='偏差²（实心）'),
               matplotlib.patches.Patch(facecolor='white', edgecolor=GREY, hatch='///', label='方差（斜线）'),
               matplotlib.lines.Line2D([], [], marker='D', color=BLACK, ls='none', label='无噪声对照误差²'),
               matplotlib.lines.Line2D([], [], marker='o', mfc='none', mec=PURPLE, mew=1.4, ls='none', label='归档单次实现误差²')]
    ax.legend(handles=proxies, fontsize=7.5, ncols=2, loc='upper left')
    ax.set_title('每个估计量 R=100 个噪声实现：均方误差 = 偏差² + 方差（蓝=FullCalRefit，橙=QScan）', fontsize=9)
    save(fig, 'bias_variance_decomposition',
         '实心柱为无偏校正后的偏差²，斜线柱为方差，二者相加严格等于经验 MSE；均除以目标能量。菱形为无噪声 D0、相同 λ 的确定性误差²；'
         '空心圆为 bae_q 归档单次实现的误差²。R=100 个新噪声种子，配对设计。')

    # Figure 2: factor sweep.
    fig, axs = plt.subplots(2, 2, figsize=(9.2, 6.4), layout='constrained', sharex=True)
    for ax, (model, cond) in zip(axs.flat, SCENES):
        d = dec[(dec.model == model) & (dec.condition == cond) & (dec.estimator != SCAN)].sort_values('factor')
        p = pair[(pair.model == model) & (pair.condition == cond)].iloc[0]
        ax.plot(d.factor, d.mse_rel, color=BLACK, marker='o', ms=3.5, label='MSE')
        ax.plot(d.factor, d.bias2_rel, color=BLUE, ls='--', marker='s', ms=3.5, label='偏差²')
        ax.plot(d.factor, d.var_rel, color=VERM, ls=':', marker='^', ms=4, label='方差')
        ax.plot(d.factor, d.noise_free_rel_error ** 2, color=GREY, lw=1, label='无噪声误差²')
        ax.axvline(1, color=GREY, lw=.8)
        ax.axvline(1.2803, color=GREEN, lw=1.2, ls='-.')
        ax.text(1.2803, .97, f' QScan 选中 {p.fraction_selected_1p2803:.0%}', transform=ax.get_xaxis_transform(),
                color=GREEN, fontsize=7.5, va='top')
        ax.set_xscale('log')
        ax.set_ylim(0, min(.6, d.mse_rel.max() * 1.05))
        ax.set_title(f'{model}/{cond}', fontsize=9)
        ax.set_xticks([.6, .8, 1, 1.28, 1.6, 2.2], ['0.6', '0.8', '1', '1.28', '1.6', '2.2'])
        ax.minorticks_off()
    for ax in axs[1]:
        ax.set_xlabel('施加的 Q 倍率 c（对数轴）')
    for ax in axs[:, 0]:
        ax.set_ylabel('相对误差能量')
    axs[0, 0].legend(fontsize=7.5)
    save(fig, 'factor_sweep',
         '同一批 100 个噪声实现上，FullCalRefit 在固定 c 下的 MSE、偏差²与方差（不做选择）。灰线为无噪声 D0、相同 λ 的误差²。'
         '竖实线 c=1 为公开 Q0，绿色点划线为 QScan 所选 1.2803，文字为选中比例。纵轴截断于 0.6，只截掉极端过补偿端。')

    # Figure 3: paired per-replicate ratio.
    fig, ax = plt.subplots(figsize=(7.2, 3.6), layout='constrained')
    rng = np.random.default_rng(20261024)
    for i, (model, cond) in enumerate(SCENES):
        r = rep[(rep.model == model) & (rep.condition == cond)]
        ratio = r.rel_QScan / r.rel_c1p0000
        ax.boxplot(ratio, positions=[i], widths=.45, showfliers=False, medianprops=dict(color=BLACK))
        ax.scatter(i + rng.uniform(-.17, .17, len(ratio)), ratio, s=7, color=BLUE if cond == 'Moderate' else VERM, alpha=.6,
                   marker='o' if cond == 'Moderate' else '^')
        arch = pair[(pair.model == model) & (pair.condition == cond)].archived_ratio.item()
        ax.plot(i + .32, arch, marker='*', color=PURPLE, ms=10, ls='none', label='归档单次实现' if i == 0 else None)
    ax.axhline(1, color=GREY, lw=.8)
    ax.set_xticks(range(4), [f'{m}/{c}' for m, c in SCENES])
    ax.set_ylabel('配对误差比  rel(QScan) / rel(FullCalRefit)')
    ax.legend(fontsize=8, loc='upper right')
    save(fig, 'paired_ratio',
         '每点为一个噪声实现上 QScan 与 FullCalRefit 的相对误差比（同一观测）；箱线为四分位，须为 1.5 IQR，未画离群标记但点全部保留。'
         '紫色星为 bae_q 归档实现的比值。>1 表示 QScan 更差。')

    # Figure 4: where the change happens in time.
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.6), layout='constrained', sharey=True)
    for ax, (model, cond) in zip(axs, SCENES[:2]):
        p = prof[(prof.model == model) & (prof.condition == cond)]
        b, s = p[p.estimator == BASE].set_index('start_s'), p[p.estimator == SCAN].set_index('start_s')
        mid = (b.index.to_numpy() + b.end_s.to_numpy()) / 2
        ax.bar(mid - .022, s.bias2_rel - b.bias2_rel, width=.04, color=BLUE, label='Δ偏差²')
        ax.bar(mid + .022, s.var_rel - b.var_rel, width=.04, color='white', edgecolor=VERM, hatch='///', label='Δ方差')
        ax.plot(mid, (s.bias2_rel + s.var_rel) - (b.bias2_rel + b.var_rel), color=BLACK, marker='o', ms=3, label='ΔMSE')
        ax.axhline(0, color=GREY, lw=.8)
        ax.set_title(f'{model}/{cond}：QScan − FullCalRefit（100 ms 窗，全部 99 道）', fontsize=9)
        ax.set_xlabel('双程时间 (s)')
    axs[0].set_ylabel('变化量 / 全道目标能量')
    axs[0].legend(fontsize=8)
    save(fig, 'time_profile',
         '按 100 ms 时间窗汇总全部 99 道的偏差²与方差变化（QScan 减 FullCalRefit），各窗贡献相加等于全道 ΔMSE。数据：analysis/time_profile.csv。')

    # Figure 5: spectral transfer of the replicate-mean output.
    bands = [(.3, .8), (.8, 1.2), (1.2, 1.6)]
    fig, axs = plt.subplots(2, 3, figsize=(10.4, 5.6), layout='constrained', sharex=True, sharey=True)
    for row, (model, cond) in enumerate(SCENES[:2]):
        for col, (b0, b1) in enumerate(bands):
            ax = axs[row, col]
            s = spec[(spec.model == model) & (spec.condition == cond) & np.isclose(spec.band_start_s, b0)]
            for est, color, ls in [(BASE, BLUE, '-'), (SCAN, VERM, '--')]:
                z = s[(s.source == 'replicate_mean') & (s.estimator == est)]
                ax.plot(z.frequency_hz, z.transfer_abs, color=color, ls=ls, label=f'{NAME[est]} 平均输出')
            z = s[(s.source == 'noise_free') & (s.estimator == 'c=1.2803')]
            ax.plot(z.frequency_hz, z.transfer_abs, color=VERM, lw=.8, alpha=.6, label='c=1.28 无噪声')
            z = s[(s.source == 'noise_free') & (s.estimator == BASE)]
            ax.plot(z.frequency_hz, z.transfer_abs, color=BLUE, lw=.8, alpha=.6, label='c=1 无噪声')
            ax.axhline(1, color=GREY, lw=.8)
            ax.axvspan(10, 60, color=GREY, alpha=.08)
            ax.set_xlim(0, 90)
            ax.set_ylim(0, 1.6)
            ax.set_title(f'{model}/{cond}  {b0:.2f}–{b1:.2f} s', fontsize=9)
    for ax in axs[1]:
        ax.set_xlabel('频率 (Hz)')
    for ax in axs[:, 0]:
        ax.set_ylabel('|H(f)|：目标→输出')
    axs[0, 0].legend(fontsize=7)
    save(fig, 'spectral_transfer',
         '|H(f)| 为目标到输出的最小二乘传递函数幅度（Hann 窗，全部活道合并）：=1 为振幅谱恢复正确，<1 为补偿不足，>1 为过补偿。'
         '粗线为 100 个实现的平均输出（只含偏差），细线为无噪声对照。灰带 10–60 Hz 为残余吸收斜率拟合区间。')

    # Figure 6: local sensitivity index versus empirical std and |bias| (FullCalRefit windows).
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0), layout='constrained')
    for ax, (model, cond) in zip(axs, SCENES[:2]):
        g = w[(w.model == model) & (w.condition == cond)]
        floor = 1e-6 * g.U_r.median()
        ok, deg = g[g.U_r >= floor], g[g.U_r < floor]
        x0 = .02 * g.U_r.median()   # plotting position for windows whose U_r is numerically zero
        ax.scatter(ok.U_r, ok.std_base, s=8, color=BLUE, marker='o', alpha=.55,
                   label=f'经验标准差  ρ={rho(g.U_r, g.std_base):.2f}')
        ax.scatter(ok.U_r, ok.absbias_base.clip(lower=1e-5), s=9, color=VERM, marker='^', alpha=.55,
                   label=f'|偏差|  ρ={rho(g.U_r, g.absbias_base):.2f}')
        ax.scatter(np.full(len(deg), x0), deg.std_base, s=14, facecolors='none', edgecolors=BLUE, marker='o',
                   label=f'U_r≈0 的窗（{len(deg)} 个）的标准差')
        ax.axvline(x0 * 1.6, color=GREY, lw=.6, ls=':')
        hi = g.U_r.max()
        ax.plot([floor * 1e4, hi], [floor * 1e4, hi], color=GREY, lw=.8, label='y = U_r')
        ax.set(xscale='log', yscale='log', xlabel='归档局部敏感性 U_r（窗 RMS；虚线左侧为 U_r≈0）', ylabel='窗 RMS（波形单位）')
        ax.set_xlim(x0 / 2, hi * 2)
        ax.set_ylim(1e-5, None)
        ax.set_title(f'{model}/{cond}：FullCalRefit，910 窗', fontsize=9)
        ax.legend(fontsize=7, loc='lower right')
    save(fig, 'index_vs_empirical',
         '横轴为 theory_predict_diagnostic 封存的仅反射系数敏感性 U_r；纵轴为 100 个实现的逐窗经验标准差（蓝）与无偏偏差幅度（橙，'
         '负的无偏估计截为 1e-5 只为画在对数轴上，统计未截断）。U_r 小于场景中位数 1e-6 倍的窗（中心解附近无支持）集中画在虚线左侧，'
         '其经验标准差并不为零。ρ 为全部 910 窗的 Spearman。')

    # Figure 7: QScan objective terms with and without noise.
    fig, axs = plt.subplots(1, 4, figsize=(11.2, 3.4), layout='constrained', sharex=True)
    for ax, (model, cond) in zip(axs, SCENES):
        c = cobj[(cobj.model == model) & (cobj.condition == cond)].sort_values('factor')
        ax.plot(c.factor, c.objective - c.objective.min(), color=BLACK, marker='o', ms=3, label='无噪声 J')
        ax.plot(c.factor, c.misfit - c.misfit.min(), color=BLUE, ls='--', label='无噪声 拟合项')
        ax.plot(c.factor, c.penalty - c.penalty.min(), color=VERM, ls=':', label='无噪声 ℓ1 项')
        r = robj[(robj.model == model) & (robj.condition == cond)].copy()
        r['rel'] = r.objective - r.groupby('replicate').objective.transform('min')
        q = r.groupby('factor').rel.quantile([.05, .5, .95]).unstack()
        ax.fill_between(q.index, q[.05], q[.95], color=GREEN, alpha=.18, label='含噪 J 的 5–95%')
        ax.plot(q.index, q[.5], color=GREEN, lw=1.2, label='含噪 J 中位数')
        ax.axvline(1, color=GREY, lw=.8)
        ax.set(xscale='log', title=f'{model}/{cond}', xlabel='c（对数轴）')
        ax.set_yscale('symlog', linthresh=1)
        ax.set_ylim(0, None)
        ax.set_xticks([.6, 1, 1.28, 2.2], ['0.6', '1', '1.28', '2.2'])
        ax.minorticks_off()
    axs[0].set_ylabel('目标值减去各自最小值（symlog）')
    axs[0].legend(fontsize=6.8)
    save(fig, 'qscan_objective',
         'QScan 选择目标 J=(½‖D−Bx‖²+λ‖x‖₁)/σ² 及其两项随 c 的变化，均减去各自最小值。黑/蓝/橙为无噪声 D0（归档 σ），'
         '绿带为 100 个含噪实现。纵轴 symlog（线性阈 1），零值可显示。')
    probe = A / 'snr_probe.csv'
    if probe.is_file():
        snr_figure(dec, pd.read_csv(probe), decisions)
        (A / 'decisions.json').write_text(json.dumps(decisions, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in decisions.items() if not isinstance(v, dict)}, indent=2))


def snr_figure(dec, probe, decisions):
    # Addendum A: QScan minus FullCalRefit across SNR; 8/14 dB from the main run, 20/26 dB from the probe.
    rows = []
    for model in ['F1', 'F2']:
        for cond, snr in [('Strong', 8), ('Moderate', 14)]:
            d = dec[(dec.model == model) & (dec.condition == cond)].set_index('estimator')
            rows.append(dict(model=model, snr_db=snr, R=int(d.loc[BASE, 'R']), **{f'{k}_{e}': d.loc[n, f'{k}_rel']
                        for e, n in [('base', BASE), ('fixed', 'c=1.2803'), ('qscan', SCAN)] for k in ['mse', 'bias2', 'var']}))
        for snr in [20, 26]:
            p = probe[(probe.model == model) & (probe.snr_db == snr)].set_index('estimator')
            rows.append(dict(model=model, snr_db=snr, R=int(p.loc[BASE, 'R']), **{f'{k}_{e}': p.loc[n, f'{k}_rel']
                        for e, n in [('base', BASE), ('fixed', 'c=1.2803'), ('qscan', SCAN)] for k in ['mse', 'bias2', 'var']}))
    t = pd.DataFrame(rows)
    for k in ['mse', 'bias2', 'var']:
        t[f'delta_{k}'] = t[f'{k}_qscan'] - t[f'{k}_base']
    t['delta_mse_over_base'] = t.delta_mse / t.mse_base
    t.to_csv(A / 'snr_trend.csv', index=False)
    f2 = t[t.model == 'F2'].set_index('snr_db')
    decisions['PA_F2_degrades_at_20_and_26dB'] = bool(f2.loc[20, 'delta_mse'] > 0 and f2.loc[26, 'delta_mse'] > 0)
    fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.6), layout='constrained', sharey=True)
    for ax, model in zip(axs, ['F1', 'F2']):
        z = t[t.model == model].sort_values('snr_db')
        ax.plot(z.snr_db, z.delta_mse, color=BLACK, marker='o', label='ΔMSE')
        ax.plot(z.snr_db, z.delta_bias2, color=BLUE, ls='--', marker='s', label='Δ偏差²')
        ax.plot(z.snr_db, z.delta_var, color=VERM, ls=':', marker='^', label='Δ方差')
        ax.axhline(0, color=GREY, lw=.8)
        ax.axvspan(17, 29, color=GREY, alpha=.08)
        ax.set(xlabel='Frobenius SNR (dB)', title=f'{model}：QScan − FullCalRefit', xticks=[8, 14, 20, 26])
    axs[0].set_ylabel('变化量 / 目标能量')
    axs[0].legend(fontsize=8)
    save(fig, 'snr_trend',
         '8/14 dB 来自主分析（R=100），灰底 20/26 dB 来自补充登记 A 的探针（R=40）。ΔMSE=Δ偏差²+Δ方差，均除以目标能量。'
         'ΔMSE>0 表示 QScan 更差。')


if __name__ == '__main__':
    main()
