"""Descriptive evaluation of sealed predictions; never tunes the predictor."""
from pathlib import Path
import hashlib, json, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'], 'axes.unicode_minus': False,
                     'font.size': 10, 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
FEATURES = ['depth_time', 'inverse_observed_snr', 'q_spread', 'U_r', 'U_Q', 'U_Q_f_phase', 'U_cal', 'U_expanded']
LABELS = {'depth_time':'双程时间', 'inverse_observed_snr':'观测信噪比倒数', 'q_spread':'Q 扰动离散度',
          'U_r':'仅反射系数', 'U_Q':'反射系数 + Q', 'U_Q_f_phase':'+ 主频与相位',
          'U_cal':'+ 独立校准', 'U_expanded':'加入遗漏支持候选'}
COLORS = ['#555555', '#946A00', '#009E73', '#0072B2', '#6A3D9A', '#B34700', '#000000', '#CC79A7']

def rho(x, y):
    x, y = pd.Series(x).rank(), pd.Series(y).rank()
    if x.std() == 0 or y.std() == 0: return float('nan')
    return float(np.corrcoef(x, y)[0, 1])

def auc(score, label):
    score, label = np.asarray(score), np.asarray(label, bool)
    valid = ~np.isnan(score)
    score, label = score[valid], label[valid]
    n1, n0 = label.sum(), (~label).sum()
    if n1 == 0 or n0 == 0: return float('nan')
    ranks = pd.Series(score).rank().to_numpy()
    return float((ranks[label].sum()-n1*(n1+1)/2)/(n1*n0))

def energy_top(score, error, fraction=.2):
    # Fractional membership of a tied cutoff avoids arbitrary row-order wins.
    score, error = np.asarray(score), np.asarray(error)
    order = np.argsort(-score, kind='stable'); score, error = score[order], error[order]
    budget = fraction*len(score); total = 0.
    for value in np.unique(score)[::-1]:
        m = score == value; weight = min(1., budget/m.sum())
        total += weight*error[m].sum(); budget -= weight*m.sum()
        if budget <= 1e-10: break
    return float(total/error.sum()) if error.sum() else float('nan')

def save(fig, name, caption):
    for ext in ['png', 'pdf']: fig.savefig(ROOT/'figures'/f'{name}.{ext}', dpi=180)
    (ROOT/'figures'/f'{name}.txt').write_text(caption, encoding='utf-8')
    plt.close(fig)

def main():
    seal = json.loads((ROOT/'predictions/SEALED_SHA256.json').read_text())
    for rel, digest in seal['hashes'].items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == digest, f'Sealed artifact changed: {rel}'
    d = pd.read_csv(ROOT/'analysis/scored_windows.csv')
    d['scene'] = d.model+'_'+d.condition
    normal = d[d.live.astype(bool)].copy()
    per=[]
    for subset, frame in [('all',d), ('live',normal)]:
        for (model, cond), g in frame.groupby(['model','condition']):
            for f in FEATURES:
                defined=not (f.startswith('U_') and not g.live.astype(bool).all())
                per.append(dict(subset=subset, model=model, condition=cond, feature=f, n=len(g),
                                spearman=rho(g[f],g.error_rms) if defined else np.nan,
                                top20_energy=energy_top(g[f],g.base_mse) if defined else np.nan,
                                valid_domain=defined,
                                infinite_fraction=float(np.isinf(g[f]).mean())))
    per=pd.DataFrame(per);per.to_csv(ROOT/'analysis/predictor_by_scene.csv',index=False)
    domain=d.groupby(['model','condition']).agg(total_windows=('live','size'),valid_windows=('live','sum')).reset_index()
    domain['invalid_reason']='zero observation with inherited sigma=realmin; local uncertainty invalid'
    domain.to_csv(ROOT/'analysis/prediction_domain.csv',index=False)
    macro=per.groupby(['subset','feature'],sort=False)[['spearman','top20_energy','infinite_fraction']].mean().reset_index()
    for idx,row in macro.iterrows():
        coverage=per[(per.subset==row.subset)&(per.feature==row.feature)]
        if not coverage.valid_domain.all():macro.loc[idx,['spearman','top20_energy','infinite_fraction']]=np.nan
    macro.to_csv(ROOT/'analysis/predictor_comparison.csv',index=False)
    # Leave a whole geology out, including both of its noise conditions.
    for f in FEATURES:
        normal[f+'_rank']=normal.groupby('scene')[f].rank(pct=True)
    normal['actual_rank']=normal.groupby('scene').error_rms.rank(pct=True)
    basic=['depth_time','inverse_observed_snr','q_spread'];cv=[];pred=[]
    variants=[('simple',basic),('simple_plus_r',basic+['U_r']),('simple_plus_Q',basic+['U_Q']),
              ('simple_plus_f_phase',basic+['U_Q_f_phase']),('simple_plus_theory',basic+['U_cal']),
              ('simple_plus_expanded',basic+['U_expanded'])]
    for variant,features in variants:
        cols=[f+'_rank' for f in features]
        for held in sorted(normal.model.unique()):
            train=normal[normal.model!=held];test=normal[normal.model==held]
            X=np.column_stack([np.ones(len(train)),train[cols].to_numpy()]);y=train.actual_rank.to_numpy()
            # Fixed numerical stabilizer, no hyperparameter selection.
            reg=np.eye(X.shape[1])*1e-8;reg[0,0]=0
            coef=np.linalg.solve(X.T@X+reg,X.T@y)
            Xtest=np.column_stack([np.ones(len(test)),test[cols].to_numpy()]);p=Xtest@coef
            yy=test.actual_rank.to_numpy();r2=1-np.sum((p-yy)**2)/np.sum((yy-yy.mean())**2)
            cv.append(dict(variant=variant,held_model=held,spearman=rho(p,yy),r2=r2,n=len(test)))
            z=test[['model','condition','trace','window','actual_rank']].copy();z['variant']=variant;z['prediction']=p;pred.append(z)
    cv=pd.DataFrame(cv);cv.to_csv(ROOT/'analysis/leave_model_out.csv',index=False)
    pd.concat(pred).to_csv(ROOT/'analysis/leave_model_out_predictions.csv',index=False)
    cases=pd.read_csv(ROOT/'analysis/case_scores.csv');casepred=[];window_scan=[]
    for (model,condition),g in normal.groupby(['model','condition']):
        finite=np.isfinite(g.U_cal)&np.isfinite(g.U_scan)&(g.U_cal>0)&(g.U_scan>0)
        ratio=np.full(len(g),np.nan);ratio[finite]=np.log(g.loc[finite,'U_scan']/g.loc[finite,'U_cal'])
        change=(g.scan_mse-g.base_mse).to_numpy()
        casepred.append(dict(model=model,condition=condition,predicted_log_risk_ratio=float(np.nanmedian(ratio)) if finite.any() else np.nan,
                             finite_ratio_coverage=float(finite.mean()),observed_window_degraded_fraction=float((change>0).mean()),
                             delta_spearman=rho(ratio[finite],change[finite]) if finite.any() else np.nan,
                             degradation_auc=auc(ratio,change>0),qspread_degradation_auc=auc(g.q_spread,change>0)))
        z=g[['model','condition','trace','window']].copy();z['predicted_log_risk_ratio']=ratio;z['actual_mse_change']=change;window_scan.append(z)
    c=cases.merge(pd.DataFrame(casepred),on=['model','condition']);c.to_csv(ROOT/'analysis/qscan_diagnosis.csv',index=False)
    pd.concat(window_scan).to_csv(ROOT/'analysis/qscan_windows.csv',index=False)
    # Preserve original scores and check agreement with archived old live-trace rule.
    checks=[]
    for _, row in c.iterrows():
        split='dev' if row.model in ['primary','complex'] else 'ind'
        old=pd.read_csv(ROOT.parent/'bae_q_20261024/analysis'/f'{split}_metrics_all.csv')
        sub=old[(old.model==row.model)&(old.condition==row.condition)&(old.Q_percent==0)]
        for method,col in [('FullCalRefit','base_old_live'),('QScanRefit','scan_old_live')]:
            val=float(sub.loc[sub.method==method,'relative_error'].iloc[0]);diff=abs(row[col]-val)
            assert diff<1e-8,(row.model,method,diff)
            checks.append(dict(model=row.model,condition=row.condition,method=method,absolute_difference=diff))
    pd.DataFrame(checks).to_csv(ROOT/'analysis/old_score_reproduction.csv',index=False)
    derivative=pd.read_csv(ROOT/'analysis/derivative_checks.csv');audit=pd.read_csv(ROOT/'analysis/numerical_audit.csv')
    hcheck=derivative[derivative.step==1e-4]
    assert hcheck.observation_relative.max()<1e-5
    assert hcheck.target_relative.max()<1e-5
    assert derivative[derivative.parameter.str.startswith('logQ')].target_absolute.max()<1e-10
    # Figure 1: all derivative curves (maximum across models/conditions).
    fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for parameter,g in derivative.groupby('parameter',sort=False):
        for ax,field in zip(axs,['observation_relative','target_relative']):
            if field=='target_relative' and parameter.startswith('logQ'):continue
            z=g.groupby('step')[field].max();ax.loglog(z.index,np.maximum(z.to_numpy(),1e-14),marker='o',label=parameter)
    for ax,title in zip(axs,['观测导数：所有场景最大相对误差','目标导数：所有场景最大相对误差']):
        ax.set(xlabel='中心差分步长',ylabel='相对误差',title=title);ax.grid(alpha=.2);ax.axhline(1e-5,color='gray',ls='--')
    axs[0].legend(fontsize=8,ncols=2);axs[1].legend(fontsize=8,ncols=2)
    save(fig,'derivative_convergence','中心差分；每条线取四个模型、两个条件的最大值。logQ 的目标导数严格为零，另见 CSV 绝对误差。虚线 1e-5 为预定数值核验标准。')
    # Figure 2: compare rankings on common axes, show all model/noise pairs.
    fig,axs=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    s=per[per.subset=='live']
    for i,f in enumerate(FEATURES):
        z=s[s.feature==f]
        for ax,field in zip(axs,['spearman','top20_energy']):
            ax.scatter(np.full(len(z),i),z[field],s=24,facecolors='none',edgecolors=COLORS[i])
            ax.plot(i,z[field].mean(),marker='_',ms=18,mew=2,color=COLORS[i])
    for ax in axs:
        ax.set_xticks(range(len(FEATURES)),[LABELS[f] for f in FEATURES],rotation=32,ha='right');ax.grid(axis='y',alpha=.2)
    axs[0].set(ylabel='Spearman 相关',title='风险排序与实际窗口 RMS 误差');axs[0].axhline(0,color='gray',lw=.8)
    axs[1].set(ylabel='实际误差能量占比',title='风险最高 20% 窗口捕获的误差');axs[1].axhline(.2,color='gray',ls='--',label='无信息参考 20%');axs[1].legend()
    save(fig,'predictor_comparison','每个空心点是一个已有模型/噪声场景（8 个场景，4 个地质模型）；横线为等权场景平均。正常道，固定 50 ms 窗，未用真值选择窗口。20% 截止并列分数按比例分配。无穷风险保留并列，不以有限伪逆值代替。')
    # Figure 3: rank maps avoid pretending uncertainty and error have same units.
    g=d[(d.model=='complex')&(d.condition=='Moderate')].copy()
    for col in ['U_cal','U_expanded','error_rms']:
        value=g[col].where(g.live.astype(bool)) if col.startswith('U_') else g[col]
        g[col+'_percentile']=value.rank(pct=True)
    fig,axs=plt.subplots(1,3,figsize=(12,5),layout='constrained',sharey=True)
    for ax,col,title in zip(axs,['U_cal','U_expanded','error_rms'],['联合校准风险排序','扩展支持风险排序','实际误差排序']):
        mat=g.pivot(index='window',columns='trace',values=col+'_percentile')
        cmap=plt.colormaps['cividis'].with_extremes(bad='#C0C0C0')
        im=ax.imshow(mat.to_numpy(),origin='upper',aspect='auto',extent=[0,len(mat.columns),1.6,.3],vmin=0,vmax=1,cmap=cmap,interpolation='nearest')
        ids=[1,23,40,61,78,97];ticks=[list(mat.columns).index(i) for i in ids]
        ax.set_xticks(np.asarray(ticks)+.5,[str(i) for i in ids]);ax.set(title=title,xlabel='固定抽取道号（非等间距）')
    axs[0].set_ylabel('双程时间 / s');fig.colorbar(im,ax=axs,label='场景内分位排序，含并列值',shrink=.8)
    save(fig,'risk_error_maps','complex/Moderate；35 条固定抽取道，横轴为样本道序，含额外坏道23、78，非连续等距剖面。风险图灰色表示零观测/继承sigma=realmin，局部风险不适用，封存计算的零占位不能解释为可靠。实际误差图保留这两道。三图分位色标0–1；风险在有效域内排名，误差按全部抽样道排名。无穷风险按最高并列秩显示。原始数值见scored_windows.csv，有效域见prediction_domain.csv。')
    # Figure 4: QScan sign prediction at case level and within-case AUC.
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.6),layout='constrained')
    for i,row in c.iterrows():
        if np.isfinite(row.predicted_log_risk_ratio):
            axs[0].scatter(row.predicted_log_risk_ratio,row.q0_error_ratio,marker='o' if row.condition=='Moderate' else '^',s=50)
            axs[0].annotate(f'{row.model}/{row.condition[0]}',(row.predicted_log_risk_ratio,row.q0_error_ratio),xytext=(3,3),textcoords='offset points',fontsize=8)
    axs[0].axhline(1,color='gray',ls='--');axs[0].axvline(0,color='gray',ls='--')
    axs[0].set(xlabel='预测 log(U_scan / U_base) 的窗中位数',ylabel='原评分 QScan / FullCalRefit 误差比',title='Q0 案例：风险变化能否预测退化')
    xx=np.arange(len(c));axs[1].plot(xx,c.degradation_auc,'o',label='联合风险变化');axs[1].plot(xx,c.qspread_degradation_auc,'s',mfc='none',label='Q 扰动离散度')
    axs[1].axhline(.5,color='gray',ls='--');axs[1].set_xticks(xx,[f'{r.model}/{r.condition[0]}' for _,r in c.iterrows()],rotation=35,ha='right')
    axs[1].set(ylabel='窗口退化识别 AUROC',ylim=(0,1),title='实际 MSE 增加作为退化标签');axs[1].legend()
    save(fig,'qscan_failure_prediction','Q0 是旧实验施加 Q 扰动为零。只有 primary 为给定真Q；complex为独立标定Q，F1/F2为上尺度宏观Q。横轴只在两种风险都有限时定义；覆盖率逐场景列在 qscan_diagnosis.csv，缺失不填零。每个噪声条件使用同一地质模型。')
    # Figure 5: actual held-geology prediction, no window-random split.
    fig,axs=plt.subplots(1,2,figsize=(10.7,4.2),layout='constrained')
    marks=['o','s','^','v','D','x']
    for i,(variant,g) in enumerate(cv.groupby('variant',sort=False)):
        for ax,col in zip(axs,['spearman','r2']):ax.plot(g.held_model,g[col],marker=marks[i],label=variant)
    axs[0].set(ylabel='留出模型的排序相关',title='额外预测能力');axs[1].set(ylabel='留出模型的秩回归 R²',title='同一地质的两种噪声条件同时留出');axs[1].legend(fontsize=7,ncols=2)
    for ax in axs:ax.grid(alpha=.2);ax.axhline(0,color='gray',lw=.6)
    save(fig,'held_model_prediction','固定线性秩回归；所有输入先在场景内转为分位秩，实际误差秩只在评分阶段计算。基础输入：双程时间、观测SNR倒数、Q扰动离散度。新增U_cal或U_expanded，不搜索参数。4个留出地质模型；属于已有数据诊断，不能作为最终新模型验证。')
    # Direct correspondence: pre-scored risk rank deciles, empirical MSE.
    bins=[]
    for f in ['q_spread','U_r','U_cal']:
        for scene,g in normal.groupby('scene'):
            rank=g[f].rank(pct=True);decile=np.minimum((rank*10).astype(int),9)
            for b in range(10):
                m=decile==b
                bins.append(dict(feature=f,scene=scene,decile=b+1,n=int(m.sum()),
                                 relative_mse=float(g.loc[m,'base_mse'].mean()/g.base_mse.mean())))
    bins=pd.DataFrame(bins);bins.to_csv(ROOT/'analysis/risk_error_deciles.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained',sharey=True)
    for ax,f in zip(axs,['q_spread','U_r','U_cal']):
        z=bins[bins.feature==f];agg=z.groupby('decile').relative_mse.agg(['mean','min','max'])
        for scene,g in z.groupby('scene'):ax.scatter(g.decile,g.relative_mse,s=12,color='#888888',alpha=.6)
        ax.plot(agg.index,agg['mean'],color='#0072B2',marker='o',label='8 场景等权均值')
        ax.axhline(1,color='gray',ls='--');ax.set(xlabel='预测风险十分位（低 → 高）',title=LABELS[f]);ax.grid(alpha=.2)
    axs[0].set_ylabel('真实窗口 MSE / 该场景平均 MSE');axs[-1].legend(fontsize=8)
    save(fig,'risk_error_correspondence','风险分位在真值评分前固定；真值只用于纵轴。灰点为每场景每风险十分位的平均真实MSE，除以该场景全部正常道窗口的平均MSE；蓝线为8场景等权均值。点的差异为已有场景间差异，不是置信区间。未作跨窗独立性推断；全部数值见risk_error_deciles.csv。')
    result=dict(derivative_observation_max_at_1e4=float(hcheck.observation_relative.max()),
                derivative_target_max_at_1e4=float(hcheck.target_relative.max()),
                frozen_base_max=float(audit.frozen_base_difference.max()),frozen_scan_max=float(audit.frozen_scan_difference.max()),
                old_score_max_difference=max(x['absolute_difference'] for x in checks),
                macro_live=macro[macro.subset=='live'].to_dict(orient='records'),
                leave_model_out=cv.groupby('variant')[['spearman','r2']].mean().to_dict(orient='index'),
                cases=c.to_dict(orient='records'),rows=len(d),live_rows=len(normal),geological_models=4)
    (ROOT/'analysis/summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=True),encoding='utf-8')
    print(json.dumps(result,indent=2,ensure_ascii=True))

if __name__=='__main__':main()
