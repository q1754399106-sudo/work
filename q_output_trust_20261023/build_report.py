from pathlib import Path
import json,csv,hashlib,datetime,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parent
def rows(name):
 with (ROOT/'analysis'/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
 spec=json.loads((ROOT/'SPEC.json').read_text());dev=json.loads((ROOT/'analysis/development_summary.json').read_text());val=json.loads((ROOT/'analysis/validation_summary.json').read_text());sens=json.loads((ROOT/'analysis/sensitivity_summary.json').read_text());selected=json.loads((ROOT/'SELECTED.json').read_text())
 vr=rows('validation_metrics.csv');mr=rows('validation_modes.csv');sr=rows('sensitivity_metrics.csv');methods=spec['methods'];cases=[v['name'] for v in spec['variants']];pairs=[(c,k) for c in cases for k in ['Moderate','Strong']]
 means=val['means'];candidate='QAwareOutputTrust';gain=lambda a,b:100*(1-a/b);overall=gain(means[candidate]['LiveTraceError'],means['FullCalRefit']['LiveTraceError']);low=gain(means[candidate]['LowEnergyRegionError'],means['FullCalRefit']['LowEnergyRegionError'])
 out=ROOT/'figures';plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 colors=['#666666','#999933','#009E73','#56B4E9','#0072B2','#D55E00'];markers=['s','v','P','D','^','o']
 fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True,layout='constrained')
 for ax,key,title in zip(axs,['LiveTraceError','LowEnergyRegionError'],['Waveform reconstruction error','Low-energy region waveform error']):
  for m,col,mk in zip(methods,colors,markers):ax.plot(range(6),[float(next(r[key] for r in vr if r['Case']==c and r['Condition']==k and r['Method']==m)) for c,k in pairs],color=col,marker=mk,lw=1.2,label=m)
  ax.set(ylabel='Relative error (lower is better)',title=title,ylim=(0,None));ax.grid(axis='y',alpha=.2)
 axs[0].legend(ncol=3,fontsize=8);axs[-1].set_xticks(range(6),[f'{c}\n{k}' for c,k in pairs],fontsize=8);fig.suptitle('Round 2: frozen validation on three new same-family geologies');fig.supxlabel('22 fixed traces per condition; conditions/traces are dependent. No confidence intervals implied.',fontsize=9)
 fig.savefig(out/'validation_errors.png',dpi=170);fig.savefig(out/'validation_errors.pdf');plt.close(fig)
 new=[r for r in mr if r['Method']==candidate];fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
 for i,(c,k) in enumerate(pairs):
  sub=[r for r in new if r['Case']==c and r['Condition']==k];axs[0].plot(range(22),[float(r['MaxProtectedChangeSD'])/selected['kappa'] for r in sub],marker=markers[i],markersize=3,color=colors[i],lw=.8,label=f'{c}/{k}');axs[1].scatter([float(r['RankBefore']) for r in sub],[float(r['RankAfter']) for r in sub],s=20,marker=markers[i],color=colors[i],label=f'{c}/{k}')
 axs[0].axhline(1,color='black',ls='--',lw=1);axs[0].set(xlabel='Registered spatial sample index',ylabel='Cumulative protected change / budget',title='Output constraint use (not truth accuracy)',ylim=(0,1.08));axs[0].legend(fontsize=6)
 lim=max(float(r['RankBefore']) for r in new)*1.05;axs[1].plot([0,lim],[0,lim],'k--',lw=1);axs[1].set(xlabel='Initial support size',ylabel='Retained modes',title='Actual mode compression',xlim=(0,lim),ylim=(0,lim));axs[1].grid(alpha=.2)
 fig.suptitle(f'Output trust constraints, selected kappa={selected["kappa"]:g}');fig.savefig(out/'constraint_and_rank.png',dpi=170);fig.savefig(out/'constraint_and_rank.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 ks=[.5,1,2];names=['Trust05','Trust10','Trust20']
 for ax,key in zip(axs,['LiveTraceError','LowEnergyRegionError']):
  vals=[dev['means'][m][key] for m in names];upper=1.08*max(vals+[dev['means']['FullCalRefit'][key],dev['means']['QAwareRidge'][key]])
  label='Waveform relative error' if key=='LiveTraceError' else 'Low-energy region relative error'
  ax.plot(ks,vals,'o-',color='#D55E00',label='Output trust candidates');ax.axhline(dev['means']['FullCalRefit'][key],color='#666666',ls='--',label='FullCalRefit');ax.axhline(dev['means']['QAwareRidge'][key],color='#009E73',ls=':',label='QAwareRidge');ax.set(xlabel='kappa (development only)',ylabel=label,ylim=(0,upper));ax.set_xticks(ks);ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.suptitle('All preregistered development candidates retained');fig.savefig(out/'development_selection.png',dpi=170);fig.savefig(out/'development_selection.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
 for m,col,mk in zip(methods,colors,markers):
  ss=sorted([r for r in sr if r['Method']==m],key=lambda r:float(r['QPointScale']));x=[float(r['QPointScale']) for r in ss]
  axs[0].plot(x,[float(r['LiveTraceError']) for r in ss],color=col,marker=mk,label=m);axs[1].plot(x,[float(r['OutputChange']) for r in ss],color=col,marker=mk,label=m)
 axs[0].set(ylabel='Relative waveform error',title='Estimated-Q perturbation',ylim=(0,None));axs[1].set(ylabel='Norm change from nominal output',title='Sensitivity, not an accuracy score',ylim=(0,None))
 for ax in axs:ax.set_xlabel('Multiplier on estimated layer Q');ax.set_xticks([.8,1,1.25]);ax.grid(alpha=.2)
 axs[0].legend(fontsize=7);fig.suptitle('Same observation, three Q-point settings; not independent validation cases');fig.savefig(out/'Q_sensitivity.png',dpi=170);fig.savefig(out/'Q_sensitivity.pdf');plt.close(fig)
 run=Path(val['run']);case=cases[0];truth=np.load(ROOT/'truth'/f'{case}.npz');t=truth['t'];xx=truth['x']/1000;target=truth['Moderate_target'];result=json.loads((run/case/'Moderate.json').read_text());arrays=[target]+[np.array(result[m]) for m in ['FullCalRefit','QAwareCluster',candidate]];limit=max(float(np.abs(a).max()) for a in arrays)
 fig,axs=plt.subplots(1,4,figsize=(13,6),sharey=True,layout='constrained')
 for ax,a,title in zip(axs,arrays,['Clean target','FullCalRefit','Previous cluster','Output trust']):
  im=ax.pcolormesh(xx,t,a,shading='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True);ax.set(title=title,xlabel='Position (km)',ylim=(t[-1],t[0]))
 axs[0].set_ylabel('Time (s)');fig.colorbar(im,ax=axs,label='Shared amplitude scale; no trace normalization',shrink=.7);fig.suptitle(f'Fixed example: first registered model {case}, Moderate');fig.supxlabel('Actual 22 sample positions; same amplitude scale, no display gain or clipping.',fontsize=9);fig.savefig(out/'fixed_gathers.png',dpi=170);fig.savefig(out/'fixed_gathers.pdf');plt.close(fig)
 manifest=[]
 for p in sorted(out.glob('*.png')):
  with Image.open(p) as im:dims=list(im.size);assert min(dims)>500
  manifest.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),pixels=dims))
 (out/'manifest.json').write_text(json.dumps(dict(plots=manifest,python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,gather_shared_limit=limit,clipping='none',selection='first registered validation model/Moderate, never chosen by error',uncertainty='no confidence intervals; dependent conditions',source_csv_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'analysis').glob('*.csv')}),indent=2))
 changed=sum(float(r['Removed'])>0 for r in new)/len(new);prot=np.array([float(r['ProtectedModes']) for r in new]);refchange=np.array([float(r['ReferenceOutputChange']) for r in new]);wins=sum(float(next(r['LiveTraceError'] for r in vr if r['Case']==c and r['Condition']==k and r['Method']==candidate))<float(next(r['LiveTraceError'] for r in vr if r['Case']==c and r['Condition']==k and r['Method']=='FullCalRefit')) for c,k in pairs)
 lines=['# 第二轮：输出保幅信赖域实验结果','','## Material Passport','', '- Origin Skill: academic-research-suite / experiment-agent','- Origin Mode: validate',f'- Origin Date: {datetime.date.today().isoformat()}','- Verification Status: ANALYZED; actual MATLAB runs and numerical controls completed; no independent reproduction','- Version Label: q_output_trust_v2','','## 结论','',
 f'预登记验收：**{"通过" if val["gates"]["all_pass"] else "未通过"}**。正常道波形误差相对FullCalRefit改善 {overall:+.2f}%，低能量区误差改善 {low:+.2f}%。负值表示退化。候选在6个模型×条件组合胜出 {wins}/6；这些条件只来自3个地质配置。',
 f'开发按既定规则从0.5/1/2三个kappa中选定 **{selected["kappa"]:g}**，新验证与Q扰动全程固定。实际发生模式删除的验证道占 {changed:.1%}；平均受保护输出方向数 {prot.mean():.2f}；相对完整Q加权参考的输出变化均值 {refchange.mean():.4f}。',
 '', '## 机制与先前诊断','',
 '先前4个已查看模型的8个条件中，簇秩删除相对QAwareRidge均使真值均方误差增加。用delta=S_cluster-S_ref核对误差变化恒等式 2<error_ref,delta>+||delta||²，完整逐条件结果见analysis/previous_round_output_diagnosis.csv；仅作开发诊断。',
 '本轮保持Q误差白化与簇内模式表示，但每次删除还需满足输出空间的累计信赖域：在参考投影超过2倍近似标准差的方向，整个候选相对完整参考的变化不能超过kappa倍该标准差。协方差通过支持到补偿波形的映射传播，约束相对完整参考，防止多步漂移。数学定义见PLAN_ZH.md和core/qtrust_refit.m。',
 '参考可能有偏；支持选择、噪声假设、校准协方差、固定先导和尺度截断均为近似，因此不是经覆盖率验证的置信区间，也不保证真实保幅。自适应参数化、SVD和误差传播各自已有研究；本实现的首创性尚未确认。',
 '', '## 冻结验证均值','','| 方法 | 波形误差 | 低能量区误差 | 深部误差 | 相关系数 |','|---|---:|---:|---:|---:|']
 for m in methods:
  a=means[m];lines.append(f"| {m} | {a['LiveTraceError']:.6f} | {a['LowEnergyRegionError']:.6f} | {a['DeepLiveError']:.6f} | {a['CorrelationLive']:.6f} |")
 lines+=['','## 逐模型与配对消融','','| 模型 | FullCalRefit误差 | 候选误差 | 相对改善 | 校准Q触边 |','|---|---:|---:|---:|---|']
 for case in cases:
  av=lambda m:float(np.mean([float(r['LiveTraceError']) for r in vr if r['Case']==case and r['Method']==m]));bc=av('FullCalRefit');cc=av(candidate);cal=json.loads((run/case/'calibration.json').read_text());lines.append(f"| {case} | {bc:.6f} | {cc:.6f} | {gain(cc,bc):+.2f}% | {cal['estimate'].get('q_boundary',False)} |")
 lines+=['','本轮全部模型均纳入评分。','']
 for m in ['QAwareCluster','QAwareRidge','OutputTrustNoQ']:lines.append(f"- 相对 {m} 的波形误差改善：{gain(means[candidate]['LiveTraceError'],means[m]['LiveTraceError']):+.2f}%。")
 lines+=['','## 预登记门槛','','| 门槛 | 结果 |','|---|---|']
 for k,v in val['gates'].items():lines.append(f'| {k} | {"PASS" if v else "FAIL"} |')
 lines+=['','## Q点估计扰动','','第一新模型的同一Moderate观测，Q点估计乘0.8/1/1.25；保持名义校准协方差来源，不使用真Q。','','| 方法 | 三点最差误差 | 最大输出变化 |','|---|---:|---:|']
 for m in methods:
  a=sens['summary'][m];lines.append(f"| {m} | {a['worst_error']:.6f} | {a['max_output_change']:.6f} |")
 covcheck=json.loads((ROOT/'analysis/covariance_dense_symmetric.json').read_text())
 lines+=['','## 补充条件协方差核验','',
 f'固定支持、固定先导、已知高斯噪声与模型误差的独立稠密计算，512次模拟方差相对解析方差的偏差中位数为 {100*covcheck["median_relative_variance_error"]:.2f}%，最大为 {100*covcheck["maximum_relative_variance_error"]:.2f}%，通过原定数值要求。该检查验证条件线性误差传播，不验证真实Q误差分布或选择后的置信覆盖率。',
 '补充检验首次失败：稠密协方差约9e-17的非对称舍入误差使重根特征分解进入一般非对称路径，得到非正交特征基；第二次增加诊断输出时暴露微小复数导致JSON写入失败。随后仅在补充检验中显式对称化协方差，保留原样本数和判据，通过核验。主估计器使用薄SVD白化，未改动。失败日志与空诊断文件均保留，详见analysis/covariance*。',
 '', '## 范围、核查与后续边界','',
 '- 开发44条波形条件×7方法；新模型验证132条×6方法；额外Q扰动44条×6方法。共220条观测/扰动条件，1364条方法波形；其中真实不同的输入观测只有176条。',
 '- 三个新地质配置与上一轮均不同，改变层密度、Q趋势、主频、相位、偏移距以及校准/观测种子；仍为同族准1D属性驱动正演，不能称为跨正演或实测验证。',
 '- 固定22道空间抽样不含原第23/78坏道；不据此推断坏道鲁棒性。低能量区分位定义只在真值评分程序中使用，不解释为单个弱反射系数误差。',
 '- 初始支持及残差补入沿用FullCalRefit；本轮不解决已遗漏事件、任意子波、未知层位或不可靠Q校准。',
 '- 所有基线共用观测和校准；真值从求解输入中隔离。原MATLAB依赖哈希保持一致，开发/验证/扰动均保存实际执行源码和输入快照。',
        '- 在独立MATLAB进程中从准备缓存重放一个选定开发道，输出相对差为0；该检查只核对缓存和估计器一致性。',
 '- 在开发快照后只新增了预登记Q扰动入口；名义开发估计器qtrust_refit.m未改变。验证前完整冻结，未根据新模型结果改参数。',
 '- 完整记录正规方程残差、保护方向数、累计预算使用量、删除秩、被保护规则拒绝的提案和所有开发候选。数值通过与科学效果晋级分开判断。',
 '- 11/11类统计谬误按本描述性范围核查：分模型展示避免聚合掩盖；不从区域均值推断单事件；披露生成模型选择与先导依赖；不按结果删案例或控制后验筛选变量；基率不适用；不把开发变化当重复测量收益；保留失败；多候选仅在开发选择、不作显著性；完整冻结分析分支；归因限配对机制；受控计算不涉及观察性反向因果。',
 '- 若继续利用本轮结果改算法，本轮三个新配置应转为开发/回归资料，后续需要新的确认资料。',
 '', '## 复现与文件','',
 'MATLAB R2024b；Python 3.12 + NumPy/Matplotlib/Pillow。run_round2.py development/validation/sensitivity 各自新建运行目录；analyze_round2.py development/validation 评分；prepare_sensitivity.py（输入已有时跳过）和summarize_sensitivity.py处理扰动；build_report.py生成报告图件；audit_completion.py核查交付。',
 'RESULTS_ZH.md为本报告，SELECTED.json保留开发选择，FROZEN_CANDIDATE.json保留验证前源码/输入/真值哈希。analysis/含全部逐道和逐条件CSV。runs/保留日志、退出码、精确快照、原始输出以及不含真值的prepared.mat缓存。图件统一尺度、不逐道归一化、不截幅。']
 (ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
 html='<html lang="zh-CN"><meta charset="utf-8"><title>输出保幅信赖域实验</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px}img{width:100%}</style><h1>第二轮：输出保幅信赖域</h1>'+f'<p>验收：{"通过" if val["gates"]["all_pass"] else "未通过"}；波形误差改善 {overall:+.2f}%，低能量区 {low:+.2f}%。</p><p><a href="RESULTS_ZH.md">完整报告</a> · <a href="analysis/validation_metrics.csv">全部指标</a></p>'
 for p in manifest:html+=f'<h2>{p["path"]}</h2><img src="{p["path"]}" alt="完整配对实验指标、数值约束与固定选择的道集；数值见CSV">'
 (ROOT/'index.html').write_text(html+'</html>',encoding='utf-8');print(json.dumps(dict(passed=val['gates']['all_pass'],overall_gain=overall,low_energy_gain=low,kappa=selected['kappa'],changed_trace_fraction=changed,report=str(ROOT/'RESULTS_ZH.md')),indent=2))
if __name__=='__main__':main()
