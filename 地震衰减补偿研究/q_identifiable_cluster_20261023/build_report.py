from pathlib import Path
import json,csv,hashlib,datetime,sys,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

ROOT=Path(__file__).resolve().parent
def readrows(name):
    with (ROOT/'analysis'/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    spec=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'))
    dev=json.loads((ROOT/'analysis/development_summary.json').read_text(encoding='utf-8'))
    val=json.loads((ROOT/'analysis/validation_summary.json').read_text(encoding='utf-8'))
    sens=json.loads((ROOT/'analysis/sensitivity_summary.json').read_text(encoding='utf-8'))
    methods=spec['methods'];vrows=readrows('validation_metrics.csv');mrows=readrows('validation_modes.csv');srows=readrows('sensitivity_metrics.csv')
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    colors=['#666666','#999933','#0072B2','#56B4E9','#009E73','#D55E00'];markers=['s','v','^','D','P','o']
    cases=[v['name'] for v in spec['variants']];conditions=['Moderate','Strong'];pairs=[(c,k) for c in cases for k in conditions];x=np.arange(len(pairs))
    fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True,layout='constrained')
    for ax,key,title in zip(axs,['LiveTraceError','LowEnergyRegionError'],['Waveform reconstruction error','Low-energy region waveform error']):
        for m,col,mk in zip(methods,colors,markers):
            y=[float(next(r[key] for r in vrows if r['Case']==c and r['Condition']==k and r['Method']==m)) for c,k in pairs]
            ax.plot(x,y,marker=mk,color=col,lw=1.2,label=m)
        ax.set(ylabel='Relative error (lower is better)',title=title,ylim=(0,None));ax.grid(axis='y',alpha=.2)
    axs[0].legend(ncol=3,fontsize=8)
    axs[-1].set_xticks(x,[f'{c}\n{k}' for c,k in pairs],fontsize=8)
    fig.suptitle('Frozen validation: three new same-family geological configurations')
    fig.supxlabel('22 preselected traces per condition. Conditions/traces are dependent; no confidence intervals are implied.',fontsize=9)
    fig.savefig(out/'validation_errors.png',dpi=170);fig.savefig(out/'validation_errors.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for m,col,mk in zip(methods,colors,markers):
        r=sorted([r for r in srows if r['Method']==m],key=lambda r:float(r['QPointScale']))
        axs[0].plot([float(a['QPointScale']) for a in r],[float(a['LiveTraceError']) for a in r],marker=mk,color=col,label=m)
        axs[1].plot([float(a['QPointScale']) for a in r],[float(a['OutputChangeFromNominal']) for a in r],marker=mk,color=col,label=m)
    axs[0].set(ylabel='Relative waveform error',title='Recovery under point-Q perturbation',ylim=(0,None))
    axs[1].set(ylabel='Norm change relative to nominal output',title='Output sensitivity (not an accuracy score)',ylim=(0,None))
    for ax in axs:ax.set_xlabel('Multiplier on estimated layer Q');ax.set_xticks([.8,1,1.25]);ax.grid(alpha=.2)
    axs[0].legend(fontsize=8)
    fig.suptitle('One observation, artificial perturbations: not additional independent validation')
    fig.savefig(out/'Q_sensitivity.png',dpi=170);fig.savefig(out/'Q_sensitivity.pdf');plt.close(fig)
    run=Path(val['run']);case=cases[0];condition='Moderate'
    truth=np.load(ROOT/'truth'/f'{case}.npz');t=truth['t'];xx=truth['x']/1000;target=truth[condition+'_target']
    result=json.loads((run/case/(condition+'.json')).read_text(encoding='utf-8'))
    arrays=[target]+[np.array(result[m]) for m in ['FullCalRefit','PrunedFullCalRefit','QAwareCluster']]
    limit=max(float(np.max(np.abs(a))) for a in arrays)
    fig,axs=plt.subplots(1,4,figsize=(13,6),sharey=True,layout='constrained')
    for ax,a,title in zip(axs,arrays,['Clean target','FullCalRefit','Group pruning','QAwareCluster']):
        im=ax.pcolormesh(xx,t,a,shading='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
        ax.set(title=title,xlabel='Position (km)',ylim=(t[-1],t[0]))
    axs[0].set_ylabel('Time (s)');fig.colorbar(im,ax=axs,label='Amplitude: same limits, no trace normalization',shrink=.7)
    fig.suptitle(f'Fixed example: {case} / {condition} (first registered validation case)')
    fig.supxlabel('22 preselected original traces shown at their actual spatial coordinates; no missing traces are reconstructed for display.',fontsize=9)
    fig.savefig(out/'fixed_example_gathers.png',dpi=170);fig.savefig(out/'fixed_example_gathers.pdf');plt.close(fig)
    # First and middle registered traces, never chosen using reconstruction error.
    fig,axs=plt.subplots(2,1,figsize=(12,6.5),sharex=True,layout='constrained')
    for ax,j in zip(axs,[0,len(spec['indices_zero_based'])//2]):
        ax.plot(t,target[:,j],color='black',lw=1.1,label='Clean target')
        for m,col in [('FullCalRefit','#666666'),('PrunedFullCalRefit','#0072B2'),('QAwareCluster','#D55E00')]:
            ax.plot(t,np.array(result[m])[:,j],color=col,lw=.9,label=m)
        ax.set(ylabel='Amplitude',title=f'Original trace {spec["indices_zero_based"][j]+1}',xlim=(.25,1.65));ax.grid(alpha=.2)
    axs[0].legend(ncol=4,fontsize=8);axs[-1].set_xlabel('Time (s)')
    fig.suptitle('Fixed traces: amplitudes unchanged; no independent per-trace display gains')
    fig.savefig(out/'fixed_example_waveforms.png',dpi=170);fig.savefig(out/'fixed_example_waveforms.pdf');plt.close(fig)
    plots=[]
    for p in out.glob('*.png'):
        with Image.open(p) as im:assert im.width>1000 and im.height>500;dims=list(im.size)
        plots.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=digest(p),pixels=dims))
    (out/'manifest.json').write_text(json.dumps(dict(plots=plots,python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,error_bars='none: no independent replicates for each configuration',gather_amplitude_limit=limit,gather_clipping='none',example_selection='first registered validation geology/condition; traces first and middle registered indices',source_csvs={p.name:digest(p) for p in (ROOT/'analysis').glob('*metrics.csv')}),indent=2),encoding='utf-8')
    def improve(a,b):return 100*(1-a/b)
    v=val['means'];d=dev['means'];candidate=v['QAwareCluster'];base=v['FullCalRefit'];passed=val['gates']['all_pass']
    normal_gain=improve(candidate['LiveTraceError'],base['LiveTraceError']);weak_gain=improve(candidate['LowEnergyRegionError'],base['LowEnergyRegionError'])
    lines=['# Q 不确定性感知反射簇实验结果','','## Material Passport','', '- Origin Skill: academic-research-suite / experiment-agent','- Origin Mode: validate',f'- Origin Date: {datetime.date.today().isoformat()}','- Verification Status: ANALYZED; actual MATLAB execution and numerical limiting controls passed; no independent reproduction', '- Version Label: q_cluster_v1_completed','','## 结论','',
        f'首轮预登记工程验收：**{"通过" if passed else "未通过"}**。新模型等权正常道误差相对 FullCalRefit 改善 {normal_gain:+.2f}%，低能量区误差改善 {weak_gain:+.2f}%。负数表示退化。',
        '这是新机制的实际首轮检验，不是通过重命名 FullCalRefit 得到的方法。即使某一项指标改善，也不能跳过整体、局部与消融门槛。未证明首创性。',
        '', '## 方法与创新假设','',
        '保持原独立校准和 L1 支持发现，使用校准 logQ 的局部误差尺度构造固定先导模型误差白化；对每个反射簇分解组合模式，再按固定复杂度目标降低簇秩并联合重估振幅。目标是保留可辨识的组合波形，而非整簇删除。',
        '创新候选集中在 Q 误差驱动的局部自由度选择。SVD、固定先导误差协方差、自适应参数化和稀疏反演本身均不是本次首创。完整规则见 PLAN_ZH.md，全部候选源码、开发版本快照和冻结哈希保留。',
        '', '## 实际执行范围','',
        '开发为原复杂模型×2条件×22道=44条波形；验证为3个新地质配置×2条件×22道=132条波形；Q点估计扰动另44条波形。共220条观测/扰动处理条件，每条比较6种恢复方法，共1320条方法波形；重复扰动不是新观测或独立样本。数值控制另计。',
        '三个新配置分别改变层位密度、Q深度趋势、主频、相位、偏移距、地质随机种子、校准噪声与观测噪声。依然继承同族准1D属性驱动正演，不是完整二维粘弹性正演或实测验证。',
        '预先从99道中均匀选22道。本次选中位置恰好不含原生成器的第23/78坏道，因此全道和正常道评分相同，不能从本轮推断坏道鲁棒性。',
        '', '## 冻结验证均值（每模型、每条件等权）','', '| 方法 | 波形误差 | 低能量区误差 | 深部误差 | 相关系数 |','|---|---:|---:|---:|---:|']
    for m in methods:
        a=v[m];lines.append(f"| {m} | {a['LiveTraceError']:.6f} | {a['LowEnergyRegionError']:.6f} | {a['DeepLiveError']:.6f} | {a['CorrelationLive']:.6f} |")
    lines+=['','## 逐模型核对与校准诊断','','| 模型 | FullCalRefit误差 | 候选误差 | 候选相对改善 | Q点估计触边 |','|---|---:|---:|---:|---|']
    for case in cases:
        av=lambda m:float(np.mean([float(r['LiveTraceError']) for r in vrows if r['Case']==case and r['Method']==m]))
        bc=av('FullCalRefit');cc=av('QAwareCluster');cal=json.loads((run/case/'calibration.json').read_text(encoding='utf-8'))
        lines.append(f"| {case} | {bc:.6f} | {cc:.6f} | {improve(cc,bc):+.2f}% | {cal['estimate'].get('q_boundary',False)} |")
    lines+=['','其中 Q 触边记录仍纳入全部评分；没有按校准或处理效果删除不利模型。触边提示校准受限，不等于候选已自动修正偏差。']
    lines+=['','## 配对消融','']
    for other in ['PrunedFullCalRefit','GlobalTSVD','ClusterModesNoQ','QAwareRidge']:
        lines.append(f"- QAwareCluster 相对 {other} 的波形误差改善：{improve(candidate['LiveTraceError'],v[other]['LiveTraceError']):+.2f}%。")
    lines+=['','Q误差机制是否有收益，应同时查看 QAwareRidge 对 FullCalRefit，以及 QAwareCluster 对 ClusterModesNoQ；簇秩机制应查看 QAwareCluster 对 QAwareRidge。不能把组合结果都归功于 Q。',
        '', '## 预登记门槛','', '| 门槛 | 结果 |','|---|---|']
    for key,value in val['gates'].items():lines.append(f'| {key} | {"PASS" if value else "FAIL"} |')
    nwins=sum(float(next(r['LiveTraceError'] for r in vrows if r['Case']==c and r['Condition']==k and r['Method']=='QAwareCluster'))<float(next(r['LiveTraceError'] for r in vrows if r['Case']==c and r['Condition']==k and r['Method']=='FullCalRefit')) for c,k in pairs)
    lines+=['',f'候选在6个模型×条件组合中的整体误差胜出 {nwins}/6；条件来自3个模型，不把该胜率解释成独立成功概率。']
    lines+=['','## Q点估计扰动','', '第一新模型的同一 Moderate 观测，名义 Q 与0.8/1.25倍估计 Q；使用相同的名义校准协方差来源，不使用真实 Q。以下最差误差与输出变化分别衡量准确性和敏感性，稳定不等于正确。',
        '', '| 方法 | 三点最差波形误差 | 最大输出变化 |','|---|---:|---:|']
    for m in methods:
        a=sens['summary'][m];lines.append(f"| {m} | {a['worst_error']:.6f} | {a['maximum_output_change']:.6f} |")
    allchecks=val['checks']+sens['checks'];maxnr=max(c['max_normal'] for c in allchecks)
    lines+=['','## 数值、流程与解释边界','',
        f'- 开发、验证和敏感性 MATLAB 运行均实际完成；验证/敏感性最大正规方程相对残差 {maxnr:.3e}。全部数值检查详见各运行 checks.json。',
        '- 校准协方差为带固定弱先验和尺度截断的局部线性化估计，不是已验证的贝叶斯后验/置信区间；第四层可能缺乏校准信息。',
        '- 使用单个固定先导计算 Q 模型误差，忽略更高阶项、源子波空间变化和有色噪声协方差；不能称为完整模型误差似然。',
        '- 模式方向、支持和误差权重依赖观测先导，复杂度项不是经校准的显著性检验。更小拟合目标并不保证真值误差更低。',
        '- 普通截断、NoQ消融和候选共用已有 L1 支持。候选无法恢复支持发现阶段已遗漏且未补入的反射。',
        '- validation_metrics.csv 中 ObservedReplayResidual 字段实际为 ||输出-输入||/||输入||，仅表示处理变化量，不能解释为正演重放残差或去噪改善。本报告不将其用于晋级。',
        '- 模式路径中的数据拟合项可由 objective_path 减去 sigma² log(nt)×rank 得到；白化不同的方法其数值不能直接当作相同似然比较。',
        '- 首轮开发已显示整体退化，仍按已登记规则完成验证以检验取舍；未根据新模型结果调参。开发快照之后仅增加缺字段防护、数值控制与预登记Q扰动入口，未改变已有开发条件下估计器。',
        '- 本次冻结验证的三个地质配置现已被查看。若利用本轮结果修改算法，这些配置应转为开发/回归资料，不能继续充当后续独立确认集。',
        '- 11/11类统计谬误检查覆盖见 analysis/validation_RESULTS_ZH.md；不做显著性、独立样本置信区间或跨地区泛化声明。',
        '', '## 后续定位（尚未实验）','',
        '当前模式选择目标控制的是加权观测拟合与自由度，没有直接保护补偿输出 W r 的振幅。应先在开发资料上量化被删模式对输出波形的贡献，再研究带输出保真约束的模式选择。该解释是下一步假设；本轮不能证明某一个阈值或物理因素是唯一失败原因。',
        '', '## 可复现入口','',
        '依赖 MATLAB R2024b（未使用额外工具箱）；数据生成/评分/绘图使用 Python 3.12、NumPy、Matplotlib、Pillow。各运行 SOURCE_MANIFEST.json 保存准确源码快照，FROZEN_CANDIDATE.json 保存验证前源码、输入和真值哈希。',
        '重新执行已冻结输入时，运行 python run_experiment.py development 或 python run_experiment.py validation，会创建新的运行目录并保留历史。运行后 python analyze.py development / validation 生成对应评分。新一轮重复实验不要覆盖历史真值或再次运行有防覆盖保护的生成器。',
        '先运行 prepare_sensitivity.py，再运行 run_experiment.py sensitivity 和 summarize_sensitivity.py 可复核已登记 Q 扰动；本次输入已存在，重复时直接使用现有输入运行 sensitivity，不再准备覆盖。最后运行 build_report.py 生成本报告和图件。',
        '', '## 文件','', '- RESULTS_ZH.md：本报告。','- analysis/validation_metrics.csv：逐模型、逐条件、逐方法指标。','- analysis/validation_trace_metrics.csv：全部逐道误差。','- analysis/validation_modes.csv：支持数、保留秩、复杂度目标和数值残差。','- analysis/sensitivity_metrics.csv：所有 Q 扰动指标。','- figures/：全部PNG/PDF、来源与显示尺度清单。','- runs/：原始MATLAB结果、校准诊断、日志、PID、退出码和源码快照。']
    (ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
    index='<html lang="zh-CN"><meta charset="utf-8"><title>Q反射簇实验</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px}img{width:100%;height:auto}a{margin-right:20px}</style><h1>Q 不确定性感知反射簇实验</h1>'+f'<p>首轮验收：{"通过" if passed else "未通过"}。波形误差相对 FullCalRefit 改善 {normal_gain:+.2f}%；低能量区改善 {weak_gain:+.2f}%。完整不利结果保留。</p><p><a href="RESULTS_ZH.md">完整报告</a><a href="analysis/validation_metrics.csv">验证指标</a><a href="analysis/sensitivity_metrics.csv">Q扰动指标</a></p>'
    for p in plots:index+=f'<h2>{p["path"]}</h2><img src="{p["path"]}" alt="实验原始指标或固定选择波形，完整数值见CSV">'
    (ROOT/'index.html').write_text(index+'</html>',encoding='utf-8')
    print(json.dumps(dict(passed=passed,normal_gain_percent=normal_gain,low_energy_gain_percent=weak_gain,report=str(ROOT/'RESULTS_ZH.md'),figures=plots),indent=2),flush=True)

if __name__=='__main__':main()
