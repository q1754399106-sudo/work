"""All paired cases: common-scale seismic gathers, wiggles, traces and report."""
from pathlib import Path
import csv,json,hashlib,statistics,math,html,platform,shutil,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def avg(rs,key):
 x=[float(r[key]) for r in rs if math.isfinite(float(r[key]))]
 return statistics.mean(x) if x else float('nan')
def writecsv(p,rows):
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
SCENES={'SmoothAVO':'连续AVO','CoherentNoiseAVO':'相干噪声AVO','FaultSection':'断层','WeakSection':'弱反射','PinchoutSection':'尖灭'}
def main():
 run=Path(read(ROOT/'LATEST_RUN.json')['run']);assert read(run/'execution.json')['status']=='COMPLETED'
 checks=read(run/'matlab/checks.json');assert checks['completed'] and checks['observations']==60 and checks['rows']==120
 manifest=read(run/'SOURCE_MANIFEST.json');assert all(sha(run/'source/core_p1'/k)==v for k,v in manifest.items())
 src=run/'source/core_p1';helpers={}
 for name in ['ray_absorption.m','attenuation_matrix.m','ridge_operator.m']:
  assert (src/'matlab_src'/name).read_bytes()==(src/'original_v1'/name).read_bytes();helpers[name]=sha(src/'matlab_src'/name)
 dest=ROOT/'figures';dest.mkdir(exist_ok=True);dataout=ROOT/'analysis';dataout.mkdir(exist_ok=True)
 rows=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')))
 keys=['Model','Scenario','SNR_dB','Seed','Method'];assert len({tuple(r[k] for k in keys) for r in rows})==120
 bycase={}
 for r in rows:bycase.setdefault(r['CaseID'],{})[r['Method']]=r
 assert len(bycase)==60 and all(set(v)=={'FixedSplit','DecoupledWeak'} for v in bycase.values())
 for v in bycase.values():
  assert all(v['FixedSplit'][k]==v['DecoupledWeak'][k] for k in ['Q1','Q2','Q3','Q4','VelocityScale','EstimatedFrequency'])
 pair=[]
 for case,rr in bycase.items():
  old=rr['FixedSplit'];new=rr['DecoupledWeak'];pair.append({k:old[k] for k in ['CaseID','Model','Scenario','SNR_dB','Seed']}|{
   'FixedSplitError':float(old['RelativeError']),'NewError':float(new['RelativeError']),
   'ImprovementPercent':100*(1-float(new['RelativeError'])/float(old['RelativeError']))})
 writecsv(dataout/'paired_comparison.csv',pair)
 groups=[]
 for model in read(run/'matlab/executed_plan.json')['models']:
  for snr in [18,12]:
   a=[r for r in rows if r['Model']==model and float(r['SNR_dB'])==snr and r['Method']=='FixedSplit']
   b=[r for r in rows if r['Model']==model and float(r['SNR_dB'])==snr and r['Method']=='DecoupledWeak']
   groups.append({'Model':model,'SNR_dB':snr,'n':len(a),'FixedSplit':avg(a,'RelativeError'),'DecoupledWeak':avg(b,'RelativeError'),
    'ImprovementPercent':100*(1-avg(b,'RelativeError')/avg(a,'RelativeError')),
    'OldAboveOne':sum(float(r['RelativeError'])>1 for r in a),'NewAboveOne':sum(float(r['RelativeError'])>1 for r in b)})
 writecsv(dataout/'model_snr_summary.csv',groups)
 local=[]
 for scene in SCENES:
  for method in ['FixedSplit','DecoupledWeak']:
   rs=[r for r in rows if r['Scenario']==scene and r['Method']==method]
   local.append({'Scenario':scene,'Method':method,'n':len(rs),**{k:avg(rs,k) for k in ['RelativeError','EventAmplitudeError','AVOGradientError','WeakEventError','FaultNeighborhoodError','PinchoutTailError']}})
 writecsv(dataout/'scenario_local_summary.csv',local)
 files=[];gallery=[];rawmanifest=[]
 def export(fig,stem,caption):
  for ext in ['png','pdf']:
   p=dest/(stem+'.'+ext);fig.savefig(p,dpi=180,facecolor='white')
   item={'file':p.name,'sha256':sha(p)}
   if ext=='png':
    with Image.open(p) as im:item|={'pixels':list(im.size),'dpi':list(im.info.get('dpi',[]))}
   files.append(item)
  (dest/(stem+'.txt')).write_text(caption,encoding='utf-8');plt.close(fig)
 style={'font.family':'sans-serif','font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'font.size':9,'pdf.fonttype':42,
  'axes.spines.top':False,'axes.spines.right':False,'savefig.transparent':False}
 with plt.rc_context(style):
  for number,(case,rr) in enumerate(bycase.items(),1):
   p=run/'matlab/cases'/(case+'.json');d=read(p);m=d['meta'];t=np.asarray(d['t']).reshape(-1)
   arrays=[np.asarray(d[k],float) for k in ['truth','observed','fixed','new']];assert all(a.shape==(601,11) and np.isfinite(a).all() for a in arrays)
   x=np.asarray(d['offsets'] if m['domain']=='Offset' else d['lateral']).reshape(-1)
   dx=x[1]-x[0];extent=[x[0]-dx/2,x[-1]+dx/2,t[-1],t[0]]
   lim=max(float(np.max(np.abs(a))) for a in arrays);errors=[arrays[2]-arrays[0],arrays[3]-arrays[0]]
   elim=max(float(np.max(np.abs(a))) for a in errors);gain=.42/lim
   xlabel='炮检距 / Offset (m)' if m['domain']=='Offset' else '横向位置 / Lateral position (m)'
   title=f"{m['model']} | {SCENES[m['scenario']]} | {m['snr_db']} dB | seed {m['seed']}"
   old=float(rr['FixedSplit']['RelativeError']);new=float(rr['DecoupledWeak']['RelativeError'])
   labels=['无吸收真值 / Truth','含吸收观测 / Observed',f'原方法 / FixedSplit   error={old:.4f}',f'新方法 / DecoupledWeak   error={new:.4f}', '原方法误差 / FixedSplit - truth','新方法误差 / New - truth']
   fig,axs=plt.subplots(3,2,figsize=(11,10),layout='constrained',sharex=True,sharey=True)
   for i,(ax,a,label) in enumerate(zip(axs.flat,arrays+errors,labels)):
    bound=lim if i<4 else elim
    im=ax.imshow(a,cmap='RdBu_r',vmin=-bound,vmax=bound,extent=extent,aspect='auto',interpolation='nearest')
    ax.set(title=label,xlabel=xlabel,ylabel='时间 / Time (s)',xticks=x[::2]);ax.set_ylim(t[-1],t[0])
    if i==0:mainim=im
    if i==4:errim=im
   fig.colorbar(mainim,ax=axs[:2,:].ravel().tolist(),label='共享振幅 / Shared amplitude',shrink=.85)
   fig.colorbar(errim,ax=axs[2,:].ravel().tolist(),label='共享误差振幅 / Shared error amplitude',shrink=.85)
   fig.suptitle(title,fontsize=11)
   cap=f'{title}。上四图共同色标±{lim:.8g}，下两幅误差图共同色标±{elim:.8g}；两类色标分别标示。全时间范围，无逐道归一化、滤波、截幅或插值。Section为固定750m炮检距的横向剖面。真值只用于评分和画图。'
   export(fig,case+'_gather',cap)
   fig,axs=plt.subplots(1,4,figsize=(15,7),layout='constrained',sharey=True)
   for ax,a,label in zip(axs,arrays,labels[:4]):
    for j in range(11):
     xx=j+gain*a[:,j];ax.plot(xx,t,color='#222222',lw=.45)
     ax.fill_betweenx(t,j,xx,where=a[:,j]>=0,color='#222222',linewidth=0)
    ax.set(title=label,xlim=(-.65,10.65),ylim=(t[-1],t[0]),xticks=np.arange(0,11,2),xticklabels=[str(int(z)) for z in x[::2]],xlabel=xlabel)
    ax.grid(axis='y',alpha=.15)
   axs[0].set_ylabel('时间 / Time (s)');fig.suptitle(title+f'\n同组统一增益 / Common gain = {gain:.5g}; no trace normalization',fontsize=11)
   export(fig,case+'_wiggle',f'{title}。四面板共同波形增益{gain:.10g}=0.42/同组四道集最大绝对振幅；黑色填充表示正振幅。无逐道或逐方法归一化，无振幅裁剪，完整时间范围。')
   fig,axs=plt.subplots(1,2,figsize=(12,4),layout='constrained')
   for ax in axs:
    for a,label,color,lw,ls in zip(arrays,['Truth','Observed','FixedSplit','DecoupledWeak'],['black','#999999','#A35000','#0072B2'],[1.5,.65,1.1,1.1],['-','-', '--','-.']):
     ax.plot(t,a[:,5],label=label,color=color,lw=lw,ls=ls)
    ax.set(xlabel='时间 / Time (s)',ylabel='原始振幅 / Amplitude');ax.grid(alpha=.15)
   axs[0].set_title('第6道全记录 / Trace 6, full record');axs[1].set(title='固定时间细节 / Fixed time detail',xlim=(.45,1.0))
   axs[0].legend(fontsize=8);fig.suptitle(title,fontsize=10)
   export(fig,case+'_trace',f'{title}。预定第6道，原始绝对振幅。左图全部采样时间，右图统一取0.45–1.0s；未重新缩放曲线。')
   gallery.append({'id':case,'model':m['model'],'scenario':SCENES[m['scenario']],'snr':m['snr_db'],'seed':m['seed'],'old':old,'new':new,'limit':lim,'error_limit':elim,'wiggle_gain':gain})
   rawmanifest.append({'file':case+'.json','sha256':sha(p),'observed_array_sha256':hashlib.sha256(arrays[1].astype('<f8').tobytes()).hexdigest()})
   if number%10==0:print(f'Figures {number}/60',flush=True)
  for model in read(run/'matlab/executed_plan.json')['models']:
   firstseed=read(run/'matlab/executed_plan.json')['seeds'][0]
   fig,axs=plt.subplots(5,4,figsize=(14,12),layout='constrained',sharey=True)
   rowlimits=[]
   for k,scene in enumerate(SCENES):
    chosen=next(g for g in gallery if g['model']==model and g['snr']==18 and g['seed']==firstseed and g['scenario']==SCENES[scene])
    d=read(run/'matlab/cases'/(chosen['id']+'.json'));t=np.asarray(d['t']).reshape(-1)
    x=np.asarray(d['offsets'] if d['meta']['domain']=='Offset' else d['lateral']).reshape(-1);dx=x[1]-x[0]
    aa=[np.asarray(d[key],float) for key in ['truth','observed','fixed','new']];lim=chosen['limit'];rowlimits.append(lim)
    for j,(ax,a) in enumerate(zip(axs[k],aa)):
     im=ax.imshow(a,cmap='RdBu_r',vmin=-lim,vmax=lim,extent=[x[0]-dx/2,x[-1]+dx/2,t[-1],t[0]],aspect='auto',interpolation='nearest')
     ax.set(xticks=x[::5],ylim=(t[-1],t[0]));ax.tick_params(labelsize=7)
     if k==0:ax.set_title(['无吸收真值 / Truth','含吸收观测 / Observed','最初 / FixedSplit','新 / DecoupledWeak'][j])
     if j==0:ax.set_ylabel(SCENES[scene]+'\nTime (s)')
     ax.set_xlabel('Offset (m)' if d['meta']['domain']=='Offset' else 'Lateral position (m)',fontsize=7)
    fig.colorbar(im,ax=axs[k].tolist(),label='Amplitude',shrink=.8)
   fig.suptitle(f'{model} | 18dB | seed {firstseed} | 五场景总览 / All five scenarios',fontsize=12)
   export(fig,'overview_'+model+'_18dB',f'预定主图：{model}，18dB，第一个种子{firstseed}，五场景全部展示，不按误差挑选。各行四幅同一振幅色标，五行对称幅值分别为{rowlimits}，每行均画出色标。完整时间范围，无归一化或滤波；前两行横轴炮检距，后三行横轴剖面位置（固定750m炮检距）。')
  fig,ax=plt.subplots(figsize=(7,6),layout='constrained');colors=['#0072B2','#A35000','#007A45']
  names=read(run/'matlab/executed_plan.json')['models']
  for model,color in zip(names,colors):
   for snr,marker in [(18,'o'),(12,'s')]:
    pp=[r for r in pair if r['Model']==model and float(r['SNR_dB'])==snr]
    ax.scatter([r['FixedSplitError'] for r in pp],[r['NewError'] for r in pp],label=f'{model}, {snr}dB',color=color,marker=marker,s=35,alpha=.85)
  vmax=max(max(r['FixedSplitError'],r['NewError']) for r in pair)*1.04
  ax.plot([0,vmax],[0,vmax],color='black',ls='--',lw=1);ax.set(xlim=(0,vmax),ylim=(0,vmax),xlabel='原方法相对误差 / FixedSplit error',ylabel='新方法相对误差 / New error',title='全部60组配对 / All 60 pairs; below diagonal favors new')
  ax.set_aspect('equal');ax.legend(fontsize=7);ax.grid(alpha=.15)
  export(fig,'all_paired_errors','全部60组观测，无删样本。每点为相同观测与估计模型下的两方法误差；虚线下方新方法更好。种子在场景间复用，不把60点视为独立地质重复，无显著性检验。')
 (dest/'manifest.json').write_text(json.dumps({'source_run':str(run),'raw_data':rawmanifest,'figures':files,'case_scales':gallery,'python':platform.python_version(),'matplotlib':matplotlib.__version__,'numpy':np.__version__},ensure_ascii=False,indent=2),encoding='utf-8')
 (dataout/'fairness_checks.json').write_text(json.dumps({'helpers_byte_identical':helpers,'paired_observations':60,'same_calibration_parameters_verified':True,'source_snapshot_changed_files':0,'numerical_checks':checks},ensure_ascii=False,indent=2),encoding='utf-8')
 build_gallery(gallery)
 a=[r for r in rows if r['Method']=='FixedSplit'];b=[r for r in rows if r['Method']=='DecoupledWeak'];overall=100*(1-avg(b,'RelativeError')/avg(a,'RelativeError'))
 lines=['# 最初方法与新方法：同观测多场景比较','','## Material Passport','','- Origin Skill: academic-research-suite / experiment-agent','- Origin Mode: run','- Verification Status: ANALYZED; paired inputs and numerical checks verified','- Origin Date: 2026-10-22','- Version Label: fixed_vs_new_gathers_v1','',
  f'已完成60组观测、120条方法结果。FixedSplit与DecoupledWeak(0.3)共用每组同一个D及同一份校准所得Q、速度和频率。整体等权误差：原方法{avg(a,"RelativeError"):.6f}，新方法{avg(b,"RelativeError"):.6f}；新方法相对改善{overall:+.2f}%。新方法在{sum(r["NewError"]<r["FixedSplitError"] for r in pair)}/60组误差更低。','',
  '## 分模型、噪声结果','','| 模型 | SNR | 原方法误差 | 新方法误差 | 相对改善 | 原/新误差>1 |','|---|---:|---:|---:|---:|---:|']
 for g in groups:lines.append(f'|{g["Model"]}|{g["SNR_dB"]}|{g["FixedSplit"]:.6f}|{g["DecoupledWeak"]:.6f}|{g["ImprovementPercent"]:+.2f}%|{g["OldAboveOne"]}/{g["NewAboveOne"]}|')
 lines+=['','## 逐场景结果（模型、SNR、种子等权）','','| 场景 | 原方法 | 新方法 | 相对改善 |','|---|---:|---:|---:|']
 for scene in SCENES:
  x=next(r for r in local if r['Scenario']==scene and r['Method']=='FixedSplit');y=next(r for r in local if r['Scenario']==scene and r['Method']=='DecoupledWeak')
  lines.append(f'|{SCENES[scene]}|{x["RelativeError"]:.6f}|{y["RelativeError"]:.6f}|{100*(1-y["RelativeError"]/x["RelativeError"]):+.2f}%|')
 lines+=['','## 局部指标（仅适用场景）','','| 指标 | 原方法 | 新方法 |','|---|---:|---:|']
 for key in ['EventAmplitudeError','AVOGradientError','WeakEventError','FaultNeighborhoodError','PinchoutTailError']:lines.append(f'|{key}|{avg(a,key):.6f}|{avg(b,key):.6f}|')
 lines+=['','## 数值与显示核验','',f'新方法全部选参候选最终标志失败数{checks["new_final_solver_flags"]}，最大法方程残差{checks["max_new_normal_residual"]:.3e}；原两步最大法方程残差{checks["max_old_normal_residual"]:.3e}，原两步矩阵复合重放差{checks["max_old_composition_replay"]:.3e}。相同D断言通过60次，逐组校准参数相等，三个基础数值函数与最初归档逐字节相同，冻结源码哈希未变。','',
  '每组均提供色彩图、黑白波形道集及第6道波形。主图按计划第一模型、第一种子、18dB展示；图库包含全部模型、场景、噪声与种子。色彩图上四图同振幅尺度，下面两个误差图采用另外一个共用误差尺度；波形道集采用同组统一增益。没有逐道归一化或挑选有利样例。','',
  '## 比较范围','',
  'FixedSplit使用最初归档的两步固定参数公式，但接收与新方法相同的校准估计模型；本轮没有给旧方法真实Q，也不复现缺失的毕设原代码。结果评价的是共享校准信息下两套恢复流程，不是只比较新正则项，也不代表分步流程优于一步法。','新方法还包含自适应选参、噪声预处理与多道结构信息。其额外模块、计算预算均未与旧法逐项消融；不将全部收益归为权重解耦。','合成采用Shuey二项、理想NMO和既有矩阵模拟器；横向Section为固定750m采集几何，非真实迁移成像剖面。每模型一份校准、两个观测种子，场景间复用种子；均值和胜率为描述统计，不作显著性或实测泛化声明。','',
  '## 文件入口','','- index.html：离线道集图库，可切换全部60组。','- figures/：全部PNG、PDF及逐图说明；manifest.json记录幅值尺度、原始数据和图像SHA256。','- analysis/：逐观测配对、分模型和局部指标、输入公平性检查。','- runs/ 下本批 matlab/cases/：全部原始绘图数组JSON与MAT；metrics.csv为全部120行。','- core_p1/run_fixed_vs_new.m与inputs/fixed_vs_new_plan.json：试验入口和冻结设计。','',f'实际运行：{run}；MATLAB耗时{checks["elapsed_seconds"]:.1f}s。独立新建进程，不操作用户原有MATLAB会话。']
 (ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 add_interpretation(rows)
 print(json.dumps({'overall_old':avg(a,'RelativeError'),'overall_new':avg(b,'RelativeError'),'improvement_percent':overall,'new_wins':sum(r['NewError']<r['FixedSplitError'] for r in pair),'figure_count':len(files),'report':str(ROOT/'RESULTS_ZH.md')},ensure_ascii=False),flush=True)
def add_interpretation(rows):
 a=[r for r in rows if r['Method']=='FixedSplit'];b=[r for r in rows if r['Method']=='DecoupledWeak']
 weak=100*(avg(b,'WeakEventError')/avg(a,'WeakEventError')-1)
 avo=100*(avg(b,'AVOGradientError')/avg(a,'AVOGradientError')-1)
 with (ROOT/'RESULTS_ZH.md').open('a',encoding='utf-8') as f:
  f.write('\n## 评价与对照区别\n\n')
  f.write(f'本轮整体波形改善并未覆盖全部保幅指标：弱反射振幅误差平均增加{weak:.2f}%，AVO梯度误差平均增加{avo:.2f}%。全部60组胜出只指整体波形误差，不表示每个事件、频段或局部指标都改善。\n')
  f.write('这次基准是最初固定参数FixedSplit，上一轮0.75%的比较基准则是已包含多项改进的Legacy；两者不是同一个问题，也不能跨不同模型、噪声设置直接比较误差均值。两个Dev模型此前用于开发，本轮的新校准与新观测种子不能使其成为完全独立的模型确认。\n')
  f.write('尖灭末端三个道的第二、第三事件重合，既有模板评分将每组9个系数标记为不可辨识，并从EventAmplitudeError中排除；具体计数在metrics.csv的UnidentifiableAmplitudes列。整体波形和尖灭尾部评分仍包含全部道集。弱反射场景没有不可辨识标记，上述弱反射误差变化不受这项排除影响。\n')

def build_gallery(items):
 content='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>最初方法与新方法｜道集对比</title><style>body{font:16px system-ui,sans-serif;background:#f4f6fa;color:#202d40;margin:0}main{max-width:1400px;margin:auto;padding:28px}h1{font-size:26px}select,a{padding:9px;margin:4px;border:1px solid #bfcada;border-radius:6px;background:white}img{width:100%;background:white;margin-top:15px;border-radius:8px}p{line-height:1.6}#info{padding:12px;background:white;border-radius:6px}footer{font-size:13px;color:#536174}</style><main><h1>最初 FixedSplit 与新 DecoupledWeak：60组相同观测</h1><p>两方法共用校准估计的Q和速度。真值只用于合成与评分。色彩图上四幅同振幅尺度，下两幅另用共同误差尺度；波形道集不逐道归一化。所有种子完整保留。</p><label>场景 <select id="case"></select></label><label>显示 <select id="kind"><option value="gather">色彩道集与误差</option><option value="wiggle">黑白波形道集</option><option value="trace">固定第6道波形</option></select></label><a id="pdf">下载当前PDF</a><div id="info"></div><img id="figure" alt="相同观测下原方法和新方法的道集对比"><footer>详细结果见 RESULTS_ZH.md；图件数据、显示尺度及来源见 figures/manifest.json。此为有限合成研究比较。</footer></main><script>const items=ITEMS;const c=document.getElementById('case'),k=document.getElementById('kind');items.forEach((x,i)=>c.add(new Option(`${x.model} · ${x.scenario} · ${x.snr}dB · ${x.seed}`,i)));function show(){const x=items[+c.value],base='figures/'+x.id+'_'+k.value;document.getElementById('figure').src=base+'.png';document.getElementById('pdf').href=base+'.pdf';document.getElementById('info').textContent=`原方法误差 ${x.old.toFixed(5)} → 新方法 ${x.new.toFixed(5)}；相对改善 ${(100*(1-x.new/x.old)).toFixed(2)}%。同组主图幅值±${x.limit.toPrecision(4)}，误差图±${x.error_limit.toPrecision(4)}。`;}c.onchange=k.onchange=show;show();</script></html>'''
 (ROOT/'index.html').write_text(content.replace('ITEMS',json.dumps(items,ensure_ascii=False)),encoding='utf-8')
if __name__=='__main__':main()
