from pathlib import Path
import json,csv,hashlib,sys,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
truth=np.load(ROOT/'truth/model_truth.npz');x=truth['x'];z=truth['z'];t=truth['t'];live=truth['live'];plan=json.loads((ROOT/'PLAN.json').read_text(encoding='utf-8'))
payload=json.loads((ROOT/'core_p1/inputs/complex_input.json').read_text(encoding='utf-8'))
fig,axs=plt.subplots(3,1,figsize=(14,12),sharex=True,sharey=True,layout='constrained')
for ax,key,title,cmap,unit in zip(axs,['vp','rho','Q'],['纵波速度模型','密度模型','真实 Q 模型'],['turbo','viridis','viridis'],['m/s','g/cm³','Q']):
 a=truth[key];im=ax.imshow(a,extent=[x[0]/1000,x[-1]/1000,z[-1],z[0]],aspect='auto',interpolation='nearest',cmap=cmap,vmin=float(a.min()),vmax=float(a.max()))
 ax.set_title(title);ax.set_ylabel('深度 / m');fig.colorbar(im,ax=ax,label=unit,shrink=.9)
 for xpos in [1.26,2.28]:ax.axvline(xpos,color='white',lw=.6,ls='--',alpha=.65)
axs[0].contour(x/1000,z,truth['channel_mask'].astype(float),levels=[.5],colors='white',linewidths=.8)
axs[0].contour(x/1000,z,truth['gas_mask'].astype(float),levels=[.5],colors='black',linewidths=.8)
axs[0].annotate('侵蚀河道',xy=(1.71,930),xytext=(1.6,500),color='white',arrowprops=dict(arrowstyle='->',color='white'),bbox=dict(facecolor='black',alpha=.6,edgecolor='none'))
axs[0].annotate('低 Q 透镜体',xy=(2.08,1320),xytext=(2.35,1750),color='white',arrowprops=dict(arrowstyle='->',color='white'),bbox=dict(facecolor='black',alpha=.6,edgecolor='none'))
axs[-1].set_xlabel('横向位置 / km');fig.suptitle('复杂合成模型：密集层系、褶皱、双断层、河道、薄层尖灭和低 Q 透镜体',fontsize=16)
fig.supxlabel('99 道 × 30 m；深度网格 5 m；白色虚线标记断层横向位置。此图为真实地质属性，不是补偿结果。',fontsize=9)
fig.savefig(OUT/'complex_geology.png',dpi=170);fig.savefig(OUT/'complex_geology.pdf');plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(16,7),sharey=True,layout='constrained')
target=truth['Moderate_target'];obs=[np.asarray(d['D']) for d in payload['datasets']];limit=max(np.max(np.abs(a)) for a in [target]+obs)
for ax,a,title in zip(axs,[target]+obs,['14 dB 条件的匹配纯净参考','复杂观测：14 dB 随机噪声','复杂观测：8 dB 随机噪声']):
 im=ax.imshow(a,extent=[x[0]/1000,x[-1]/1000,t[-1],t[0]],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_title(title,fontsize=11);ax.set_xlabel('横向位置 / km')
axs[0].set_ylabel('双程时间 / s');fig.colorbar(im,ax=axs,label='三个面板共用振幅色标',shrink=.85)
fig.suptitle('原始资料外观：非平稳衰减、混合相位、多次波、相干干扰、相关噪声和坏道',fontsize=15)
fig.supxlabel('随机噪声 SNR 不包含相干干扰、多次波和坏道损失；两种观测的增益/静校正不同，各自有匹配纯净目标。',fontsize=9)
fig.savefig(OUT/'clean_and_observed.png',dpi=160);fig.savefig(OUT/'clean_and_observed.pdf');plt.close(fig)
if '--model-only' in sys.argv:print('Complex geology and input figures complete.',flush=True);raise SystemExit(0)

RUN=Path(json.loads((ROOT/'run.json').read_text(encoding='utf-8'))['run']);DATA=RUN/'matlab'
checks=json.loads((DATA/'checks.json').read_text(encoding='utf-8'));assert checks['complete'] and checks['blocks']==18
expected=json.loads((ROOT/'FROZEN_SOLVER_HASHES.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((RUN/'source/core_p1'/p).read_bytes()).hexdigest()==h for p,h in expected.items())
rows=[];blockrows=[];manifest=[];methods=dict(Observed='observed',FixedSplit='fixed',FixedSplitSharedCal='fixed_shared',FullCalRefit='unpruned',PrunedFullCalRefit='optimized')
def score(s,target,weakmask,live_mask=live):
 e=s-target;deep=(t>=1.05)&(t<=1.6)
 a=target[:,live_mask].ravel();b=s[:,live_mask].ravel();a=a-a.mean();b=b-b.mean()
 return dict(RelativeError=float(np.linalg.norm(e)/np.linalg.norm(target)),LiveTraceError=float(np.linalg.norm(e[:,live_mask])/np.linalg.norm(target[:,live_mask])),CorrelationLive=float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))) if np.linalg.norm(b)>0 else 0,
  AmplitudeProjectionGain=float(np.sum(s[:,live_mask]*target[:,live_mask])/np.sum(target[:,live_mask]**2)),DeepLiveError=float(np.linalg.norm(e[deep,:][:,live_mask])/np.linalg.norm(target[deep,:][:,live_mask])),LowEnergyRegionError=float(np.linalg.norm(e[weakmask])/np.linalg.norm(target[weakmask])),PreArrivalErrorRMS=float(np.sqrt(np.mean(e[t<.25,:][:,live_mask]**2))))
for cfg in plan['observations']:
 name=cfg['name'];path=DATA/(name+'.json');d=json.loads(path.read_text(encoding='utf-8'));target=truth[name+'_target']
 power=np.apply_along_axis(lambda a:np.convolve(a*a,np.ones(15)/15,mode='same'),0,target);local=np.sqrt(power)
 region=(t[:,None]>=.30)&(t[:,None]<=1.6)&live[None,:];nonzero=region&(local>.001*local.max());lo,hi=np.quantile(local[nonzero],[.20,.45]);weakmask=region&(local>=lo)&(local<=hi)
 np.savez_compressed(ROOT/'analysis'/(name+'_evaluation_masks.npz'),low_energy=weakmask,live=live,deep=(t>=1.05)&(t<=1.6),local_rms=local,low_energy_thresholds=np.array([lo,hi]))
 for method,key in methods.items():
  s=np.asarray(d[key]);rows.append(dict(Dataset=name,NominalRandomSNR_dB=cfg['snr_db'],Method=method,**score(s,target,weakmask)))
  for k in range(9):
   idx=slice(k*11,(k+1)*11);blockrows.append(dict(Dataset=name,Block=k+1,XStart_m=x[k*11],XEnd_m=x[(k+1)*11-1],Method=method,**score(s[:,idx],target[:,idx],weakmask[:,idx],live[idx])))
 fixed=np.asarray(d['fixed']);new=np.asarray(d['optimized']);observed=np.asarray(d['observed']);arrays=[target,fixed,new,observed,fixed-target,new-target]
 amp=max(float(np.abs(a).max()) for a in arrays[:4]);err=max(float(np.abs(a).max()) for a in arrays[4:])
 fig=plt.figure(figsize=(16,10),layout='constrained');gs=fig.add_gridspec(2,4,width_ratios=[1,1,1,.045])
 for k,(a,title) in enumerate(zip(arrays,['匹配纯净参考','最初分步法','优化方法','复杂观测','分步法误差','优化方法误差'])):
  ax=fig.add_subplot(gs[k//3,k%3]);limit=amp if k<4 else err
  im=ax.imshow(a,extent=[0,x[-1]/1000,t[-1],0],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_title(title);ax.set_xlabel('横向位置 / km')
  if k%3==0:ax.set_ylabel('双程时间 / s')
  if k<4:ampim=im
  else:errim=im
 fig.colorbar(ampim,cax=fig.add_subplot(gs[0,3]),label='四个道集共用振幅色标');fig.colorbar(errim,cax=fig.add_subplot(gs[1,3]),label='两个误差面板共用误差色标')
 fig.suptitle(f'复杂模型对比 · {name} · 随机噪声 {cfg["snr_db"]} dB',fontsize=16)
 fig.supxlabel('同一观测、同一公开几何、独立校准；不逐道归一化。第23、78道为坏道，保留在全道误差中，另报正常道误差。',fontsize=9)
 dest=OUT/(name+'_comparison.png');fig.savefig(dest,dpi=155);fig.savefig(OUT/(name+'_comparison.pdf'));plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(15,8),sharex=True,sharey=True,layout='constrained');idx=np.arange(0,len(x),3);gain=.40/max(float(np.abs(a[:,idx]).max()) for a in [target,fixed,new])
 for ax,a,title in zip(axs,[target,fixed,new],['纯净参考','最初分步法','优化方法']):
  for k,j in enumerate(idx):
   line=k+1+gain*a[:,j];ax.plot(line,t,'k',lw=.4);ax.fill_betweenx(t,k+1,line,where=a[:,j]>0,color='k',alpha=.65)
  ax.set_xlim(0,len(idx)+1);ax.set_ylim(t[-1],0);ax.set_title(title);ax.set_xlabel('每3道抽取1道的序号')
 axs[0].set_ylabel('双程时间 / s');fig.suptitle(name+' · 同增益波形对比');fig.supxlabel(f'全99道见热图；此处固定每3道取1道，共33道；全部面板与各道的显示增益相同：{gain:.5g}。',fontsize=9)
 fig.savefig(OUT/(name+'_wiggles.png'),dpi=150);fig.savefig(OUT/(name+'_wiggles.pdf'));plt.close(fig)
 freqs=np.fft.rfftfreq(len(t),plan['dt_s']);fig,ax=plt.subplots(figsize=(9,5),layout='constrained');ref=np.mean(np.abs(np.fft.rfft(target[:,live],axis=0))**2,axis=1);refmax=ref.max()
 for a,label,color in [(target,'纯净参考','black'),(observed,'观测','#009E73'),(fixed,'分步法','#777777'),(new,'优化方法','#D55E00')]:
  ps=np.mean(np.abs(np.fft.rfft(a[:,live],axis=0))**2,axis=1);ax.plot(freqs,10*np.log10(np.maximum(ps/refmax,1e-12)),label=label,color=color)
 ax.set_xlim(0,120);ax.set_ylim(-80,20);ax.set_xlabel('频率 / Hz');ax.set_ylabel('平均功率谱 / dB（相对纯净参考峰值）');ax.legend();ax.grid(alpha=.2);ax.set_title(name+' · 正常道平均功率谱')
 fig.savefig(OUT/(name+'_spectra.png'),dpi=160);fig.savefig(OUT/(name+'_spectra.pdf'));plt.close(fig)
 manifest.append(dict(dataset=name,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),gather=str(dest),amplitude_limit=amp,error_limit=err,sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))
with (ROOT/'analysis/metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
with (ROOT/'analysis/block_metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(blockrows[0]));writer.writeheader();writer.writerows(blockrows)
block_wins={}
fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True,layout='constrained')
for ax,name in zip(axs,['Moderate','Strong']):
 position=np.arange(1,10);subset=[r for r in blockrows if r['Dataset']==name]
 for delta,method,label,color in [(-.18,'FixedSplit','原分步法','#777777'),(.18,'PrunedFullCalRefit','优化方法','#D55E00')]:
  vals=[next(r['LiveTraceError'] for r in subset if r['Method']==method and r['Block']==k) for k in position];ax.bar(position+delta,vals,.35,label=label,color=color)
 ax.set_title(name+'：逐块正常道误差');ax.set_ylabel('相对波形误差');ax.set_ylim(bottom=0);ax.legend();ax.grid(axis='y',alpha=.2)
 block_wins[name]=sum(next(r['LiveTraceError'] for r in subset if r['Method']=='PrunedFullCalRefit' and r['Block']==k)<next(r['LiveTraceError'] for r in subset if r['Method']=='FixedSplit' and r['Block']==k) for k in position)
axs[-1].set_xticks(np.arange(1,10),[f'{x[k*11]/1000:.2f}–{x[(k+1)*11-1]/1000:.2f}' for k in range(9)],rotation=20);axs[-1].set_xlabel('固定11道区块的横向范围 / km');fig.supxlabel('区块仅用于同一模型内的空间诊断，不作为独立地质模型。',fontsize=9)
fig.savefig(OUT/'block_comparison.png',dpi=160);fig.savefig(OUT/'block_comparison.pdf');plt.close(fig)
summary=dict(checks=checks,metrics=rows,block_wins_out_of_9=block_wins,source_files_unchanged=len(expected),manifest=manifest,scope='One complex geological model, two noise conditions; no independent real-data validation.')
(ROOT/'analysis/summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 复杂模型实际对比结果','','## Material Passport','','类型：复杂合成地震补偿实验。原有求解代码冻结，实际MATLAB运行状态和日志保留在runs/。','','## 模型与正演','','横向2.94 km、99道、深度2.2 km、5 m地质网格；时间1.7 s、2 ms采样。60组不规则密集层位、褶皱、双断层、侵蚀河道、薄层尖灭、低Q透镜体和微小非均质。','真实Vp、Vs、密度、Q均有横向变化。主频36Hz随位置变化±7%，源相位24±7度，叠加12ms/25ms源端延迟分量；包括近似绕射、多次波、相关噪声、面波形相干干扰、突发噪声、静校正、接收增益及两条坏道。','正演采用逐列二维属性驱动的局部一维反射/卷积，累计dt/Q和RMS直射线路径因子构建非平稳衰减与频散；绕射和多次波为运动学近似。没有模拟完整二维弹性波散射、模式转换或自由表面多次波。不能当作实测资料或完整弹性/粘弹性正演。','','## 公平性与评价目标','','原分步法、同校准分步法、未剪枝与最终优化使用完全相同观测。真实Q/属性模型不进入求解器。独立18dB标定采用已知井约束的粗四层几何，估计的单组Q/主频迁移到整条剖面，并不逐位置使用真实Q。','99道按固定相邻11道分为9块，复用既有准备模块；两种观测共18块。分步法固定双岭参数1e-3，优化稀疏和组筛选规则不变。固定偏移距450m剖面不进行跨横向AVO消除。','纯净目标保留对应观测的真实源相位、接收增益和静校正；不包含吸收衰减、额外多次波、噪声、坏道。全道评分保留两条坏道的损失，正常道评分按事先固定的第23和78道质量标记排除。','低能量区域按真值局部RMS在信号窗内的20%—45%分位预定义，只用于事后评价，不传入算法；它是低能量区域波形误差，不是可单独辨识的弱反射系数误差。','14/8dB仅表示随机噪声相对吸收后主波的SNR，额外相干干扰和多次波未包含在该数值中。两个条件来自同一地质模型，不当作两个独立模型。','','## 实际指标','','| 条件 | 方法 | 全道误差 | 正常道误差 | 正常道相关系数 | 振幅投影增益 | 深部误差 | 低能量区误差 |','|---|---|---:|---:|---:|---:|---:|---:|']
for r in rows:lines.append('| '+r['Dataset']+' | '+r['Method']+' | '+' | '.join(f'{r[k]:.5g}' for k in ['RelativeError','LiveTraceError','CorrelationLive','AmplitudeProjectionGain','DeepLiveError','LowEnergyRegionError'])+' |')
lines+=['','## 结果解释','']
for name in ['Moderate','Strong']:
 a=next(r for r in rows if r['Dataset']==name and r['Method']=='FixedSplit');b=next(r for r in rows if r['Dataset']==name and r['Method']=='PrunedFullCalRefit');lines.append(f'{name}：优化方法相对原分步法的全道误差变化为{100*(b["RelativeError"]/a["RelativeError"]-1):+.2f}%；9个固定区块中，正常道误差降低{block_wins[name]}/9。低能量区误差原分步法{a["LowEnergyRegionError"]:.4g}，优化方法{b["LowEnergyRegionError"]:.4g}。')
 u=next(r for r in rows if r['Dataset']==name and r['Method']=='FullCalRefit');lines.append(f'同一观测在反射组筛选前后的全道误差为{u["RelativeError"]:.4g}→{b["RelativeError"]:.4g}，低能量区误差为{u["LowEnergyRegionError"]:.4g}→{b["LowEnergyRegionError"]:.4g}。这次组筛选牺牲了恢复精度，不能把更干净的背景等同于更好保幅。')
lines+=['','粗四层几何由真实属性的空间平均构造，作为假定已知的速度/层厚先验；没有模拟测井或速度建模误差，不代表实际测井资料。独立标定采用36Hz名义源及模型中点相位，同族但不严格等于每个位置的真实源波形。','本轮为接近实测困难的受控复杂合成试验，不能代替真实数据或完整二维波动方程验证。比较保持模型和算法冻结，所有结果包括不利区域均保留。']
cal=json.loads((DATA/'calibration.json').read_text(encoding='utf-8'));e=cal['estimate'];lines+=['','## 标定与数值检查','',f'估计Q：{e.get("layer_Q",[e["effective_Q"]]*4)}；估计主频：{e["estimated_frequency"]:.4g} Hz。',f'稀疏阶段未收敛标志总数：{checks["optimization_flags"]}；振幅重估最大正规方程相对残差：{checks["max_normal"]:.4g}。',f'{len(expected)}个原有MATLAB文件与上一轮哈希相同。冻结源码、观测输入、独立标定、原始结果和真值分目录保留。','未根据这次结果调整参数。参考、观测和补偿道集采用统一色标；误差面板共用色标；波形图固定每3道取1道且使用同一显示增益。','完整参数在PLAN.json，生成脚本scripts/build_model.py，执行入口core_p1/run_complex_model.m，绘图与评分scripts/plot_and_report.py。']
report='\n'.join(lines)+'\n'
report=report.replace('独立18dB标定采用已知井约束的粗四层几何','独立18dB标定采用假定已知的粗四层几何（由属性空间平均构造）')
report+='\nPLAN中的coherent_energy_fraction字段实际表示相干干扰与基准观测的Frobenius振幅范数比（0.18/0.30），不是功率或能量比。\n'
(ROOT/'RESULTS_ZH.md').write_text(report,encoding='utf-8')
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>复杂模型补偿对比</title><style>body{font:16px system-ui;max-width:1500px;margin:24px auto;padding:0 20px}img{width:100%}p{line-height:1.7}</style><h1>复杂模型补偿对比</h1><p>同一复杂地质模型的两种噪声条件。完整结果均保留，没有按结果选择。<a href="RESULTS_ZH.md">实验报告</a> · <a href="analysis/metrics.csv">全部指标</a></p><img src="figures/complex_geology.png" alt="复杂地质属性模型"><img src="figures/clean_and_observed.png" alt="参考与原始观测">'''
page+='<img src="figures/block_comparison.png" alt="逐区块误差对比">'
for name in ['Moderate','Strong']:page+=f'<h2>{name}</h2><img src="figures/{name}_comparison.png" alt="道集与误差比较"><img src="figures/{name}_wiggles.png" alt="同增益波形比较"><img src="figures/{name}_spectra.png" alt="平均功率谱">'
(ROOT/'index.html').write_text(page+'</html>',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=True,indent=2),flush=True)
