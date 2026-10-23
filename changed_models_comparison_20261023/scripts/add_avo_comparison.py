"""Add auditable two-term AVO removal to frozen, already-computed gathers.

No seismic compensation solver or calibration is changed or rerun.
The fitted slope is removed, and the non-AVO residual is retained.
"""
from pathlib import Path
import json,csv,hashlib,zipfile,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'avo_comparison'
for folder in [OUT,OUT/'figures',OUT/'data']:folder.mkdir(exist_ok=True)
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
PLAN=json.loads((ROOT/'PLAN.json').read_text(encoding='utf-8'))
RUN=Path(json.loads((ROOT/'run.json').read_text(encoding='utf-8'))['run'])
NAMES=['深层低 Q','薄层极性反转','低阻抗差','高频 Q 反转']

def ray_basis(t,offsets,vp,h):
 """Normalized p(t,x)^2, matching ray_absorption.m; no Q needed."""
 t=np.asarray(t).reshape(-1);vp=np.asarray(vp);offsets=np.asarray(offsets)
 bounds=np.r_[0,np.cumsum(2*np.asarray(h)/vp[:3]),np.inf]
 tau=np.maximum(0,np.minimum(t[:,None],bounds[1:])-bounds[:-1]);active=tau>0
 vmax=np.max(np.where(active,vp,0),axis=1);vmax=np.where(vmax>0,vmax,vp[0])
 low=np.zeros((len(t),len(offsets)));high=np.repeat(((1-1e-12)/vmax)[:,None],len(offsets),axis=1)
 for _ in range(70):
  mid=(low+high)/2
  distance=np.sum(tau[:,None,:]*vp[None,None,:]**2*mid[:,:,None]/np.sqrt(np.maximum(1-(mid[:,:,None]*vp[None,None,:])**2,np.finfo(float).tiny)),axis=2)
  take=distance<offsets[None,:];low=np.where(take,mid,low);high=np.where(take,high,mid)
 p=(low+high)/2;p[:,offsets==0]=0;p[t<=0,:]=0
 q=p*p;scale=np.max(q,axis=1);return q/np.where(scale>0,scale,1)[:,None]

def remove_avo(gather,q):
 """Fit s(t,x)=I(t)+G(t)q(t,x)+e(t,x); return s-Gq.

 Preserves e instead of replacing every trace with the fitted intercept.
 """
 s=np.asarray(gather,dtype=float);assert s.shape==q.shape and np.isfinite(s).all()
 qc=q-q.mean(axis=1,keepdims=True);den=np.sum(qc*qc,axis=1)
 slope=np.divide(np.sum(qc*s,axis=1),den,out=np.zeros(len(s)),where=den>1e-24)
 intercept=np.mean(s,axis=1)-slope*np.mean(q,axis=1)
 removed=slope[:,None]*q;return s-removed,intercept,slope

def zero_angle_truth(model,t):
 """Analytic model reference ONLY, not used by remove_avo."""
 vp=np.asarray(model['vp']);rho=np.asarray(model['rho'])
 A=.5*(np.diff(vp)/((vp[:-1]+vp[1:])/2)+np.diff(rho)/((rho[:-1]+rho[1:])/2))
 z=np.pi*model['fpeak']*(t[:,None]-np.asarray(model['event_t'])[None,:]);w=(1-2*z*z)*np.exp(-z*z)
 return np.repeat((w@A)[:,None],11,axis=1)

rows=[];records=[];overviews=[];basis_checks={}
for model,name in zip(PLAN['models'],NAMES):
 calpath=RUN/'matlab'/('calibration_'+model['name']+'.json');cal=json.loads(calpath.read_text(encoding='utf-8'));public=cal['corrected']
 t=np.asarray(public['t']).reshape(-1);q=ray_basis(t,public['offset'],public['vp'],public['h'])
 # Check physical ray equations independently of the AVO fit.
 h=np.asarray(public['h']);vp=np.asarray(public['vp']);bounds=np.r_[0,np.cumsum(2*h/vp[:3]),np.inf]
 # Basic operator invariants: exact synthetic two-term removal; constants
 # preserved; residuals retained; applying the transform twice is unchanged.
 I=np.sin(3*t);G=np.cos(7*t);synthetic=I[:,None]+G[:,None]*q
 result,_,_=remove_avo(synthetic,q);exact=float(np.max(np.abs(result-I[:,None])))
 assert exact<1e-12
 const=np.repeat(I[:,None],11,axis=1);assert np.max(np.abs(remove_avo(const,q)[0]-const))<1e-12
 ideal=zero_angle_truth(model,t)
 paths=sorted((RUN/'matlab/cases').glob(model['name']+'_*AVO_*.json'));assert len(paths)==10
 for path in paths:
  d=json.loads(path.read_text(encoding='utf-8'));meta=d['meta'];case=meta['id'];assert meta['domain']=='Offset'
  original={k:np.asarray(d[k]) for k in ['truth','fixed','optimized']};flat={};intercepts={};slopes={};idempotence=0
  for k,s in original.items():
   flat[k],intercepts[k],slopes[k]=remove_avo(s,q)
   idempotence=max(idempotence,float(np.max(np.abs(remove_avo(flat[k],q)[0]-flat[k]))))
  assert idempotence<1e-12
  # Identical linear processing of reference and both methods. The separate
  # analytic zero-angle target is shown to reveal transformation error.
  for method,key in [('Reference','truth'),('FixedSplit','fixed'),('PrunedFullCalRefit','optimized')]:
   rows.append(dict(CaseID=case,Model=model['name'],Scenario=meta['scenario'],Suite=meta['suite'],SNR_dB=meta['snr_db'],Seed=meta['seed'],Method=method,
    BeforeErrorToAVOTruth=float(np.linalg.norm(original[key]-original['truth'])/np.linalg.norm(original['truth'])),
    AfterErrorToTransformedReference=float(np.linalg.norm(flat[key]-flat['truth'])/np.linalg.norm(flat['truth'])),
    AfterErrorToIdealZeroAngle=float(np.linalg.norm(flat[key]-ideal)/np.linalg.norm(ideal))))
  np.savez_compressed(OUT/'data'/(case+'.npz'),t=t,q=q,ideal_zero_angle=ideal,truth=original['truth'],fixed=original['fixed'],optimized=original['optimized'],truth_noavo=flat['truth'],fixed_noavo=flat['fixed'],optimized_noavo=flat['optimized'],**{k+'_intercept':v for k,v in intercepts.items()},**{k+'_slope':v for k,v in slopes.items()})
  arrays=[original['truth'],original['fixed'],original['optimized'],flat['truth'],flat['fixed'],flat['optimized']]
  limit=max(float(np.abs(a).max()) for a in arrays)
  fig,axs=plt.subplots(2,3,figsize=(13,9),sharex=True,sharey=True,layout='constrained')
  titles=['纯净参考：保留 AVO','原分步法：保留 AVO','优化方法：保留 AVO','参考：拟合去 AVO','原分步法：拟合去 AVO','优化方法：拟合去 AVO']
  for ax,a,title in zip(axs.flat,arrays,titles):
   im=ax.imshow(a,extent=[-75,1575,t[-1],t[0]],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_title(title,fontsize=11);ax.set_xlabel('偏移距 / m');ax.set_xticks([0,750,1500])
  axs[0,0].set_ylabel('双程时间 / s');axs[1,0].set_ylabel('双程时间 / s');fig.colorbar(im,ax=axs,label='六幅图共用振幅色标',shrink=.85)
  fig.suptitle(f'{name} · {meta["scenario"]} · {meta["suite"]} · {meta["snr_db"]} dB · 种子 {meta["seed"]}',fontsize=13)
  fig.supxlabel('去 AVO = 减去拟合角度项，保留残差；不逐道归一化，不使用真实反射系数。时间轴已对齐。',fontsize=9)
  png=OUT/'figures'/(case+'.png');fig.savefig(png,dpi=140);plt.close(fig)
  selected=(meta['scenario']=='SmoothAVO' and meta['suite']=='normal' and meta['snr_db']==12 and meta['seed']==20262411)
  if selected:
   overviews.append((name,model['fpeak'],t,original['truth'],ideal,flat['truth']))
   # Predetermined waveform panels retain a single gain across all 6 gathers.
   fig,axs=plt.subplots(2,3,figsize=(12,9),sharex=True,sharey=True,layout='constrained');gain=.42/limit
   for ax,a,title in zip(axs.flat,arrays,titles):
    for j in range(11):
     line=j+1+gain*a[:,j];ax.plot(line,t,'k',lw=.5);ax.fill_betweenx(t,j+1,line,where=a[:,j]>0,color='k',alpha=.65)
    ax.set_title(title,fontsize=10);ax.set_xlim(0,12);ax.set_ylim(t[-1],t[0]);ax.set_xlabel('道号')
   axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s');fig.suptitle(name+' · 12 dB · 固定首个种子')
   fig.supxlabel(f'全部道与面板使用同一显示增益 {gain:.5g}；去 AVO 后保留拟合残差。',fontsize=9)
   fig.savefig(OUT/'figures'/(model['name']+'_wiggles.png'),dpi=150);fig.savefig(OUT/'figures'/(model['name']+'_wiggles.pdf'));plt.close(fig)
  records.append(dict(case=case,meta=meta,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),figure=png.relative_to(OUT).as_posix(),amplitude_limit=limit,figure_sha256=hashlib.sha256(png.read_bytes()).hexdigest()))
 basis_checks[model['name']]=dict(synthetic_two_term_max_error=exact,max_idempotence_error=idempotence,calibration_source=str(calpath),calibration_sha256=hashlib.sha256(calpath.read_bytes()).hexdigest())
 print('Completed',model['name'],flush=True)

fig,axs=plt.subplots(4,3,figsize=(13,15),layout='constrained',sharex=True,sharey=True)
limit=max(float(np.abs(a).max()) for _,_,_,a,b,c in overviews for a in [a,b,c])
for i,(name,freq,t,original,ideal,estimated) in enumerate(overviews):
 for j,a in enumerate([original,ideal,estimated]):
  im=axs[i,j].imshow(a,extent=[-75,1575,t[-1],t[0]],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit)
  if i==0:axs[i,j].set_title(['原始参考：保留 AVO','模型理想零角度参考','从数据拟合去 AVO'][j],fontsize=12)
  if j==0:axs[i,j].set_ylabel(name+f' · {freq} Hz\n双程时间 / s')
  axs[i,j].set_xticks([0,750,1500]);axs[i,j].set_xlabel('偏移距 / m')
fig.colorbar(im,ax=axs,label='全部参考面板共用振幅色标',shrink=.6)
fig.suptitle('消除 AVO 前后的参考道集对比',fontsize=17)
fig.supxlabel('中列由模型真值计算，仅作核对；右列只使用公开几何和参考道集拟合。薄层干涉与真实弱反射不会因去 AVO 自动消失。',fontsize=9)
fig.savefig(OUT/'reference_before_after.png',dpi=155);fig.savefig(OUT/'reference_before_after.pdf');plt.close(fig)

with (OUT/'metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(OUT/'figure_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'checks.json').write_text(json.dumps(dict(cases=len(records),model_checks=basis_checks,original_compensation_unchanged=True),ensure_ascii=False,indent=2),encoding='utf-8')
options=''.join(f'<option value="{html.escape(r["figure"])}">{html.escape(r["case"])}</option>' for r in records)
(OUT/'index.html').write_text('''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>消除AVO前后对比</title><style>body{font:16px system-ui;max-width:1400px;margin:24px auto;padding:0 20px}img{width:100%}select{max-width:100%;padding:10px}p{line-height:1.7}</style><h1>消除AVO前后对比</h1><p>40个实际已计算案例（4模型×2种AVO场景×5种噪声组合），使用同一两项AVO去除规则。保留拟合残差，没有把所有道直接替换成截距。原衰减补偿算法保持不变。</p><img src="reference_before_after.png" alt="参考道集保留AVO、理想零角度和拟合去AVO对比"><p><a href="METHOD_ZH.md">处理与局限说明</a> · <a href="metrics.csv">完整指标</a></p><select id="cases">'''+options+'''</select><img id="gather" alt="原分步法和优化方法去AVO前后"><script>const s=document.getElementById('cases'),im=document.getElementById('gather');s.onchange=()=>im.src=s.value;s.onchange();</script></html>''',encoding='utf-8')
(OUT/'METHOD_ZH.md').write_text('''# 消除AVO效应的对比

基于已经运行保存的40组偏移距道集，增加相同的AVO后处理；未重新估计Q，未修改两种补偿算法。断层、弱反射和尖灭Section场景的各道固定偏移距750米，不构成角度道集，因此不在横向位置之间拟合AVO。

## 实际处理

由独立校准后公开的速度、层厚和偏移距计算射线参数p(t,x)，构造q=p²/max_x(p²)。每个时间样点拟合

S(t,x) = I(t) + G(t)q(t,x) + e(t,x)

输出 S_noAVO = S - Gq = I + e。保留非AVO残差，不把各道替换成同一条拟合截距，不进行逐道振幅归一化，也不除以可能过零的反射系数。

同一几何与同一线性算子用于参考道集、原分步法和优化方法；拟合不读取真实Q、真实反射系数、Vs、密度或真实主频。0偏移距道q=0，因此该道保持原样；截距由全部11道共同估计。

## 参考与评分

reference_before_after的中列用真实模型参数计算Shuey截距A并与真实Ricker子波合成，是理想零角度参考，只用于核对。右列为实际拟合去AVO后的参考，并不保证与理想参考完全相同。参考事件内部采用固定事件角度，而后处理使用逐时间射线参数；波形重叠和薄层干涉也会造成近似误差。

metrics.csv分别给出：处理前对原始参考的误差；处理后对同样变换参考的误差；处理后对理想零角度参考的误差。三个指标的参考与分母不同，不能把处理前后误差比直接当作补偿改善比例。去除后拟合AVO梯度为零是操作目标，不能作为补偿算法保幅提升的证据。

所有40例均完整展示，包含先前深层低Q相干干扰失败案例。保留残差意味着干扰、噪声与补偿误差可能仍在；真实弱反射、薄层抵消、横向地质变化不会被这一操作可靠区分。两项AVO模型、大角度、横向变化及噪声均可能造成偏差。

全部六面板使用同一振幅色标，固定12dB首种子的四组波形图共用各案例显示增益。原数据和原补偿结果不覆盖。checks记录常量保持、已知两项AVO精确去除和重复处理幂等性检查；原始源码哈希另外核对。
''',encoding='utf-8')
# Preserve exact executable source with the deliverable.
import shutil
shutil.copy2(Path(__file__),OUT/'add_avo_comparison_source.py')
expected=json.loads((ROOT/'FROZEN_SOLVER_HASHES.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((ROOT/'core_p1'/p).read_bytes()).hexdigest()==h for p,h in expected.items())
from PIL import Image
for p in OUT.rglob('*.png'):
 with Image.open(p) as im:im.verify()
archive=ROOT.parent/'avo_comparison_20261023.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():z.write(p,Path(OUT.name)/p.relative_to(OUT))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps(dict(cases=len(records),output=str(OUT),archive=str(archive),all_frozen_solver_hashes_match=True),ensure_ascii=True),flush=True)
