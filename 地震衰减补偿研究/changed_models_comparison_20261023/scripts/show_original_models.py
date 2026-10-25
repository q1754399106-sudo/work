from pathlib import Path
import json,csv,hashlib,zipfile,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'original_models';OUT.mkdir(exist_ok=True)
(OUT/'data').mkdir(exist_ok=True)
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
plan=json.loads((ROOT/'PLAN.json').read_text(encoding='utf-8'));models=plan['models']
run=Path(json.loads((ROOT/'run.json').read_text(encoding='utf-8'))['run'])
names=['深层低 Q','薄层极性反转','低阻抗差','高频 Q 反转']
scenes=plan['scenarios'];scene_names=['正常 AVO','相干干扰场景的纯净参考','断层剖面','弱反射剖面','尖灭剖面']
truths={};sources=[]
for m in models:
 for scene in scenes:
  path=run/'matlab/cases'/f'{m["name"]}_{scene}_normal_18dB_20262411.json'
  d=json.loads(path.read_text(encoding='utf-8'));t=np.asarray(d['t']).reshape(-1);a=np.asarray(d['truth']);assert a.shape==(601,11)
  truths[m['name'],scene]=(t,a)
  np.savetxt(OUT/'data'/f'{m["name"]}_{scene}.csv',np.column_stack([t,a]),delimiter=',',header='time_s,'+','.join('trace_'+str(j) for j in range(1,12)),comments='',fmt='%.17g')
  np.savez_compressed(OUT/'data'/f'{m["name"]}_{scene}.npz',t=t,truth=a)
  sources.append(dict(model=m['name'],scene=scene,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),field='truth'))
shutil.copy2(ROOT/'PLAN.json',OUT/'model_parameters.json')
with (OUT/'layer_parameters.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.writer(f);w.writerow(['model','layer','vp_m_s','vs_m_s','rho_g_cm3','Q','top_depth_m','thickness_m'])
 for m in models:
  h=np.diff([0]+m['event_t'])*np.array(m['vp'][:3])/2;z=np.r_[0,np.cumsum(h)]
  for j in range(4):w.writerow([m['name'],j+1,m['vp'][j],m['vs'][j],m['rho'][j],m['Q'][j],z[j],h[j] if j<3 else 'halfspace'])

# Base models are the one-dimensional layered media actually used by the forward code.
fig=plt.figure(figsize=(16,9),layout='constrained');gs=fig.add_gridspec(2,5,height_ratios=[3.3,1.25],width_ratios=[1,1,1,1,.045]);norm=Normalize(20,130);cmap=plt.colormaps['viridis'];depth_max=1300
for i,(m,name) in enumerate(zip(models,names)):
 ax=fig.add_subplot(gs[0,i]);h=np.diff([0]+m['event_t'])*np.array(m['vp'][:3])/2;z=np.r_[0,np.cumsum(h),depth_max]
 for j in range(4):
  ax.add_patch(Rectangle((0,z[j]),500,z[j+1]-z[j],facecolor=cmap(norm(m['Q'][j])),edgecolor='white',lw=.8))
  if z[j+1]-z[j]>70:ax.text(250,(z[j]+z[j+1])/2,f'第 {j+1} 层\nQ = {m["Q"][j]}',ha='center',va='center',color='white' if m['Q'][j]<90 else 'black',fontsize=11)
 for v in z[1:-1]:ax.axhline(v,color='black',lw=.6)
 if h[2]<60:
  ax.annotate(f'第 3 层：{h[2]:.1f} m',xy=(400,(z[2]+z[3])/2),xytext=(440,z[2]+95),ha='right',color='black',bbox=dict(facecolor='white',alpha=.92,edgecolor='none'),arrowprops=dict(arrowstyle='->',color='white',lw=1.3))
 ax.set_xlim(0,500);ax.set_ylim(depth_max,0);ax.set_xticks([0,250,500]);ax.set_xlabel('横向距离示意 / m');ax.set_title(f'{name}\n主频 {m["fpeak"]} Hz',fontsize=12)
 if i==0:ax.set_ylabel('深度 / m')
 tabax=fig.add_subplot(gs[1,i]);tabax.axis('off')
 cells=[[str(j+1),str(m['vp'][j]),str(m['vs'][j]),f'{m["rho"][j]:.2f}',str(m['Q'][j]),f'{h[j]:.1f}' if j<3 else '半空间'] for j in range(4)]
 table=tabax.table(cellText=cells,colLabels=['层','Vp','Vs','ρ','Q','厚度'],loc='center',cellLoc='center',colWidths=[.09,.20,.20,.16,.13,.22]);table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.7)
 tabax.set_title('Vp、Vs：m/s；ρ：g/cm³；厚度：m',fontsize=8,pad=2)
fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_subplot(gs[0,4]),label='真实 Q 值（全部模型共用色标）')
fig.suptitle('本轮实验的四个基础地层模型',fontsize=17)
fig.supxlabel('界面深度按 h = Δt × Vp / 2 计算；第四层为半空间，图中仅画至 1300 m。横向均匀，图中宽度用于示意。',fontsize=10)
fig.savefig(OUT/'base_layer_models.png',dpi=180);fig.savefig(OUT/'base_layer_models.pdf');plt.close(fig)

# Pure reference gathers: read exact stored truth, no regeneration or gain.
limit=max(float(np.abs(truths[m['name'],'SmoothAVO'][1]).max()) for m in models)
reference_limit=limit
fig,axs=plt.subplots(2,2,figsize=(12,10),layout='constrained')
for ax,m,name in zip(axs.flat,models,names):
 t,a=truths[m['name'],'SmoothAVO'];im=ax.imshow(a,extent=[-75,1575,t[-1],t[0]],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,interpolation='nearest')
 ax.set_title(f'{name} · {m["fpeak"]} Hz');ax.set_xlabel('偏移距 / m');ax.set_ylabel('双程时间 / s');ax.set_xticks([0,750,1500])
fig.colorbar(im,ax=axs,label='参考振幅（四个模型共用色标）',shrink=.85)
fig.suptitle('模型原始参考道集：未衰减、未加噪、未补偿',fontsize=16)
fig.supxlabel('11 道，偏移距 0–1500 m；已按反射时间对齐。低阻抗差模型本身振幅较弱，不做逐道归一化。',fontsize=10)
fig.savefig(OUT/'clean_reference_gathers.png',dpi=170);fig.savefig(OUT/'clean_reference_gathers.pdf');plt.close(fig)

# Show every exact scene reference separately, with one amplitude scale per model.
for m,name in zip(models,names):
 limit=max(float(np.abs(truths[m['name'],s][1]).max()) for s in scenes)
 fig,axs=plt.subplots(1,5,figsize=(16,6),sharey=True,layout='constrained')
 for ax,s,title in zip(axs,scenes,scene_names):
  t,a=truths[m['name'],s];extent=[-75,1575,t[-1],t[0]] if s.endswith('AVO') else [-25,525,t[-1],t[0]]
  im=ax.imshow(a,extent=extent,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,interpolation='nearest');ax.set_title(title,fontsize=10);ax.set_xlabel('偏移距 / m' if s.endswith('AVO') else '横向位置 / m')
 axs[0].set_ylabel('双程时间 / s');fig.colorbar(im,ax=axs,label='本模型各场景共用振幅色标',shrink=.9)
 fig.suptitle(name+'：五种场景的原始参考道集',fontsize=15)
 fig.supxlabel('纯净参考不含相干干扰：相干噪声仅加入观测数据，因此前两幅纯净参考相同。断层/尖灭等是参考反射事件的几何变化。',fontsize=9)
 fig.savefig(OUT/(m['name']+'_all_scenes.png'),dpi=160);fig.savefig(OUT/(m['name']+'_all_scenes.pdf'));plt.close(fig)

notes='''# 四个模型本来的样子

- base_layer_models.png / pdf：本轮实际使用的四层基础介质、真实Q、速度、密度和层厚。
- clean_reference_gathers.png / pdf：四个模型正常AVO的精确参考道集，未衰减、未加噪、未补偿，统一色标。
- Changed_*_all_scenes.png / pdf：每模型五种场景的精确纯净参考，模型内统一色标。
- layer_parameters.csv：完整分层参数；model_parameters.json：原实验参数计划。
- data/：20个模型/场景组合的参考数据，CSV和NPZ，均直接取自已经运行保存的truth字段；CSV第一列为秒，随后为11道。

基础介质为一维四层模型，界面深度由实验采用的 h=Δt×Vp/2 得到。第四层为半空间，地层图只截画至1300米；横向宽度用于示意。
断层、弱反射和尖灭场景在参考反射事件的时间或系数上作变化，再通过一维衰减算子产生观测；不是含断层或尖灭体的完整二维波动方程正演。
参考道集采用Shuey两项AVO与Ricker子波，已按反射时间对齐；不包含几何扩散、透射、多次波或NMO拉伸。
WeakSection第二反射系数乘0.08；FaultSection有8毫秒正弦起伏，第7至11道的深部两反射延迟40毫秒；PinchoutSection第二、第三反射异号并随横向靠拢，最后三道完全抵消。
薄层极性反转模型WeakSection中央道的第二、第三反射重合，不能分别恢复振幅。
'''
(OUT/'README_ZH.md').write_text(notes,encoding='utf-8')
(OUT/'provenance.json').write_text(json.dumps(dict(source_field='truth',processing='None; exact source arrays exported. Colors/gain described in figures.',sources=sources),ensure_ascii=False,indent=2),encoding='utf-8')
archive=ROOT.parent/'original_models_20261023.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():z.write(p,Path(OUT.name)/p.relative_to(OUT))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps(dict(folder=str(OUT),zip=str(archive),reference_arrays=len(truths),common_reference_amplitude=reference_limit),ensure_ascii=True))
