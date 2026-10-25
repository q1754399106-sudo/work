from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fk_baseline import soft_fk
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'field_data';FIG=ROOT/'figures'
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False})
records={}
for name in ['wz25','alaska_111','alaska_120','alaska_130']:
 d=np.load(OUT/(name+'_raw.npz'));D=d['data'];dt=float(d['dt']);x=d['x'];t=d['t'];filtered,noise,diag=soft_fk(D,dt,float(np.median(np.diff(x))))
 np.savez_compressed(OUT/(name+'_fk.npz'),raw=D,filtered=filtered,removed=noise,t=t,x=x,dt=dt)
 earlier=OUT/(name+'_processed.npz');adaptive=np.load(earlier)['denoised'] if earlier.exists() else D
 shown=[a*t[:,None]**2 for a in [D,adaptive,filtered,noise]];limit=float(np.quantile(np.abs(np.concatenate([a.ravel() for a in shown])),.995))
 fig,axs=plt.subplots(1,4,figsize=(17,8),sharex=True,sharey=True,layout='constrained')
 for ax,a,title in zip(axs,shown,['真实输入','长孔径分频减除','常规软 f–k 基线','f–k 减掉的分量']):
  im=ax.imshow(a,extent=[x[0]/1000,x[-1]/1000,t[-1],t[0]],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_title(title,fontsize=11);ax.set_xlabel('偏移距 / km')
 axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='共同 t² 显示增益后的振幅',shrink=.8);fig.suptitle(name+'：实测资料两种去噪尝试',fontsize=16)
 fig.supxlabel('同一t²显示增益与共同99.5%分位色标；不进行逐道归一化。去除能量不等于去除噪声，右图需检查真实反射泄漏。',fontsize=9)
 fig.savefig(FIG/(name+'_field_comparison.png'),dpi=150);fig.savefig(FIG/(name+'_field_comparison.pdf'));plt.close(fig)
 records[name]=dict(**diag,adaptive_change_fraction=float(np.linalg.norm(adaptive-D)/np.linalg.norm(D)),display_limit=limit,display_clip_fractions=[float(np.mean(np.abs(a)>limit)) for a in shown])
(OUT/'fk_diagnostics.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(records,ensure_ascii=True,indent=2))
