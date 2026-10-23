from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiband_filter import suppress_multiband
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'field_data';FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False})
diagnostics={}
for name in ['wz25','alaska_111','alaska_120','alaska_130']:
 p=OUT/(name+'_raw.npz');d=np.load(p);D=d['data'];t=d['t'];x=d['x'];dt=float(d['dt']);filtered,removed,diag=suppress_multiband(D,dt,x)
 np.savez_compressed(OUT/(name+'_processed.npz'),data=D,denoised=filtered,removed=removed,t=t,x=x,dt=dt)
 diagnostics[name]=dict(**diag,raw_npz_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),real_data_has_clean_truth=False,processing='Denoising prototype only; no Q compensation, no AGC or t^2 gain in processing',display_gain='t^2 common to original/output/removed; common symmetric 99.5-percentile amplitude limit; clipped display fraction recorded')
 gain=t*t;arrays=[D,filtered,removed];shown=[a*gain[:,None] for a in arrays];limit=float(np.quantile(np.abs(np.concatenate([a.ravel() for a in shown])),.995));diagnostics[name]['display_amplitude_limit']=limit;diagnostics[name]['display_clipped_fractions']=[float(np.mean(np.abs(a)>limit)) for a in shown]
 fig,axs=plt.subplots(1,3,figsize=(15,8),sharex=True,sharey=True,layout='constrained')
 for ax,a,title in zip(axs,shown,['真实观测','分频相干减除后','被减掉的分量（需检查有效波泄漏）']):
  im=ax.imshow(a,extent=[x[0]/1000,x[-1]/1000,t[-1],t[0]],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_title(title,fontsize=11);ax.set_xlabel('偏移距 / km')
 axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='共用振幅色标（统一 t² 显示增益）',shrink=.85)
 fig.suptitle(name+'：真实炮集初步去干扰试验',fontsize=16)
 fig.supxlabel('t²仅用于显示，处理保持原始振幅单位；色标为三面板合并绝对振幅99.5%分位。无纯净真值，不据此声称恢复精度。',fontsize=9)
 fig.savefig(FIG/(name+'_field_denoise.png'),dpi=155);fig.savefig(FIG/(name+'_field_denoise.pdf'));plt.close(fig)
 print('Field processed',name,'removed fraction',diag['removed_norm_fraction'],flush=True)
(OUT/'experiment_diagnostics.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2),encoding='utf-8')
