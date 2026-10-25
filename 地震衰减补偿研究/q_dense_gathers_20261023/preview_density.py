from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;d=json.loads((R/'inputs/complex.json').read_text(encoding='utf-8'));t=np.array(d['public']['t']);x=np.array(d['x'])/1000
target=np.array(next(a for a in d['datasets'] if a['name']=='Strong')['target']);idx=np.array([1,11,21,30,40,50,60,70,79,89,99])-1
lim=float(np.quantile(abs(target),.995));plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':12})
fig,axs=plt.subplots(1,2,figsize=(14,7),layout='constrained',sharey=True)
for ax,xx,a,title in [(axs[0],x[idx],target[:,idx],'之前：只抽取 11 道'),(axs[1],x,target,'现在：原模型完整 99 道')]:
    im=ax.pcolormesh(xx,t,a,shading='nearest',cmap='RdBu_r',vmin=-lim,vmax=lim);ax.set_title(title);ax.set_xlabel('横向位置 / km')
axs[0].set_ylim(t[-1],t[0]);axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='两面板共同振幅色标',shrink=.85)
fig.suptitle('相同复杂模型的纯净参考：道数造成的横向细节差异',fontsize=17)
fig.supxlabel('不是图像插值。99 道是原模型已有数据；统一采用完整参考绝对振幅 99.5% 分位色标，超过范围仅显示饱和。',fontsize=10)
fig.savefig(R/'figures/11_vs_99_reference.png',dpi=300);plt.close(fig)
