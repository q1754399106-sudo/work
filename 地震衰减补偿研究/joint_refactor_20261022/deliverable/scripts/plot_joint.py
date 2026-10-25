"""General research figures; all planned configurations, no statistical CI."""
from pathlib import Path
import csv,json,hashlib,statistics,platform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
METHODS=['Legacy','DecoupledWeak','Decoupled','SharedIsotropic','Physical']
COLORS=['#555555','#0072B2','#A35000','#007A45','#AA3377']
MARKERS=['s','o','^','D','x']
def main():
 run=Path(json.loads((ROOT/'LATEST_RUN.json').read_text(encoding='utf-8'))['run'])
 rows=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')));assert len(rows)==640
 dest=ROOT/'figures';dest.mkdir(exist_ok=True);manifest=[];plotdata=[]
 models=list(dict.fromkeys(r['Model'] for r in rows))
 def export(fig,name,caption):
  for fmt in ('png','pdf'):
   p=dest/(name+'.'+fmt);fig.savefig(p,dpi=240,facecolor='white');item={'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
   if fmt=='png':
    with Image.open(p) as im:item|={'pixels':list(im.size),'dpi':list(im.info['dpi'])}
   manifest.append(item)
  (dest/(name+'.txt')).write_text(caption,encoding='utf-8');plt.close(fig)
 with plt.rc_context({'font.family':'sans-serif','font.sans-serif':['Microsoft YaHei','DejaVu Sans'],
     'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.transparent':False}):
  fig,axes=plt.subplots(1,2,figsize=(12,5.3),layout='constrained',sharex=True,sharey=True)
  for ax,cal in zip(axes,['Old','New']):
   for ki,(method,color,marker) in enumerate(zip(METHODS,COLORS,MARKERS)):
    values=[]
    for model in models:
     rs=[r for r in rows if r['Calibration']==cal and r['Route']=='Split' and r['Method']==method and r['Model']==model]
     value=statistics.mean(float(r['RelativeError']) for r in rs);values.append(value)
     plotdata.append({'Calibration':cal,'Model':model,'Method':method,'MeanRelativeError':value,'n':len(rs)})
    ax.scatter(values,np.arange(len(models))+(ki-2)*.13,c=color,marker=marker,label=method,s=42)
   ax.set(title=('旧校准 / Old calibration' if cal=='Old' else '新校准 / New calibration'),xlabel='平均相对误差 / Mean relative error',yticks=range(len(models)),yticklabels=models)
   ax.grid(axis='x',alpha=.2);ax.set_xlim(left=0)
  axes[0].invert_yaxis();fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,fontsize=9)
  export(fig,'calibration_regularization','分步路线全部开发模型与五方法。每个点为四固定场景×两个观测种子的平均，不是独立重复8次的置信区间。两面板共享坐标，低值更好；不删去不利模型。底层数据见 model_plot_data.csv 和逐行 metrics.csv。')
  fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
  for ax,metric,label in zip(axes.flat,['AVOInterceptError','AVOGradientError','WeakEventError','PinchoutTailError'],['AVO截距 / AVO intercept','AVO梯度 / AVO gradient','弱反射 / Weak event','尖灭尾部 / Pinchout tail']):
   count=0
   for ki,(method,color,marker) in enumerate(zip(METHODS[1:],COLORS[1:],MARKERS[1:])):
    rs=[r for r in rows if r['Calibration']=='New' and r['Route']=='Split' and r['Method']==method and np.isfinite(float(r[metric]))]
    ratios=[]
    for r in rs:
     baseline=next(b for b in rows if b['Calibration']=='New' and b['Route']=='Split' and b['Method']=='Legacy' and all(b[k]==r[k] for k in ['Model','Scenario','Seed']))
     ratios.append(float(r[metric])/float(baseline[metric]))
    assert all(x>0 and np.isfinite(x) for x in ratios)
    ax.scatter(ratios,ki+np.linspace(-.18,.18,len(ratios)),color=color,marker=marker,s=25,alpha=.8)
    ax.scatter([statistics.mean(ratios)],[ki],marker='|',s=210,color='black',linewidths=2)
    count=len(ratios)
   ax.axvline(1,color='black',ls='--',lw=1);ax.set_xscale('log',base=2)
   ax.set(title=f'{label} (n={count}/method)',yticks=range(4),yticklabels=METHODS[1:],xlabel='局部误差 / 同观测 Legacy 误差 (log2 scale)');ax.grid(axis='x',alpha=.2)
  export(fig,'local_paired_ratios','新校准分步路线，所有模型/场景/种子的局部误差配对比值。小于1为改善；大于1为退化。黑竖线是配对比值的算术均值，彩色点为原始配对；垂直偏移只防重叠，不表示数值。AVO每方法16点，弱反射/尖灭各8点；场景和模型共享噪声种子，点并非相互独立。横轴为以2为底对数，未删除点，无显著性检验。')
 with (dest/'model_plot_data.csv').open('w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(plotdata[0]));w.writeheader();w.writerows(plotdata)
 (dest/'manifest.json').write_text(json.dumps({'destination':'general research report, no journal target','source':str(run/'matlab/metrics.csv'),
  'source_sha256':hashlib.sha256((run/'matlab/metrics.csv').read_bytes()).hexdigest(),
  'python':platform.python_version(),'matplotlib':matplotlib.__version__,'numpy':np.__version__,'files':manifest},indent=2),encoding='utf-8')
 print(dest)
if __name__=='__main__':main()
