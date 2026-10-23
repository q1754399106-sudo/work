"""General research figures; fixed cases, no filtering or waveform normalization."""
from pathlib import Path
import argparse,csv,hashlib,json,platform,statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
METHODS=['Frozen','CombinedUncertainty','FrozenJointCalibration','CombinedJointCalibration','CombinedUncertaintyV2']
COLORS=['#666666','#B35806','#6A3D9A','#007A45','#0066A1'];MARKERS=['s','^','D','x','o']
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();run=a.run.resolve();plan=read(run/'FROZEN_EXECUTION_PLAN.json')
 assert read(run/'analysis/coverage.json')['complete'];dest=run/'analysis/figures';dest.mkdir(exist_ok=False)
 with (run/'metrics.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
 cohorts=[('prospective_synthetic','新合成模型 / Prospective synthetic',False),('original_nine_regression','原九模型回归 / Original regression',True)]
 manifest=[]
 def export(fig,name,caption,description):
  for fmt in ['png','pdf']:
   path=dest/(name+'.'+fmt);fig.savefig(path,dpi=300,facecolor='white');item=dict(file=path.name,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
   if fmt=='png':
    with Image.open(path) as im:item.update(pixels=list(im.size),dpi=list(im.info.get('dpi',[])),mode=im.mode)
   manifest.append(item)
  (dest/(name+'.txt')).write_text(caption+'\n\n'+description+'\n',encoding='utf-8');plt.close(fig)
 with plt.rc_context({'font.family':'sans-serif','font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.transparent':False}):
  fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained',sharex=True)
  table=[]
  for ax,(cohort,title,reuse) in zip(axes,cohorts):
   names=[j['model']['name'] for j in plan['jobs'] if j['reuse_original']==reuse]
   for im,(method,color,marker) in enumerate(zip(METHODS,COLORS,MARKERS)):
    values=[]
    for name in names:
     rr=[r for r in rows if r['Model']==name and r['Suite']=='Main' and r['Route']=='Split' and r['Method']==method]
     value=statistics.mean(float(r['RelativeError']) for r in rr);values.append(value);table.append(dict(cohort=cohort,model=name,method=method,main_mean_relative_error=value,observations=len(rr)))
    ax.scatter(values,np.arange(len(names))+(im-2)*.11,label=method,c=color,marker=marker,s=32)
   ax.set(yticks=np.arange(len(names)),yticklabels=names,title=title,xlabel='主试验平均相对误差 / Mean relative error');ax.set_xlim(left=0);ax.grid(axis='x',alpha=.2);ax.invert_yaxis()
  handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8)
  export(fig,'main_all_models','图1：所有模型、五种方法的 Main Split 平均误差，越低越好。每个点平均三次校准×六观测种子×五固定场景，不是独立重复次数为90的置信区间。','两面板使用共同横轴，含全部模型和五种对照，不按效果排序。点形与颜色同时区分方法；原始逐观测数据见 metrics.csv，绘图汇总见 main_plot_data.csv。')
  with (dest/'main_plot_data.csv').open('w',encoding='utf-8-sig',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
  boot=read(run/'analysis/crossed_bootstrap.json');fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained',sharex=True)
  for ax,(cohort,title,reuse) in zip(axes,cohorts):
   for i,(baseline,z) in enumerate(boot[cohort].items()):
    point=100*z['relative_improvement'];lo,hi=[100*x for x in z['crossed_percentile_95_interval']]
    ax.plot([lo,hi],[i,i],color='#0066A1',lw=2);ax.scatter([point],[i],marker='o',color='#0066A1',s=32)
   ax.set(yticks=range(4),yticklabels=list(boot[cohort]),title=title,xlabel='V2 相对改善 / Improvement (%)');ax.axvline(0,color='black',ls='--',lw=.8);ax.grid(axis='x',alpha=.2);ax.invert_yaxis()
  export(fig,'crossed_uncertainty','图2：V2 相对四种对照的 Main 等模型平均误差改善率及交叉 bootstrap 95% 百分位区间。正数表示改善。','分别重采样模型、校准和观测种子，组内方法配对及共同随机数保持一致；场景固定。两组分别包含4或9模型、3校准、6观测种子。区间没有宣称真实实测泛化或分层Q后验覆盖。完整数值见 crossed_bootstrap.json。')
  names=[j['model']['name'] for j in plan['jobs'] if not j['reuse_original']]
  fig,axes=plt.subplots(len(names),2,figsize=(13,2.5*len(names)),layout='constrained',sharex=True)
  for i,name in enumerate(names):
   examples=[read(run/name/'matlab'/('example_'+scenario+'.json')) for scenario in ['WeakSection','CoherentNoiseAVO']]
   limit=max(abs(float(v)) for ex in examples for key in ['truth','observed','CombinedUncertainty','CombinedUncertaintyV2'] for v in ex[key])*1.05
   for ax,ex in zip(axes[i],examples):
    for key,label,color,style in [('observed','观测 / Observed','#888888','-'),('truth','真值 / Truth','#000000','-'),('CombinedUncertainty','原版 / Original','#B35806','--'),('CombinedUncertaintyV2','改进 / V2','#0066A1',':')]:
     ax.plot(ex['t'],ex[key],label=label,color=color,ls=style,lw=1.1)
    ax.set(title=name+' / '+ex['scenario'],ylim=(-limit,limit),xlabel='时间 / Time (s)',ylabel='振幅 / Amplitude');ax.grid(alpha=.15)
  handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=4)
  export(fig,'fixed_waveform_examples','图3：四个新合成模型的全部固定示例：Main、第一校准、第一观测种子、第6道；WeakSection 和 CoherentNoiseAVO。','原始振幅，无平滑、重标定或逐道归一化；每个模型两场景共享纵轴范围，显示完整时间记录。所有案例按预先封存模型顺序展示，包含原方法和V2的误差。校准与反演不读取图中的真值。原始数组位于各模型 example_*.json。')
 (dest/'FIGURE_MANIFEST.json').write_text(json.dumps(dict(files=manifest,python=platform.python_version(),matplotlib=matplotlib.__version__,numpy=np.__version__,source_metrics_sha256=hashlib.sha256((run/'metrics.csv').read_bytes()).hexdigest(),publisher='General research report; no journal compliance claimed',manual_visual_review='PENDING',transformations=['Main arithmetic mean by model and method','Previously computed crossed bootstrap intervals','Fixed raw trace arrays without smoothing or amplitude normalization']),ensure_ascii=False,indent=2),encoding='utf-8')
 print(str(dest))
if __name__=='__main__':main()
