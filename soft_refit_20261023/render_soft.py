from pathlib import Path
import hashlib,json,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from analyze_soft import analyze,ERRORS,METHODS

R=Path(__file__).resolve().parent;F=R/'figures'
summary=analyze(final=True)
assert json.loads((R/'API_CHECK.json').read_text())['maximum_relative_difference']<1e-12
assert int((R/'api_exit_code.txt').read_text(encoding='utf-8-sig'))==0
for name,h in json.loads((R/'PUBLIC_API_HASHES.json').read_text()).items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==h
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':12})
NAMES={'fixed':'最初分步法','unpruned':'未剪枝重拟合','hard':'旧硬剪枝','soft':'新软收缩'}
COLORS={'fixed':'#777777','unpruned':'#386cb0','hard':'#d78728','soft':'#078477'}
STYLES={'fixed':'--','unpruned':':','hard':'-.','soft':'-'}
lookup={(r['model'],r['condition'],r['Q_percent'],r['method']):r for r in summary['metrics']}
manifest=[];views={}

def save(fig,name,meta):
    path=F/(name+'.png');fig.savefig(path,dpi=400);plt.close(fig)
    with Image.open(path) as im:
        im.verify()
    with Image.open(path) as im:
        meta.update(width=im.width,height=im.height,dpi=im.info['dpi'][0])
    meta.update(file=name+'.png',sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    manifest.append(meta)
    return meta

for model in ['primary','complex']:
    ref=json.loads((R/'results'/model/'reference.json').read_text(encoding='utf-8'))
    x=np.array(ref['x']).ravel()/1000;t=np.array(ref['t']).ravel()
    cases={b:json.loads((R/'results'/model/f'Q_{b:+d}.json').read_text(encoding='utf-8')) for b in ERRORS}
    modelname='四层模型' if model=='primary' else '复杂模型'
    qlabel='相对真 Q' if model=='primary' else '相对既有标定 Q'
    for condition in ['Moderate','Strong']:
        data=next(z for z in ref['datasets'] if z['name']==condition)
        target=np.array(data['target']);obs=np.array(data['D'])
        records={b:next(z for z in v['records'] if z['name']==condition) for b,v in cases.items()}
        allarrays=[target,obs]+[np.array(records[b][m]) for b in ERRORS for m in METHODS]
        limit=float(np.quantile(np.concatenate([abs(a).ravel() for a in allarrays]),.995))
        absolute=max(float(abs(a).max()) for a in allarrays)
        for bias in ERRORS:
            arrays=[target]+[np.array(records[bias][m]) for m in ['unpruned','hard','soft']]
            names=['纯净参考']+[NAMES[m] for m in ['unpruned','hard','soft']]
            clips=[float(np.mean(abs(a)>limit)) for a in arrays]
            fig,axs=plt.subplots(1,4,figsize=(19,7.4),layout='constrained',sharex=True,sharey=True)
            for ax,a,label,clip in zip(axs,arrays,names,clips):
                im=ax.pcolormesh(x,t,a,shading='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
                err=np.linalg.norm(a-target)/np.linalg.norm(target)
                ax.set_title(f'{label}\n整体误差 {err:.3f} / 显示饱和 {clip*100:.2f}%',fontsize=12)
                ax.set_xlabel('横向位置 / km')
            axs[0].set_ylim(t[-1],t[0]);axs[0].set_ylabel('时间 / s')
            fig.colorbar(im,ax=axs,label='同场景全部 Q / 全方法共用细节色标',shrink=.84)
            fig.suptitle(f'{modelname} · 完整 99 道 · {condition} · {qlabel} {bias:+d}%',fontsize=18)
            fig.supxlabel('同一观测、Q、子波与噪声估计；共同 99.5% 振幅分位色标，仅影响显示；不逐道归一化。',fontsize=10)
            meta=save(fig,f'{model}_{condition}_Q{bias:+d}_detail',dict(model=model,condition=condition,Q_percent=bias,view='detail',traces=99,color_limit=limit,display_clipping_fractions=clips))
            views[f'{model}_{condition}_{bias}_detail']=meta

        # Unclipped overview of both pruning strategies for all Q settings.
        fig,axs=plt.subplots(2,6,figsize=(23,9),layout='constrained',sharex=True,sharey=True)
        panels=[(axs[0,0],target,'纯净参考'),(axs[1,0],obs,'同一观测')]
        for col,bias in enumerate(ERRORS,1):
            panels.extend([(axs[0,col],np.array(records[bias]['hard']),f'旧硬剪枝 / Q {bias:+d}%'),(axs[1,col],np.array(records[bias]['soft']),f'新软收缩 / Q {bias:+d}%')])
        for ax,a,label in panels:
            im=ax.pcolormesh(x,t,a,shading='nearest',cmap='RdBu_r',vmin=-absolute,vmax=absolute)
            ax.set_title(label);ax.set_xlabel('横向位置 / km')
        axs[0,0].set_ylim(t[-1],t[0]);axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s')
        fig.colorbar(im,ax=axs,label='同场景全部 Q / 全方法共同绝对最大振幅',shrink=.85)
        fig.suptitle(f'{modelname} / {condition} / {qlabel}：不截幅对比',fontsize=18)
        fig.supxlabel('所有面板同色标、无显示截幅；上排旧硬剪枝，下排新软收缩。请在原尺寸下查看。',fontsize=11)
        meta=save(fig,f'{model}_{condition}_absolute',dict(model=model,condition=condition,view='absolute',traces=99,color_limit=absolute,display_clipping=False))
        views[f'{model}_{condition}_absolute']=meta

        arrays=[target]+[np.array(records[0][m]) for m in ['unpruned','hard','soft']]
        names=['纯净参考']+[NAMES[m] for m in ['unpruned','hard','soft']]
        it=(t>=.98)&(t<=1.10) if model=='primary' else (t>=.95)&(t<=1.40)
        ix=np.ones(len(x),dtype=bool) if model=='primary' else (x>=.9)&(x<=2.4)
        zz=[a[np.ix_(it,ix)] for a in arrays];lim=max(float(abs(a).max()) for a in zz)
        fig,axs=plt.subplots(1,4,figsize=(19,7.4),layout='constrained',sharex=True,sharey=True)
        for ax,a,label in zip(axs,zz,names):
            im=ax.pcolormesh(x[ix],t[it],a,shading='nearest',cmap='RdBu_r',vmin=-lim,vmax=lim)
            ax.set_title(label);ax.set_xlabel('横向位置 / km')
        axs[0].set_ylim(t[it][-1],t[it][0]);axs[0].set_ylabel('时间 / s')
        fig.colorbar(im,ax=axs,label='局部各面板共同绝对最大振幅',shrink=.85)
        fig.suptitle(f'{modelname} / {condition} / Q 扰动 0%：固定局部放大',fontsize=18)
        fig.supxlabel('与上一轮相同物理窗口；不截幅、不逐道归一化。',fontsize=11)
        meta=save(fig,f'{model}_{condition}_zoom',dict(model=model,condition=condition,view='zoom',Q_percent=0,color_limit=lim,display_clipping=False))
        views[f'{model}_{condition}_zoom']=meta

        lim=max(float(abs(a).max()) for a in arrays);dx=float(np.median(np.diff(x)));gain=1.5*dx/lim
        fig,axs=plt.subplots(1,4,figsize=(19,8),layout='constrained',sharex=True,sharey=True)
        for ax,a,label in zip(axs,arrays,names):
            for j,xc in enumerate(x):
                xx=xc+gain*a[:,j];ax.plot(xx,t,color='#222222',lw=.35)
                ax.fill_betweenx(t,xc,xx,where=xx>=xc,color='#315f78',lw=0)
            ax.set_title(label);ax.set_xlabel('横向位置 / km');ax.set_xlim(x[0]-1.6*dx,x[-1]+1.6*dx)
        axs[0].set_ylim(t[-1],t[0]);axs[0].set_ylabel('时间 / s')
        fig.suptitle(f'{modelname} / {condition} / Q 扰动 0%：同增益 99 道波形',fontsize=18)
        fig.supxlabel('四面板同一全局增益；无 AGC、无逐道归一化。',fontsize=11)
        meta=save(fig,f'{model}_{condition}_wiggle',dict(model=model,condition=condition,view='wiggle',Q_percent=0,traces=99,global_gain_km_per_unit=gain))
        views[f'{model}_{condition}_wiggle']=meta

        aa=[target,np.array(records[0]['fixed']),np.array(records[0]['soft'])];lim=max(float(abs(a).max()) for a in aa)
        fig,axs=plt.subplots(1,3,figsize=(16,6.8),layout='constrained',sharex=True,sharey=True)
        for ax,a,label in zip(axs,aa,['纯净参考',NAMES['fixed'],NAMES['soft']]):
            im=ax.pcolormesh(x,t,a,shading='nearest',cmap='RdBu_r',vmin=-lim,vmax=lim)
            ax.set_title(label);ax.set_xlabel('横向位置 / km')
        axs[0].set_ylim(t[-1],t[0]);axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,shrink=.85,label='共同绝对最大振幅')
        fig.suptitle(f'{modelname} / {condition} / Q 扰动 0%：与最初分步法对比',fontsize=18)
        fig.supxlabel('同一观测与施加 Q / 子波；不截幅、不逐道归一化。',fontsize=10)
        meta=save(fig,f'{model}_{condition}_original',dict(model=model,condition=condition,view='original',Q_percent=0,traces=99,color_limit=lim,display_clipping=False))
        views[f'{model}_{condition}_original']=meta
        print('rendered',model,condition,flush=True)

fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
for row,model in enumerate(['primary','complex']):
    for col,condition in enumerate(['Moderate','Strong']):
        ax=axs[row,col]
        for m in METHODS:
            vals=[lookup[model,condition,b,m]['relative_error'] for b in ERRORS]
            ax.plot(ERRORS,vals,marker='o',ls=STYLES[m],color=COLORS[m],label=NAMES[m],lw=2)
        ax.axvline(0,color='#aaaaaa',lw=.8);ax.set_ylim(bottom=0);ax.grid(alpha=.2)
        ax.set_title(('四层模型' if model=='primary' else '复杂模型')+' / '+condition)
        ax.set_xlabel('相对真 Q 误差 / %' if model=='primary' else '相对既有标定 Q 扰动 / %');ax.set_ylabel('整体相对误差（越小越好）')
axs[0,0].legend(fontsize=10);fig.suptitle('完整 99 道：同一观测下的整体恢复误差',fontsize=17)
save(fig,'overall_metrics',dict(view='metric',metric='relative Frobenius error',replicates='One fixed observation per model/condition; no uncertainty interval'))
fig,axs=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
for ax,condition in zip(axs,['Moderate','Strong']):
    for m in METHODS:
        vals=[lookup['primary',condition,b,m]['weak_amplitude_error'] for b in ERRORS]
        ax.plot(ERRORS,vals,marker='o',ls=STYLES[m],color=COLORS[m],label=NAMES[m],lw=2)
    ax.axvline(0,color='#aaaaaa',lw=.8);ax.set_ylim(bottom=0);ax.grid(alpha=.2);ax.set_title(condition)
    ax.set_xlabel('相对真 Q 误差 / %');ax.set_ylabel('弱反射振幅相对误差（越小越好）')
axs[0].legend(fontsize=10);fig.suptitle('四层模型 4 个预定弱反射：全部 99 道',fontsize=17)
save(fig,'weak_metrics',dict(view='metric',metric='Ricker template projection amplitude relative error',events=4,traces=99))
(R/'PLOT_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')

def fmt(z):return f'{z:.3f}'
table=['| 模型 / 噪声 | 最初分步法 | 未剪枝 | 旧硬剪枝 | 新软收缩 |','|---|---:|---:|---:|---:|']
for p in ['primary','complex']:
    for c in ['Moderate','Strong']:
        table.append('| '+('四层' if p=='primary' else '复杂')+' / '+c+' | '+' | '.join(fmt(lookup[p,c,0,m]['relative_error']) for m in METHODS)+' |')
weak=['| 噪声 | 最初分步法 | 未剪枝 | 旧硬剪枝 | 新软收缩 |','|---|---:|---:|---:|---:|']
for c in ['Moderate','Strong']:
    weak.append('| '+c+' | '+' | '.join(fmt(lookup['primary',c,0,m]['weak_amplitude_error']) for m in METHODS)+' |')
report='''## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: run
- Origin Date: 2026-10-23
- Verification Status: UNVERIFIED
- Version Label: soft_group_refit_v1

## 本轮改动

从原 FullCalRefit 出发，对原有时间支撑组作连续收缩，不再像旧方案一样整组硬删除。所有模型和 Q 设置使用相同规则，没有读取模型名称、真反射位置或目标振幅来设权重。

组系数为 a，源域子波 Gram 矩阵为 M，固定支撑的岭回归映射为 K。令 E=a' M a，V=sigma² tr(M K K')，使用 gain=E/(E+V) 缩放该组系数。零传播方差时 gain=1。这里 V 是固定已选支撑下的白噪声传播代理；未包含支撑选择、相关噪声和 Q 偏差的不确定性，不能当作完整后验。

## 完整 99 道回测

两模型、两档噪声、五种 Q 设置，共 20 组。使用上一轮完全相同的观测、噪声估计、子波、Q 和连续 11 道批次。重算的未剪枝输出与旧缓存相对差异为 0；180 个分块均完成，优化异常数为 0。独立直接矩阵公式、Monte Carlo 方差传播、零噪声和死道边界检查通过。公开整幅道集入口在两个模型的 99 道强噪声 Q0 案例中与批量结果一致。这是数值与接口核对，未进行新的独立模型或噪声实现确认，故保留 UNVERIFIED 的 ARS 研究验证标签。

这些都是已知模型上的开发回测，20 组设置不是 20 个独立样本。不据此宣称统计显著或普遍泛化。

## Q 扰动 0% 的整体相对误差

小值较好。四层模型的 0% 为真 Q；复杂模型的 0% 为既有标定 Q，二者含义不同。

'''+ '\n'.join(table)+'''

## 四层模型 Q0 的弱反射振幅相对误差

同上一轮，以四个预定弱事件的 Ricker 模板投影计算。理想误差为 0。该误差不等于“丢失振幅百分比”。CSV 另有弱反射振幅回归增益，理想值为 1。

'''+ '\n'.join(weak)+f'''

新软收缩整体误差低于最初方法 {summary['overall_wins']['fixed']}/20 组，低于未剪枝 {summary['overall_wins']['unpruned']}/20 组，低于旧硬剪枝 {summary['overall_wins']['hard']}/20 组；弱反射误差低于硬剪枝 {summary['weak_wins']['hard']}/10 组。

复杂模型的 10 组整体误差均优于旧硬剪枝；四层模型 10 组整体误差均高于硬剪枝，但弱反射均改善。软收缩是本轮保弱反射与抑噪之间更温和的折中，保留旧硬剪枝作强抑噪对照。

Q 低估问题未完全解决：复杂模型 Moderate、−30% 时软收缩整体误差约 0.950，仍高于最初方法的 0.829。全部指标见 analysis/metrics.csv。

## 看图与运行

打开 index.html，默认显示复杂模型 Strong、Q 扰动 0%。细节图为 7600×2960 像素、400 DPI；点击“原尺寸”可滚动。所有同场景 Q 与方法共用 99.5% 振幅分位色标，图中标注显示饱和比例。另附共同绝对最大振幅、不截幅的硬/软对比总览；波形图为共同全局增益，不逐道归一化。

MATLAB R2024b 中运行 demo_soft，或调用 out=soft_refit_compensate(D,cm,offset,sigma)。D 为时间×道，cm 为已有标定配置，offset 为每道偏移距或一个公共值，sigma 为从观测估计的各道噪声尺度。sigma 省略时沿用原 0.02–0.25 s 窗口 MAD；若该窗口含信号应传入适合资料的噪声估计。公开函数按 11 道分块限制内存，不需要输入纯净目标。返回的 initial_normal_residual 只描述软收缩之前的初始岭重拟合，收缩后不宣称仍满足原正规方程。

下载包包含新代码、完整输入/输出、检查点及上游冻结结果。解压后两个目录保持同级。运行 run_soft_trial 会使用 cache 中已有检查点；如需重新求解，先在目录副本中将 cache 改名并建立空 cache 文件夹。render_soft.py 依赖 NumPy、Matplotlib、Pillow；analyze_soft.py --final 检查输出与冻结哈希。prepare_soft_trial.py 仅用于原工作区准备，不是解压后的运行前提。
'''
(R/'RESULTS_ZH.md').write_text(report,encoding='utf-8')
html='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>软收缩：99道验证</title><style>body{margin:0;font:17px/1.6 sans-serif;color:#16354b}header{padding:20px 28px;background:#edf5f6}h1{margin:0;font-size:28px}nav{position:sticky;top:0;background:white;z-index:2;padding:12px 24px;display:flex;gap:12px;align-items:center;flex-wrap:wrap;border-bottom:1px solid #c8dce2}select,button{font:inherit;padding:5px 8px}a{color:#007d85}#stage{overflow:auto;width:100%;background:#e7eef0}#plot{display:block;width:100%;height:auto;cursor:zoom-in}#stage.native #plot{width:auto;max-width:none;cursor:zoom-out}#info,footer{padding:12px 28px}.metrics{padding:20px 28px}.metrics img{width:100%;max-width:1500px;display:block}</style>
<header><h1>连续软收缩 · 完整 99 道</h1><div>同一观测、Q、子波与噪声估计；两模型共享一条规则。复杂模型整体误差与弱反射保留的取舍见下方指标。</div><a href="RESULTS_ZH.md">方法、结果与运行说明</a> · <a href="analysis/metrics.csv">全部数值 CSV</a></header>
<nav><label>模型 <select id="model"><option value="complex">复杂模型</option><option value="primary">四层模型</option></select></label><label>噪声 <select id="condition"><option>Strong</option><option>Moderate</option></select></label><label>Q <select id="bias"><option>-30</option><option>-10</option><option selected>0</option><option>10</option><option>30</option></select>%</label><label>视图 <select id="view"><option value="detail">剪枝消融对比</option><option value="zoom">固定局部放大（Q0）</option><option value="wiggle">99道波形（Q0）</option><option value="original">与最初方法对比（Q0）</option><option value="absolute">全Q不截幅硬/软总览</option></select></label><button id="size">原尺寸</button><a id="original" target="_blank">打开 PNG</a></nav><div id="info"></div><div id="stage"><img id="plot" alt="99道同色标方法对比"></div>
<div class="metrics"><h2>整体误差</h2><img src="figures/overall_metrics.png" alt="两模型、两档噪声、五种Q设置下四种方法的整体相对误差"><h2>弱反射误差</h2><img src="figures/weak_metrics.png" alt="四层模型四个弱事件的振幅相对误差"></div><footer>细节图共用同场景全部 Q / 四方法 / 参考 / 观测合并 99.5% 振幅分位色标，仅影响显示，未修改数据。所有曲线来自同一批开发观测，没有重复采样置信区间。</footer>
<script>const views=VIEWS;const $=id=>document.getElementById(id);function update(){let m=$('model').value,c=$('condition').value,b=$('bias').value,v=$('view').value,a=views[m+'_'+c+(v==='detail'?'_'+b+'_detail':'_'+v)];$('plot').src='figures/'+a.file;$('original').href='figures/'+a.file;$('bias').disabled=v!=='detail';$('info').textContent=a.width+' × '+a.height+' 像素，400 DPI。'+(m==='primary'?'Q 相对真值。':'Q 相对既有标定值。')+(v==='detail'?'面板：纯净参考、未剪枝、旧硬剪枝、新软收缩；显示饱和比例标于各面板。':'所有面板同一显示标度。')}for(let id of ['model','condition','bias','view'])$(id).onchange=update;function toggle(){let on=$('stage').classList.toggle('native');$('size').textContent=on?'适应窗口':'原尺寸'}$('size').onclick=toggle;$('plot').onclick=toggle;update();</script></html>'''.replace('VIEWS',json.dumps(views,ensure_ascii=False))
(R/'index.html').write_text(html,encoding='utf-8')
(R/'GALLERY_VIEWS.json').write_text(json.dumps(views,ensure_ascii=False,indent=2),encoding='utf-8')
print('Figures, report and gallery complete',len(manifest),flush=True)
