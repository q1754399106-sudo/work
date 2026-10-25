from pathlib import Path
import csv,hashlib,json,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
R=Path(__file__).resolve().parent;F=R/'figures'
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':12,'pdf.fonttype':42})
checks=json.loads((R/'checks.json').read_text(encoding='utf-8'));assert checks['complete'] and checks['blocks']==180 and checks['optimization_flags']==0
assert int((R/'matlab_exit_code.txt').read_text(encoding='utf-8-sig'))==0
assert all(hashlib.sha256((R/'source'/p).read_bytes()).hexdigest()==h for p,h in json.loads((R/'FROZEN_CORE_HASHES.json').read_text(encoding='utf-8')).items())
errors=[-30,-10,0,10,30];manifest=[];metrics=[];views={};primary=json.loads((R/'inputs/primary.json').read_text(encoding='utf-8'))
def save(fig,name,record,pdf=True):
    path=F/(name+'.png');fig.savefig(path,dpi=400)
    if pdf:fig.savefig(F/(name+'.pdf'),dpi=400)
    plt.close(fig)
    with Image.open(path) as im:im.verify()
    with Image.open(path) as im:record.update(width_pixels=im.width,height_pixels=im.height)
    record.update(file=name+'.png',pdf_file=name+'.pdf' if pdf else None,dpi=400,sha256=hashlib.sha256(path.read_bytes()).hexdigest());manifest.append(record)
def weak_error(out):
    rr=np.array(primary['reflectivity']);est=[];true=[];dt=primary['config']['dt'];freq=primary['config']['fpeak']
    for j,events in enumerate(primary['event_indices_zero_based']):
        for e in [1,3,5,7]:
            k=events[e];ii=np.arange(max(0,k-10),min(len(out),k+11));z=np.pi*freq*(ii-k)*dt;w=(1-2*z*z)*np.exp(-z*z);est.append(np.dot(out[ii,j],w)/np.dot(w,w));true.append(rr[k,j])
    return float(np.linalg.norm(np.array(est)-true)/np.linalg.norm(true))

for protocol in ['primary','complex']:
    ref=json.loads((R/'results'/protocol/'reference.json').read_text(encoding='utf-8'));t=np.array(ref['t']).ravel();x=np.array(ref['x']).ravel()/1000
    assert len(x)==99
    results={e:json.loads((R/'results'/protocol/f'Q_{e:+d}.json').read_text(encoding='utf-8')) for e in errors}
    for condition in ['Moderate','Strong']:
        data=next(d for d in ref['datasets'] if d['name']==condition);target=np.array(data['target']);obs=np.array(data['D']);assert target.shape[1]==99
        rr={e:next(d for d in z['records'] if d['name']==condition) for e,z in results.items()}
        arrays=[target,obs]+[np.array(rr[e][m]) for e in errors for m in ['fixed','optimized']]
        absolute_limit=max(float(abs(a).max()) for a in arrays)
        detail_limit=float(np.quantile(np.concatenate([abs(a).ravel() for a in arrays]),.995))
        title_model='四层模型 99 道' if protocol=='primary' else '原复杂模型完整 99 道'
        qlabel='真 Q 误差' if protocol=='primary' else '相对标定 Q 扰动'
        live=np.linalg.norm(obs,axis=0)>0
        for error in errors:
            for method in ['fixed','unpruned','optimized']:
                a=np.array(rr[error][method]);assert a.shape==target.shape and np.isfinite(a).all()
                metrics.append(dict(protocol=protocol,condition=condition,error_percent=error,method=method,traces=99,relative_error=float(np.linalg.norm(a-target)/np.linalg.norm(target)),live_error=float(np.linalg.norm((a-target)[:,live])/np.linalg.norm(target[:,live])),weak_amplitude_error=weak_error(a) if protocol=='primary' else None))
            outputs=[target,np.array(rr[error]['fixed']),np.array(rr[error]['optimized'])];clips=[float(np.mean(abs(a)>detail_limit)) for a in outputs]
            fig,axs=plt.subplots(1,3,figsize=(16,6.8),layout='constrained',sharex=True,sharey=True)
            for ax,a,label,clip in zip(axs,outputs,['纯净参考','最初分步法（同 Q / 子波）','优化剪枝法'],clips):
                im=ax.pcolormesh(x,t,a,shading='nearest',cmap='RdBu_r',vmin=-detail_limit,vmax=detail_limit,rasterized=True)
                err=float(np.linalg.norm(a-target)/np.linalg.norm(target));ax.set_title(f'{label}\n误差 {err:.3f} / 显示饱和 {100*clip:.2f}%',fontsize=12);ax.set_xlabel('横向位置 / km')
            axs[0].set_ylim(t[-1],t[0]);axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='所有 Q / 两方法共用细节色标',shrink=.86)
            fig.suptitle(f'{title_model} / {condition} / {qlabel} {error:+d}%',fontsize=17)
            caption='共同色标为同场景全部 Q、参考与观测合并绝对幅度 99.5% 分位；仅显示饱和，不改数据，不逐道归一化。'
            if protocol=='complex':caption+=' 保留原第 23、78 死道。'
            fig.supxlabel(caption,fontsize=9)
            name=f'{protocol}_{condition}_Q{error:+d}_detail'
            record=dict(protocol=protocol,condition=condition,Q_error_percent=error,view='detail',traces=99,color_limit=detail_limit,color_scope='All five Q settings, reference and observation within this protocol/condition',display_clip_fractions=clips,individual_trace_gain=False)
            save(fig,name,record);views[f'{protocol}_{condition}_{error}_detail']=record

        # Unclipped overview keeps all physical amplitudes visible at a shared scale.
        fig,axs=plt.subplots(2,6,figsize=(23,9),layout='constrained',sharex=True,sharey=True)
        panels=[(axs[0,0],target,'纯净参考'),(axs[1,0],obs,'同一观测')]
        for col,error in enumerate(errors,1):
            panels.extend([(axs[0,col],np.array(rr[error]['fixed']),f'分步法 / Q {error:+d}%'),(axs[1,col],np.array(rr[error]['optimized']),f'优化法 / Q {error:+d}%')])
        for ax,a,title in panels:
            im=ax.pcolormesh(x,t,a,shading='nearest',cmap='RdBu_r',vmin=-absolute_limit,vmax=absolute_limit,rasterized=True);ax.set_title(title,fontsize=12);ax.set_xlabel('横向位置 / km')
        axs[0,0].set_ylim(t[-1],t[0]);axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='共同绝对最大振幅色标',shrink=.85)
        fig.suptitle(f'{title_model} / {condition} / {qlabel}：不截幅总览',fontsize=18);fig.supxlabel('原始振幅同一色标，不截幅、不逐道归一化。上排分步法，下排优化剪枝法。请打开原尺寸查看细节。',fontsize=11)
        record=dict(protocol=protocol,condition=condition,view='absolute',traces=99,color_limit=absolute_limit,amplitude_clipping=False,individual_trace_gain=False)
        save(fig,f'{protocol}_{condition}_absolute',record,pdf=False);views[f'{protocol}_{condition}_absolute']=record

        # 99 true wiggle traces, one common global gain; no trace normalization.
        outputs=[target,np.array(rr[0]['fixed']),np.array(rr[0]['optimized'])];wiggle_lim=max(float(abs(a).max()) for a in outputs);dx=float(np.median(np.diff(x)));gain=1.5*dx/wiggle_lim
        fig,axs=plt.subplots(1,3,figsize=(16,7.4),layout='constrained',sharex=True,sharey=True)
        for ax,a,label in zip(axs,outputs,['纯净参考','最初分步法','优化剪枝法']):
            for j,xc in enumerate(x):
                xx=xc+gain*a[:,j];ax.plot(xx,t,color='#222222',lw=.35);ax.fill_betweenx(t,xc,xx,where=xx>=xc,color='#315f78',lw=0)
            ax.set_title(label);ax.set_xlabel('横向位置 / km');ax.set_xlim(x[0]-1.6*dx,x[-1]+1.6*dx)
        axs[0].set_ylim(t[-1],t[0]);axs[0].set_ylabel('时间 / s');fig.suptitle(f'{title_model} / {condition} / {qlabel} 0%：同增益波形道集',fontsize=17)
        fig.supxlabel('实际绘出 99 道，所有面板同一全局振幅显示倍数；无 AGC、无逐道归一化。',fontsize=10)
        record=dict(protocol=protocol,condition=condition,view='wiggle',Q_error_percent=0,traces=99,global_gain_km_per_unit=gain,individual_trace_gain=False)
        save(fig,f'{protocol}_{condition}_wiggle_0',record);views[f'{protocol}_{condition}_wiggle']=record

        # A predeclared local window, shown only for the zero perturbation case.
        it=(t>=.98)&(t<=1.10) if protocol=='primary' else (t>=.95)&(t<=1.40)
        ix=np.ones(len(x),dtype=bool) if protocol=='primary' else (x>=.9)&(x<=2.4)
        xx=x[ix];tt=t[it];zoom=[a[np.ix_(it,ix)] for a in outputs];zlim=max(float(abs(a).max()) for a in zoom)
        fig,axs=plt.subplots(1,3,figsize=(16,6.8),layout='constrained',sharex=True,sharey=True)
        for ax,a,label in zip(axs,zoom,['纯净参考','最初分步法','优化剪枝法']):
            im=ax.pcolormesh(xx,tt,a,shading='nearest',cmap='RdBu_r',vmin=-zlim,vmax=zlim,rasterized=True);ax.set_title(label);ax.set_xlabel('横向位置 / km')
        axs[0].set_ylim(tt[-1],tt[0]);axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='局部三面板共同绝对最大振幅',shrink=.85)
        fig.suptitle(f'{title_model} / {condition} / Q 0%：固定局部放大',fontsize=17);fig.supxlabel('固定物理时间窗与横向范围；仅改变展示范围。局部各面板色标一致，不截幅、不逐道归一化。',fontsize=10)
        record=dict(protocol=protocol,condition=condition,view='zoom',Q_error_percent=0,traces_shown=int(ix.sum()),time_window_s=[float(tt[0]),float(tt[-1])],x_window_km=[float(xx[0]),float(xx[-1])],color_limit=zlim,amplitude_clipping=False,individual_trace_gain=False)
        save(fig,f'{protocol}_{condition}_zoom_0',record,pdf=False);views[f'{protocol}_{condition}_zoom']=record
        print(protocol,condition,'rendered',flush=True)

with (R/'analysis/metrics_99.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
(R/'analysis/summary.json').write_text(json.dumps(dict(checks=checks,metrics=metrics),ensure_ascii=False,indent=2),encoding='utf-8')
(R/'PLOT_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
html='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>99 道高清 Q 对比</title>
<style>body{margin:0;font:17px/1.65 sans-serif;color:#16354b}header{padding:20px 28px;background:#eef5f7}h1{margin:0;font-size:28px}nav{position:sticky;top:0;z-index:3;background:white;padding:12px 24px;border-bottom:1px solid #cadce4;display:flex;gap:12px;align-items:center;flex-wrap:wrap}button,select{font:inherit;padding:5px 10px}#stage{overflow:auto;width:100%;background:#e8edef}#plot{display:block;width:100%;height:auto;cursor:zoom-in}#stage.native #plot{width:auto;max-width:none;cursor:zoom-out}#info{padding:12px 28px}.links{padding:10px 28px}a{color:#087f8c}</style>
<header><h1>完整 99 道：高清 Q 误差对比</h1><div>四层模型为重新合成的 99 道；复杂模型为原有完整 99 道，未用插值补道。细节图为 400 DPI；PDF 的坐标文字为矢量，密度部分按 400 DPI 渲染。</div></header>
<nav><label>模型 <select id="model"><option value="complex">原复杂模型</option><option value="primary">四层模型</option></select></label><label>噪声 <select id="condition"><option>Strong</option><option>Moderate</option></select></label><label>Q 偏差 <select id="bias"><option>-30</option><option>-10</option><option selected>0</option><option>10</option><option>30</option></select>%</label><label>视图 <select id="view"><option value="detail">同色标细节图</option><option value="absolute">不截幅全 Q 总览</option><option value="wiggle">99 道波形图（Q 0%）</option><option value="zoom">固定局部放大（Q 0%）</option></select></label><button id="size">显示原尺寸</button><a id="original" target="_blank">打开原图</a><a id="pdf" target="_blank">下载 PDF</a></nav>
<div id="info"></div><div id="stage"><img id="plot" alt="完整99道对比"></div><div class="links"><a href="figures/11_vs_99_reference.png" target="_blank">11 道 / 99 道密度差异</a> | <a href="README_ZH.md">说明与计算范围</a> | <a href="analysis/metrics_99.csv">完整99道指标</a></div>
<script>const views=VIEWS;const $=id=>document.getElementById(id);function update(){let m=$('model').value,c=$('condition').value,b=$('bias').value,v=$('view').value;let key=m+'_'+c+(v==='detail'?'_'+b+'_detail':'_'+v);let a=views[key];$('plot').src='figures/'+a.file;$('original').href='figures/'+a.file;$('pdf').style.display=a.pdf_file?'inline':'none';if(a.pdf_file)$('pdf').href='figures/'+a.pdf_file;$('info').textContent=a.width_pixels+' × '+a.height_pixels+' 像素；'+(m==='primary'?'Q 百分比相对真值':'Q 百分比相对既有标定值')+'。'+(v==='detail'?'细节图采用全部Q结果共同99.5%分位色标，饱和比例标于图上；需检查幅度极值时切换不截幅总览。':'各面板统一显示标度，无逐道归一化。');$('bias').disabled=v!=='detail'}for(let id of ['model','condition','bias','view'])$(id).onchange=update;function toggle(){let on=$('stage').classList.toggle('native');$('size').textContent=on?'适应窗口':'显示原尺寸'}$('size').onclick=toggle;$('plot').onclick=toggle;update();</script></html>'''.replace('VIEWS',json.dumps(views,ensure_ascii=False))
(R/'index.html').write_text(html,encoding='utf-8')
note='''# 完整 99 道高清道集

本轮针对道数过少和网页缩放导致的图像可读性问题。所有 Q 设置为 −30%、−10%、0%、+10%、+30%。原核心方法、正则化、子波和噪声尺度规则不变。

- 四层模型：在原 0-300 m 横向范围上按同一连续事件规律重新采样 99 道并生成固定噪声，未插值原来的补偿结果。因此是新观测，不应直接套用旧 11 道误差数值。
- 复杂模型：使用原模型已有完整 99 道，原观测逐值保留；第 23、78 死道仍保留。所有五个 Q 情况都重新计算，采用连续 11 道分块。
- 四层模型的 Q 偏差相对真 Q；复杂模型相对既有四层标定 Q，不能混为相对空间真值的误差。
- 分步法和优化方法共享所施加 Q 与子波，65 个核心 MATLAB 文件保持哈希不变。本轮未处理真实资料。

## 看图

打开 index.html。每次显示一个 Q 情况的三面板，避免把十几个小图挤在同一屏。点击“显示原尺寸”或图像可在原像素尺寸下滚动查看；“打开原图”单独显示 PNG。

细节 PNG 为 6400 × 2720 像素、400 DPI。所有同场景 Q / 两方法结果采用共同的 99.5% 绝对振幅分位色标；每个面板标明显示饱和比例。仅影响显示，不改变数据或指标。另有共同绝对最大幅度、不截幅的全 Q 总览，供幅度极值检查。

波形道集使用一个全局增益，未逐道归一化。PDF 保留矢量坐标文字，密度数据按 400 DPI 渲染；波形曲线为矢量。局部放大使用预定物理时间 / 空间范围，不对新方法单独增强。

## 数据与复现

本地实验目录的 inputs/ 保存完整输入，results/ 保存各 Q 结果与诊断，analysis/metrics_99.csv 为新 99 道指标。下载的 dense_gathers_hd_20261023.zip 为图件包，仅含图件、页面、指标、元数据及脚本，不含计算输入、输出和 MATLAB 核心文件；复现需使用本地完整实验目录。原始 11 道实验和结论保留，不能把新、旧不同采样下的绝对指标直接混用。

在目录副本中先将已有 results/primary 和 results/complex 改名存档，再在 MATLAB R2024b 中调用 run_dense。入口拒绝覆盖已有计算。绘图脚本 render_dense.py 需要 NumPy、Matplotlib、Pillow；首次准备脚本 prepare_dense.py 引用工作区原复杂模型。
'''
(R/'README_ZH.md').write_text(note,encoding='utf-8')
dest=R.parent/'dense_gathers_hd_20261023.zip'
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for name in ['index.html','README_ZH.md','PLAN.json','PLOT_MANIFEST.json','render_dense.py','prepare_dense.py','run_dense.m','FROZEN_CORE_HASHES.json','INPUT_PROVENANCE.json','checks.json']:
        z.write(R/name,'q_dense_gathers_20261023/'+name)
    for folder in ['figures','analysis']:
        for p in (R/folder).glob('*'):
            if p.is_file():z.write(p,'q_dense_gathers_20261023/'+str(p.relative_to(R)).replace('\\','/'))
with zipfile.ZipFile(dest) as z:
    assert z.testzip() is None
    for a in manifest:assert hashlib.sha256(z.read('q_dense_gathers_20261023/figures/'+a['file'])).hexdigest()==a['sha256']
delivery=dict(png_figures=len(manifest),computed_traces=99,Q_cases=5,conditions=2,models=2,zip=str(dest),bytes=dest.stat().st_size,all_figures_verified=True,minimum_detail_pixels=[6400,2720])
(R/'DELIVERY_CHECK.json').write_text(json.dumps(delivery,ensure_ascii=False,indent=2),encoding='utf-8')
# Leave earlier results intact, but point both old pages to the new complete gallery.
for old,link in [(R.parent/'q_error_stability_20261023/index.html','../q_dense_gathers_20261023/index.html'),(R.parent/'q_error_stability_20261023/gather_gallery/index.html','../../q_dense_gathers_20261023/index.html')]:
    text=old.read_text(encoding='utf-8');banner=f'<p style="padding:16px;background:#e1eff2;font-size:22px"><a href="{link}">已更新：点击查看完整 99 道高清版（400 DPI / 原尺寸查看）</a></p>'
    if link not in text:text=text.replace('<h1>',banner+'<h1>',1);old.write_text(text,encoding='utf-8')
print(json.dumps(delivery,ensure_ascii=True))
