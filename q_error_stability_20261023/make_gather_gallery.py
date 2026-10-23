from pathlib import Path
import csv,hashlib,json,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;OUT=R/'gather_gallery';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
summary=json.loads((R/'analysis/summary.json').read_text(encoding='utf-8'));rows=summary['metrics']
inputs={p:json.loads((R/'inputs'/f'{p}.json').read_text(encoding='utf-8')) for p in ['primary','complex']}
entries=[];manifest=[]
def read(protocol,condition,error):
    folder=R/f'matlab_{protocol}'
    for p in folder.glob('*.json'):
        if p.name in ['reference.json','progress.json','checks.json']:continue
        a=json.loads(p.read_text(encoding='utf-8'))
        if a['scenario']['error_percent']==error:return next(r for r in a['records'] if r['name']==condition)
    raise KeyError(error)
def metric(protocol,condition,method,error):return next(a for a in rows if a['protocol']==protocol and a['condition']==condition and a['method']==method and a['uniform'] and a['Q_error_percent']==error)['relative_error']
for protocol in ['primary','complex']:
    ref=json.loads((R/f'matlab_{protocol}/reference.json').read_text(encoding='utf-8'));t=np.array(ref['t']).ravel();x=np.array(inputs[protocol]['x'])/1000
    for condition in ['Moderate','Strong']:
        original=next(r for r in ref['datasets'] if r['name']==condition);target=np.array(original['target']);observed=np.array(original['D'])
        # Use a moderate-bias display to avoid a single -50% blow-up washing out every panel.
        errors=[-30,-10,0,10,30] if protocol=='primary' else [-30,-10,0,10,30]
        records={e:read(protocol,condition,e) for e in errors}
        arrays=[target,observed]+[np.array(records[e][m]) for e in errors for m in ['fixed','optimized']]
        limit=max(float(abs(a).max()) for a in arrays)
        fig,axs=plt.subplots(2,6,figsize=(22,9),layout='constrained',sharex=True,sharey=True)
        panels=[(axs[0,0],target,'纯净参考'),(axs[1,0],observed,'同一含噪观测')]
        for col,e in enumerate(errors,1):
            for ri,method,label in [(0,'fixed','最初分步法'),(1,'optimized','优化剪枝法')]:
                title=f'{label}\nQ {e:+d}% / 误差 {metric(protocol,condition,method,e):.3f}'
                panels.append((axs[ri,col],np.array(records[e][method]),title))
        for ax,a,title in panels:
            im=ax.pcolormesh(x,t,a,shading='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True);ax.set_title(title,fontsize=10);ax.set_xlabel('横向位置 / km')
        axs[0,0].set_ylim(t[-1],t[0]);axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s')
        fig.colorbar(im,ax=axs,label='所有面板共同原幅度色标',shrink=.85)
        subtitle='四层模型：真 Q=[40,60,85,110]' if protocol=='primary' else '原复杂模型：固定 11 道 / 相对标定 Q 扰动'
        fig.suptitle(f'{subtitle} / {condition} / 上：分步法，下：优化剪枝法',fontsize=17)
        fig.supxlabel('共同绝对最大幅度色标，不截幅、不逐道归一化。所有 Q 情况为同一输入；标题中的误差为整体相对波形误差。',fontsize=10)
        name=f'{protocol}_{condition}_density';fig.savefig(OUT/(name+'.png'),dpi=170);plt.close(fig)
        entries.append((f'{subtitle} / {condition}：密度道集',name))
        manifest.append(dict(file=name+'.png',protocol=protocol,condition=condition,type='density',Q_errors_percent=errors,common_amplitude_limit=limit,individual_trace_gain=False,amplitude_clipping=False))

        # Native wiggles with ONE global amplitude-to-position conversion for all panels.
        fig,axs=plt.subplots(2,5,figsize=(19,9),layout='constrained',sharex=True,sharey=True)
        spacing=float(np.median(np.diff(x)));gain=.40*spacing/limit
        for col,e in enumerate(errors):
            for ri,method,label in [(0,'fixed','最初分步法'),(1,'optimized','优化剪枝法')]:
                ax=axs[ri,col];a=np.array(records[e][method])
                for j,xc in enumerate(x):
                    xx=xc+gain*a[:,j];ax.plot(xx,t,color='#222222',lw=.55);ax.fill_betweenx(t,xc,xx,where=xx>=xc,color='#315f78',lw=0,alpha=.9)
                ax.set_title(f'{label} / Q {e:+d}%',fontsize=10);ax.set_xlim(x[0]-.55*spacing,x[-1]+.55*spacing);ax.set_xlabel('横向位置 / km')
        axs[0,0].set_ylim(t[-1],t[0]);axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s')
        fig.suptitle(f'{subtitle} / {condition}：同增益波形道集',fontsize=17)
        fig.supxlabel('全部面板使用同一个全局振幅显示倍数；无逐道归一化、无 AGC。不同横向位置为不同道，非一个炮集的偏移距轴。',fontsize=10)
        name=f'{protocol}_{condition}_wiggle';fig.savefig(OUT/(name+'.png'),dpi=170);plt.close(fig)
        entries.append((f'{subtitle} / {condition}：波形道集',name));manifest.append(dict(file=name+'.png',protocol=protocol,condition=condition,type='wiggle',Q_errors_percent=errors,common_amplitude_position_gain_km_per_unit=gain,individual_trace_gain=False))

        if protocol=='primary':
            # Predeclared weak event near 1.04 s; no search for the most flattering zoom.
            take=(t>=.98)&(t<=1.10);tz=t[take];ez=[-10,0,10];zz={e:records[e] for e in ez}
            pieces=[target[take],observed[take]]+[np.array(zz[e][m])[take] for e in ez for m in ['fixed','optimized']];lim=max(float(abs(a).max()) for a in pieces)
            fig,axs=plt.subplots(2,4,figsize=(15,7),layout='constrained',sharex=True,sharey=True)
            panels=[(axs[0,0],target[take],'纯净参考'),(axs[1,0],observed[take],'含噪观测')]
            for col,e in enumerate(ez,1):
                panels.extend([(axs[0,col],np.array(zz[e]['fixed'])[take],f'分步法 / Q {e:+d}%'),(axs[1,col],np.array(zz[e]['optimized'])[take],f'优化剪枝 / Q {e:+d}%')])
            for ax,a,title in panels:
                im=ax.pcolormesh(x,tz,a,shading='nearest',cmap='RdBu_r',vmin=-lim,vmax=lim,rasterized=True);ax.set_title(title);ax.set_xlabel('横向位置 / km')
            axs[0,0].set_ylim(tz[-1],tz[0]);axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='局部八面板共同原幅度色标',shrink=.85)
            fig.suptitle(f'四层模型 / {condition}：1.04 s 附近弱反射放大',fontsize=16)
            fig.supxlabel('时间窗固定为 0.98-1.10 s；仅改变展示范围，数据和各面板振幅标度一致，不单独增强新方法。',fontsize=9)
            name=f'primary_{condition}_weak_zoom';fig.savefig(OUT/(name+'.png'),dpi=170);plt.close(fig)
            entries.append((f'四层模型 / {condition}：弱反射局部放大',name));manifest.append(dict(file=name+'.png',protocol=protocol,condition=condition,type='weak_zoom',time_window_s=[.98,1.10],Q_errors_percent=ez,common_amplitude_limit=lim,individual_trace_gain=False))

html='<!doctype html><meta charset="utf-8"><title>Q误差道集对比</title><style>body{max-width:1900px;margin:25px auto;font:18px sans-serif;line-height:1.8;color:#16354b}img{width:100%}a{color:#087f8c}nav{position:sticky;top:0;background:#fff;padding:12px}</style><h1>不同 Q 误差下的合成道集</h1><p>两种方法均为原反射模型版本，共享 Q 和子波。四层模型的横轴为相对真 Q 误差，复杂模型则为相对原标定 Q 的扰动。所有图用同一观测。新方法通常有较低误差，但严重低估 Q 时不保证更稳定。</p><nav><a href="#complex_Strong_density">复杂模型强干扰</a> | <a href="#primary_Strong_density">四层模型强噪声</a> | <a href="#primary_Strong_weak_zoom">弱反射放大</a> | <a href="../index.html">完整稳定性报告</a></nav>'
for title,name in entries:html+=f'<h2 id="{name}">{title}</h2><p><a href="{name}.png">打开原尺寸 PNG</a></p><img src="{name}.png" alt="{title}">'
(OUT/'index.html').write_text(html,encoding='utf-8')
(OUT/'PLOT_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
dest=R.parent/'q_error_gathers_20261023.zip'
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in OUT.iterdir():
        if p.is_file():z.write(p,'q_error_gathers/'+p.name)
with zipfile.ZipFile(dest) as z:assert z.testzip() is None
print(json.dumps(dict(figures=len(manifest),gallery=str(OUT/'index.html'),zip=str(dest),bytes=dest.stat().st_size),ensure_ascii=True))
