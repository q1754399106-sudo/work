from pathlib import Path
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
checks=json.loads((R/'matlab/checks.json').read_text(encoding='utf-8'));assert checks['complete'] and checks['cases']==36 and checks['optimization_flags']==0
assert int((R/'matlab_exit_code.txt').read_text(encoding='utf-8-sig'))==0
hashes=json.loads((R/'FROZEN_CORE_HASHES.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((R/'source'/name).read_bytes()).hexdigest()==h for name,h in hashes.items())
datasets=[];controls=[];field=[];sensitivity=[]
methods=['input','fixed','refit','optimized'];labels=['输入','两步岭回归适配版','稀疏重拟合适配版','稀疏剪枝适配版'];colors=['#777777','#d47832','#547daf','#087f8c']
for path in sorted((R/'matlab').glob('20*.json')):
    d=json.loads(path.read_text(encoding='utf-8'));datasets.append(d);t=np.array(d['t']).ravel();mask=(t>=-1)&(t<=7);nominal=next(a for a in d['results'] if a['kind']=='full_path' and a['factor']==1)
    for a in d['results']:
        if a['kind']=='controlled_extra':
            target=np.array(a['target']).ravel()
            for method in methods:
                out=np.array(a[method]).ravel();den=np.linalg.norm(target[mask]);err=np.linalg.norm((out-target)[mask])/den
                controls.append(dict(station=d['metadata']['station'],event=d['metadata']['event'],split=d['metadata']['split'],extra_tstar_s=a['tstar'],added_noise_factor=a['noise_factor'],method=method,relative_error=float(err),correlation=float(np.corrcoef(out[mask],target[mask])[0,1]),amplitude_ratio=float(np.linalg.norm(out[mask])/den)))
        else:
            inp=np.array(a['input']).ravel();bg=np.array(a['background_input']).ravel()
            for method in methods[1:]:
                out=np.array(a[method]).ravel();replay=np.array(a[method+'_replay']).ravel();noise=np.array(a[method+'_background']).ravel()
                center=np.array(nominal[method]).ravel()
                sensitivity.append(dict(station=d['metadata']['station'],event=d['metadata']['event'],split=d['metadata']['split'],method=method,tstar_factor=a['factor'],tstar_s=a['tstar'],amplitude_gain=float(np.linalg.norm(out[mask])/np.linalg.norm(inp[mask])),background_gain=float(np.linalg.norm(noise[mask])/max(np.linalg.norm(bg[mask]),1e-30)),replay_relative_residual=float(np.linalg.norm((replay-inp)[mask])/np.linalg.norm(inp[mask])),relative_change_from_nominal=float(np.linalg.norm((out-center)[mask])/max(np.linalg.norm(center[mask]),1e-30))))
                if a['factor']==1:field.append(sensitivity[-1])
    # Every output is restored to the input physical amplitude scale for display.
    arrays=[np.array(nominal[m]).ravel()*d['metadata']['scale'] for m in methods]
    fig,axs=plt.subplots(3,1,figsize=(12,9),layout='constrained')
    for i in [0,1,3]:axs[0].plot(t,arrays[i],color=colors[i],lw=.8,label=labels[i])
    axs[0].set(xlim=(-1,7),xlabel='相对 P 到时 / s',ylabel='原文件幅度',title='同幅度标定的真实波形；没有无衰减真值');axs[0].legend(fontsize=9)
    f=np.fft.rfftfreq(1024,.05);taper=np.hanning(len(t))
    for i in [0,1,3]:axs[1].semilogy(f,abs(np.fft.rfft(arrays[i]*taper,n=1024)),color=colors[i],label=labels[i])
    axs[1].set(xlim=(.2,8),xlabel='频率 / Hz',ylabel='全窗 Hann 振幅谱',title='频谱变宽并不自动等于有效信息恢复');axs[1].legend(fontsize=9)
    for a in [a for a in d['results'] if a['kind']=='full_path']:
        axs[2].plot(t,np.array(a['optimized']).ravel()*d['metadata']['scale'],lw=.75,label=f"t* × {a['factor']:.2f}")
    axs[2].set(xlim=(-1,7),xlabel='相对 P 到时 / s',ylabel='原文件幅度',title='剪枝适配版对 t* ±25% 的响应');axs[2].legend(fontsize=9)
    fig.suptitle(f"{d['metadata']['event']} / {d['metadata']['station']}：作者 Q={d['metadata']['Q']:.2f}，t*={d['metadata']['tstar']:.5f} s")
    fig.supxlabel('W=I 的未知波形稀疏适配；保留核心两步岭回归 / 重拟合 / 剪枝代码。不是原分层反射模型的直接实测验证。',fontsize=9)
    fig.savefig(R/'figures'/(d['metadata']['station']+'_field_trial.png'),dpi=150);fig.savefig(R/'figures'/(d['metadata']['station']+'_field_trial.pdf'));plt.close(fig)

for name,rows in [('controlled_metrics.csv',controls),('field_metrics.csv',field),('sensitivity_metrics.csv',sensitivity)]:
    with (R/'analysis'/name).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
grouped=[]
for split in ['development','evaluation']:
    for noise in [0,1]:
        for method in methods:
            z=[a['relative_error'] for a in controls if a['split']==split and a['added_noise_factor']==noise and a['method']==method]
            grouped.append(dict(split=split,noise_factor=noise,method=method,cases=len(z),mean_relative_error=float(np.mean(z)),median_relative_error=float(np.median(z))))
with (R/'analysis/controlled_summary.csv').open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(grouped[0]));w.writeheader();w.writerows(grouped)

fig,axs=plt.subplots(2,4,figsize=(16,8),layout='constrained',sharex=True)
for col,d in enumerate(datasets):
    st=d['metadata']['station']
    for row,noise in enumerate([0,1]):
        ax=axs[row,col]
        for method,label,color in zip(methods,labels,colors):
            rr=sorted([a for a in controls if a['station']==st and a['added_noise_factor']==noise and a['method']==method],key=lambda z:z['extra_tstar_s'])
            ax.plot([a['extra_tstar_s'] for a in rr],[a['relative_error'] for a in rr],'o-',label=label,color=color,ms=4)
        ax.set(title=f'{st} / 额外噪声 {noise}σ',xlabel='额外 t* / s',ylabel='相对恢复误差',ylim=(0,None));ax.grid(alpha=.2)
        if row==0 and col==0:ax.legend(fontsize=7)
fig.suptitle('真实波形上的受控额外衰减：目标是恢复到原始已衰减记录')
fig.supxlabel('四条记录来自两个地震。LWLI 为开发记录，其余三条为评估记录；无按误差调参。不同面板纵轴范围可能不同，均从 0 起。',fontsize=9)
fig.savefig(R/'figures/controlled_recovery.png',dpi=150);fig.savefig(R/'figures/controlled_recovery.pdf');plt.close(fig)

d=next(d for d in datasets if d['metadata']['station']=='MNSI');a=next(a for a in d['results'] if a['kind']=='controlled_extra' and a['tstar']==.12 and a['noise_factor']==1);t=np.array(d['t']).ravel();target=np.array(a['target']).ravel()*d['metadata']['scale']
fig,axs=plt.subplots(3,1,figsize=(12,8),layout='constrained',sharex=True)
axs[0].plot(t,target,color='black',label='原真实记录（本试验参考）');axs[0].plot(t,np.array(a['input']).ravel()*d['metadata']['scale'],color='#999999',label='再衰减并加噪后');axs[0].legend(fontsize=9);axs[0].set_title('预先指定示例：第一条评估记录 / 额外 t*=0.12 s / 加噪 1σ')
for method,label,color in zip(methods[1:],labels[1:],colors[1:]):
    out=np.array(a[method]).ravel()*d['metadata']['scale'];axs[1].plot(t,out,label=label,color=color,lw=.8);axs[2].plot(t,out-target,label=label,color=color,lw=.8)
for ax in axs:ax.set(xlim=(-1,7),ylabel='原文件幅度');ax.grid(alpha=.15)
axs[1].plot(t,target,color='black',ls=':',lw=.6,label='参考');axs[1].legend(fontsize=8);axs[1].set_title('三种适配求解输出');axs[2].set_title('相对于原真实记录的误差');axs[2].set_xlabel('相对 P 到时 / s')
fig.savefig(R/'figures/controlled_example.png',dpi=155);fig.savefig(R/'figures/controlled_example.pdf');plt.close(fig)

summary=dict(checks=checks,controlled_summary=grouped,nominal_field_metrics=field,raw_sources_unchanged=True,solver_sources_unchanged=7,
    adaptation='Stationary author tstar + identity offset stage; W=I source waveform sparsity, not original layered Ricker reflectivity',
    no_physical_ground_truth=True,controlled_test_reverses_only_added_attenuation=True,statistics='Descriptive; repeated cases and 3 stations for one event are not independent samples. No significance or generalization claim.')
(R/'analysis/summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=True,indent=2))
