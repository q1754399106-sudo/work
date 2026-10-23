from pathlib import Path
import csv,hashlib,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;S=R.parent/'q_reference_data_20261023/cape'
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
z=np.load(S/'Qp_profile_20-200Hz.npz',allow_pickle=False);dep=z['dep_chs'];Q=z['med_Q_chs'];mad=z['mad_Q_chs'];vp=z['med_wlogVp_chs']
valid=np.isfinite(Q)&np.isfinite(mad)&np.isfinite(vp)&(Q>0)&(vp>0)
rows=[];full=[]
for p in sorted((S/'samples').glob('*.npz')):
    d=np.load(p,allow_pickle=False);chs=d['ch_list_partial'];ds=dep[chs];ok=valid[chs]
    depth=ds[ok];lo,hi=np.quantile(depth,[.25,.75]);assert hi-lo>200
    up=abs(ds-lo)<=40;down=abs(ds-hi)<=40
    up&=ok;down&=ok;assert up.sum()>=5 and down.sum()>=5
    lo=float(np.median(ds[up]));hi=float(np.median(ds[down]))
    grid=np.linspace(lo,hi,300);inside=valid&(dep>=lo-100)&(dep<=hi+100);ids=np.argsort(dep[inside]);zz=dep[inside][ids]
    vg=np.interp(grid,zz,vp[inside][ids]);qg=np.interp(grid,zz,Q[inside][ids]);mg=np.interp(grid,zz,mad[inside][ids])
    curves=[qg,qg+mg,np.maximum(qg-mg,.2*qg)]
    stars=[float(np.trapz(1/(vg*q),grid)) for q in curves]
    f=d['f_series'];upper=np.median(d['spec_list_partial'][up],axis=0);lower=np.median(d['spec_list_partial'][down],axis=0)
    nu=np.median(d['spec_n_list_partial'][up],axis=0);nd=np.median(d['spec_n_list_partial'][down],axis=0)
    band=(f>=20)&(f<=200);good=band&(upper>=3*nu)&(lower>=3*nd)&(upper>0)&(lower>0)
    limit=10.;rawg=np.exp(np.minimum(np.pi*f[:,None]*np.array(stars)[None,:],50));gains=np.minimum(rawg,limit)
    corrected=upper[:,None]*gains
    b0=b1=None
    if good.sum()>=10:
        b0=float(np.polyfit(f[good],np.log(upper[good]/lower[good]),1)[0]);b1=float(np.polyfit(f[good],np.log(corrected[good,0]/lower[good]),1)[0])
    row=dict(event=p.stem,shallow_depth_m=lo,deep_depth_m=hi,shallow_channels=int(up.sum()),deep_channels=int(down.sum()),
             vertical_interval_tstar_s=stars[0],tstar_Q_plus_MAD=stars[1],tstar_Q_minus_MAD=stars[2],
             usable_bins=int(good.sum()),raw_ratio_slope_per_hz=b0,corrected_ratio_slope_per_hz=b1,
             median_nominal_gain_in_band=float(np.median(gains[band,0])),max_gain_in_band=float(gains[band,0].max()),
             fraction_gain_capped_in_band=float(np.mean(rawg[band,0]>limit)),
             interpretation='Conditional vertical upgoing path reweighting only. DAS directional/gauge response not independently corrected. No phase/waveform reconstruction.')
    rows.append(row)
    full.append(dict(**row,f=f.tolist(),raw_upper=upper.tolist(),raw_lower=lower.tolist(),noise_upper=nu.tolist(),good=good.tolist(),corrected=corrected.tolist(),gains=gains.tolist()))
    np.savez_compressed(R/'analysis'/f'{p.stem}_spectral_trial.npz',frequency=f,shallow=upper,deep=lower,noise_shallow=nu,noise_deep=nd,corrected=corrected,gains=gains,usable=good)
    fig,axs=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    axs[0].semilogy(f,upper,label='较浅窗原始谱');axs[0].semilogy(f,lower,label='较深窗原始谱');axs[0].semilogy(f,nu,ls=':',color='gray',label='较浅窗噪声谱');axs[0].set(xlim=(20,200),xlabel='频率 / Hz',ylabel='原始振幅谱',title=f'深度窗 {lo:.0f} / {hi:.0f} m');axs[0].legend(fontsize=8)
    for j,label in enumerate(['作者 Q 中位数','Q + MAD','Q - MAD（正值截限）']):axs[1].plot(f,gains[:,j],label=label)
    axs[1].axhline(10,color='gray',ls=':');axs[1].set(xlim=(20,200),xlabel='频率 / Hz',ylabel='幅度乘数',title='条件路径增益；上限 10 倍');axs[1].legend(fontsize=8)
    axs[2].plot(f[band],np.log(np.maximum(upper[band],1e-30)/np.maximum(lower[band],1e-30)),label='原浅 / 深对数谱比',alpha=.6)
    axs[2].plot(f[band],np.log(np.maximum(corrected[band,0],1e-30)/np.maximum(lower[band],1e-30)),label='按作者 Q 重加权后',alpha=.7)
    axs[2].scatter(f[good],np.log(corrected[good,0]/lower[good]),s=8,color='black',label='两窗 SNR≥3 的频点')
    axs[2].set(xlabel='频率 / Hz',ylabel='ln(浅 / 深振幅谱)',title='谱比检查；不等同 Q 真值检验');axs[2].legend(fontsize=8)
    fig.suptitle('Cape '+p.stem.replace('spec_D-1OB_P_arr_',''))
    fig.supxlabel('只使用有实测估计的 med_Q_chs，不使用 Q_interp；区间内插值；假设竖直向上传播；未独立校正 DAS 响应。MAD 非置信区间。',fontsize=9)
    fig.savefig(R/'figures'/(p.stem+'_cape_trial.png'),dpi=155);fig.savefig(R/'figures'/(p.stem+'_cape_trial.pdf'));plt.close(fig)
    shutil.copy2(p,R/'raw'/p.name)
for name in ['Qp_profile_20-200Hz.npz','well_config.npz','README']:shutil.copy2(S/name,R/'raw'/('cape_'+name))
with (R/'analysis/cape_metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'analysis/cape_results.json').write_text(json.dumps(full,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(rows,ensure_ascii=True,indent=2))
