from pathlib import Path
import csv,hashlib,json,struct,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
(R/'figures').mkdir(exist_ok=True);(R/'analysis').mkdir(exist_ok=True)
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':11})
download=json.loads((R/'DOWNLOAD_MANIFEST.json').read_text(encoding='utf-8'))
samples=json.loads((R/'SAMPLE_DOWNLOAD_MANIFEST.json').read_text(encoding='utf-8'))
for row in download:
    assert 'error' not in row
    assert hashlib.sha256((R/row['file']).read_bytes()).hexdigest()==row['sha256']
for entry in samples:
    assert 'error' not in entry
    for row in entry['files']:assert hashlib.sha256((R/row['file']).read_bytes()).hexdigest()==row['sha256']
z=np.load(R/'cape/Qp_profile_20-200Hz.npz',allow_pickle=False)
finite=np.isfinite(z['med_Q_chs'])&(z['med_Q_chs']>0)
cape=dict(measured_finite_channels=int(finite.sum()),measured_depth_m=[float(z['dep_chs'][finite].min()),float(z['dep_chs'][finite].max())],measured_Q_range=[float(z['med_Q_chs'][finite].min()),float(z['med_Q_chs'][finite].max())],full_channel_count=1600,
    Q_type='Author-estimated interval P-wave Q, with channel-wise median and MAD; not ground truth',band_hz=[20,200],
    excluded_as_measurements=['Q_interp: includes interpolation and assumed values outside measured range'],
    spectrum_samples=[])
with (R/'analysis/cape_measured_Q.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['depth_m','measured_Q_median','measured_Q_MAD','valid_measured_Q','Vp_m_s','density_g_cm3'])
    for i in range(1600):w.writerow([z['dep_chs'][i],z['med_Q_chs'][i],z['mad_Q_chs'][i],int(finite[i]),z['med_wlogVp_chs'][i],z['med_wlogRho_chs'][i]])
for row in next(a for a in samples if a['dataset']=='cape')['files']:
    zz=np.load(R/row['file'],allow_pickle=False)
    assert not np.iscomplexobj(zz['spec_list_partial'])
    cape['spectrum_samples'].append(dict(file=row['file'],shape=list(zz['spec_list_partial'].shape),frequency_range_hz=[float(zz['f_series'].min()),float(zz['f_series'].max())],phase_available=False))

fig,axs=plt.subplots(1,3,figsize=(12,7),sharey=True,layout='constrained')
q=z['med_Q_chs'];mad=z['mad_Q_chs'];dep=z['dep_chs']
axs[0].errorbar(q[finite],dep[finite]/1000,xerr=mad[finite],fmt='none',ecolor='#aaaaaa',elinewidth=.4,alpha=.35,label='MAD（非置信区间）')
axs[0].scatter(q[finite],dep[finite]/1000,s=5,color='#087f8c',label='作者估计：中位数')
axs[0].scatter(z['Q_interp'],z['dep_interp']/1000,s=15,facecolors='none',edgecolors='#c57529',label='插值 / 假定值',zorder=4)
axs[0].set(xlabel='P 波 Q',ylabel='井深 / km',xlim=(0,max(130,float(np.nanmax(q+mad))*1.05)),title='Q 剖面与离散度');axs[0].legend(fontsize=9)
axs[1].scatter(z['med_wlogVp_chs']/1000,dep/1000,s=5,color='#426b8c');axs[1].set(xlabel='P 波速度 / (km/s)',title='配套测井速度')
axs[2].scatter(z['med_wlogRho_chs'],dep/1000,s=5,color='#426b8c');axs[2].set(xlabel='密度 / (g/cm³)',title='配套测井密度')
axs[0].set_ylim(3.65,0);fig.suptitle('Cape EGS：已下载的真实井中 DAS 衰减参考')
fig.supxlabel('Q 是作者从实测估计的参考值；橙色空心点包含假定与插值，不能作为实测标签。来源：Chang 等，GDR 1823。',fontsize=9)
fig.savefig(R/'figures/cape_Q_reference.png',dpi=160);fig.savefig(R/'figures/cape_Q_reference.pdf');plt.close(fig)

results=list(csv.DictReader((R/'sumatra/results.csv').open(encoding='utf-8-sig')))
waveforms=[];matched=[];arrays=[]
for row in next(a for a in samples if a['dataset']=='sumatra')['files']:
    p=R/row['file'];b=p.read_bytes()
    endian=next(e for e in ['<','>'] if struct.unpack_from(e+'i',b,304)[0] in [6,7])
    delta=struct.unpack_from(endian+'f',b,0)[0];npts=struct.unpack_from(endian+'i',b,316)[0];begin=struct.unpack_from(endian+'f',b,20)[0]
    a=np.frombuffer(b,dtype=endian+'f4',count=npts,offset=632).astype(float)
    assert len(a)==npts and np.isfinite(a).all() and delta>0
    ev,station,component=p.name.split('.')[:3]
    ms=[a for a in results if a['event']==ev and a['station']==station and a['component']=='Z']
    assert any(float(a['alpha'])==0 for a in ms)
    chosen=next(a for a in ms if float(a['alpha'])==0)
    waveforms.append(dict(file=row['file'],event=ev,station=station,channel=component,npts=npts,dt_s=delta,duration_s=(npts-1)*delta,
        author_constant_Q=float(chosen['Q']),author_tstar_s=float(chosen['t*']),author_snr=float(chosen['snr']),available_alpha_values=sorted(set(float(a['alpha']) for a in ms))))
    for refrow in ms:matched.append(dict(waveform_file=row['file'],**refrow))
    tt=begin+np.arange(npts)*delta;arrays.append((tt,a,chosen))
    np.savez_compressed(R/'analysis'/f'{p.stem}_parsed.npz',waveform=np.frombuffer(b,dtype=endian+'f4',count=npts,offset=632).astype(float),t=tt,dt=delta)
with (R/'analysis/sumatra_matched_waveform_Q.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(matched[0]));w.writeheader();w.writerows(matched)

fig,axs=plt.subplots(len(arrays),1,figsize=(13,9),layout='constrained')
for ax,(tt,a,c),meta in zip(axs,arrays,waveforms):
    ax.plot(tt,a,lw=.5,color='#315f7b');ax.set_title(f"{meta['event']} / {meta['station']} / Z：作者 Q={meta['author_constant_Q']:.2f}，α=0，t*={meta['author_tstar_s']:.4f} s",fontsize=10);ax.set_ylabel('原文件幅度')
axs[-1].set_xlabel('相对于 SAC 参考时刻的时间 / s')
fig.suptitle('Sumatra：已下载并与作者 Q 表逐事件 / 台站匹配的波形')
fig.supxlabel('原数据已做 0.2 Hz 高通；本图未再增益或归一化。自然地震记录，非反射炮集。Q 是路径参考估计，非介质真值。',fontsize=9)
fig.savefig(R/'figures/sumatra_matched_waveforms.png',dpi=150);fig.savefig(R/'figures/sumatra_matched_waveforms.pdf');plt.close(fig)

summary=dict(cape=cape,sumatra=dict(result_rows=len(results),sample_waveforms=waveforms,
    qualifications='Natural earthquakes, not active reflection shots. Alpha=0 subset used only as explicit constant-Q reference; frequency-dependent variants retained.'),
    poseidon=dict(status='Candidate only: published multi-well VSP Q result verified; matching raw VSP download and numerical Q table not verified',paper='https://library.hagi.or.id/wp-content/uploads/2024/01/197.pdf',portal='https://terranubis.com/datainfo/NW-Shelf-Australia-Poseidon-3D'),
    downloadable_examples_verified=True,no_experiment_run=True)
(R/'analysis/VERIFIED_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 附带 Q 参考值的真实地震资料检索与下载核验','','核验日期：2026-10-23。任务：为实测衰减研究寻找外部 Q 参考，避免对无标定炮集填入假定 Q。','','## 已核实并下载的资料','','| 数据 | 附带的 Q | 实际下载核验 | 对当前课题的用途 |','|---|---|---|---|',
'| Cape EGS / Delano 1OB 井中 DAS | 作者估计的深度 Qp、中位数与 MAD；20-200 Hz | Q 剖面、井轨迹、速度密度、事件表、读取代码及 2 个事件的振幅谱样本 | 衰减估计和频谱振幅补偿的外部参考；不支持完整时域相位比较 |',
'| Sumatra 自然地震 | 按 α / 事件 / 台站 / 分量提供 Q、t*、SNR 等 | 完整结果表及 4 条与 α=0 结果严格匹配的 SAC 波形 | 可做时域和频域算法适配试验；与反射炮集的物理场景不同 |','',
'## 1. Cape EGS：最直接的 Q 深度参考','',
f"实际文件有 1600 个通道，{cape['measured_finite_channels']} 个通道包含有限正的 med_Q_chs。有效 Q 的对应深度约 {cape['measured_depth_m'][0]:.1f}-{cape['measured_depth_m'][1]:.1f} m；Q 中位数取值约 {cape['measured_Q_range'][0]:.2f}-{cape['measured_Q_range'][1]:.2f}。这些是文件核验结果，不是对所有深度做出的测量。",'',
'med_Q_chs 为作者的估计中位数；mad_Q_chs 是测量离散度，不是 95% 置信区间。Q_interp 明确包含测量范围外的插值与假定值，本次单独标注，不混作实测标签。',
'已抽取两个事件的振幅谱与噪声谱，数据为实数幅度而非复频谱，没有波形相位。不能据此合成“真实无衰减波形”或宣称完成实测相位补偿验证。公开完整频谱 ZIP 为 16,378,110,956 字节，本次通过 HTTP Range 只抽取少量成员，未下载全包。','',
'数据：https://gdr.openei.org/submissions/1823',
'对应论文：https://doi.org/10.1029/2026JB034276',
'许可：CC BY 4.0；使用时注明 Chang 等及数据来源。','',
'## 2. Sumatra：完整波形与逐路径 Q','',
'results.csv 有 45,133 行，包含不同 α 假设的重复估计，不能把行数当独立波形数量。本次按 α=0、Z 分量、作者 SNR≥5 选择在波形 ZIP 中有严格事件与台站匹配的前四条记录。其他 α 结果保留，不混作同一常 Q。',
'每条波形均保留 SAC 原件，另提供数值解析 NPZ；已检查字节序、采样间隔、样点数、有限值和归档成员 CRC。Q 参考与波形的匹配表为 analysis/sumatra_matched_waveform_Q.csv。','','| 事件 | 台站 | 作者 Q（α=0） | t* / s | dt / s |','|---|---|---:|---:|---:|']
for a in waveforms:lines.append(f"| {a['event']} | {a['station']} | {a['author_constant_Q']:.3f} | {a['author_tstar_s']:.6f} | {a['dt_s']:.3f} |")
lines+=['',
'这些是已做 0.2 Hz 高通的自然地震位移记录，不是陆上反射炮集。样本 dt 约 0.05 s，对应采样率 20 Hz、奈奎斯特频率 10 Hz，因此不能直接套用当前面向 25-60 Hz 反射子波的参数，需要重新适配频带、波形窗和正演算子。路径 Q 不能直接当作层间 Q，更不能移植到 Alaska 测线。无真实未衰减波形参考。波形文件是从约 601 MB 的公开 ZIP 中按范围请求抽取的样本。',
'数据与作者原始说明：https://zenodo.org/records/6932666',
'许可：元数据标为 CC BY 4.0。','',
'## 3. 更接近毕设的主动源候选：Poseidon-1 / Pharos-1','',
'作者论文使用澳大利亚 Browse Basin 两口井的零偏移 VSP，给出多谱比 Q 深度剖面。数据门户说明另外提供 VSP / 测井资料；主下载本身只含全叠加及井轨迹，不能误认为已经拿到原始 VSP 和数值 Q 表。当前未核实这两项的配套下载，因此仅列候选。',
'主门户的 8-bit、AGC 叠加体不适合作为保幅实验原始输入；若继续，应核对原始 32-bit VSP、实际时深关系与原始处理报告。',
'作者论文：https://library.hagi.or.id/wp-content/uploads/2024/01/197.pdf',
'数据门户：https://terranubis.com/datainfo/NW-Shelf-Australia-Poseidon-3D','',
'## 建议使用次序','',
'先使用 Cape 的实测 Qp 及 MAD 做衰减参考和不确定性分析；若需要完整波形试验，使用 Sumatra 已匹配样本，但单独报告场景差异。反射炮集对比仍应等待核实到配套 VSP 标定的主动源资料，不将微地震或自然地震包装为炮集。',
'所有真实 Q 都是带方法假设的估计；本轮没有找到已知精确 Q 真值及未衰减参考的实测反射炮集。','',
'## 交付与复现','',
'DOWNLOAD_MANIFEST.json 保存小文件下载来源与 SHA256；SAMPLE_DOWNLOAD_MANIFEST.json 保存远程 ZIP 范围请求、CRC、成员哈希；完整中央目录单独保存。verify_sources.py 下载结果文件；fetch_samples.py 提取样本；inspect_and_package.py 核验字段与匹配并作图。未执行下载的作者代码。']
(R/'CATALOG_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
html='<!doctype html><meta charset="utf-8"><title>真实 Q 参考资料</title><style>body{max-width:1200px;margin:30px auto;font:18px sans-serif;line-height:1.7;color:#16354b}img{width:100%}a{color:#087f8c}</style><h1>附带 Q 参考值的真实资料</h1><p>已核实 Cape EGS 的 Q 深度剖面和 Sumatra 的波形 / Q 匹配样本。Q 均为作者估计参考，非精确真值；两套均非主动源反射炮集。</p><p><a href="CATALOG_ZH.md">完整检索说明、来源与限制</a></p><h2>Cape：实测估计与假定区间分开标注</h2><img src="figures/cape_Q_reference.png"><h2>Sumatra：与作者 Q 表匹配的完整波形</h2><img src="figures/sumatra_matched_waveforms.png">'
(R/'index.html').write_text(html,encoding='utf-8')

approved_samples={str((R/f['file']).resolve()) for a in samples for f in a['files']}
dest=R.parent/'q_reference_data_20261023.zip'
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as zipout:
    for p in sorted(R.rglob('*')):
        if not p.is_file() or p.name.startswith('.') or '__pycache__' in p.parts or p.suffix=='.log' or p.name=='DELIVERY_CHECK.json':continue
        if 'samples' in p.parts and str(p.resolve()) not in approved_samples:continue
        zipout.write(p,'q_reference_data_20261023/'+str(p.relative_to(R)).replace('\\','/'))
with zipfile.ZipFile(dest) as zz:assert zz.testzip() is None
check=dict(zip=str(dest),bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),verified_downloads=len(download),verified_spectrum_samples=2,verified_matched_SAC_waveforms=4)
(R/'DELIVERY_CHECK.json').write_text(json.dumps(check,indent=2),encoding='utf-8')
print(json.dumps(dict(cape=cape,sumatra_waveforms=waveforms,delivery=check),ensure_ascii=True,indent=2))
