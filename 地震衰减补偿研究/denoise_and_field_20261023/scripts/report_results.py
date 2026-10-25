from pathlib import Path
import json,csv,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'complex_fieldlike_20261023';FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
has_old=(OLD/'run.json').exists()
baseline_run=Path(json.loads((OLD/'run.json').read_text(encoding='utf-8'))['run']) if has_old else ROOT/'baseline'
TRUTH=np.load(OLD/'truth/model_truth.npz' if has_old else ROOT/'baseline/model_truth.npz');t=TRUTH['t'];x=TRUTH['x'];live=TRUTH['live']
checks={};rows=[];preprocessing=[];pairs={};hashes=json.loads((ROOT/'FROZEN_SOLVER_HASHES.json').read_text(encoding='utf-8'))
basecopy=ROOT/'baseline';basecopy.mkdir(exist_ok=True)
if has_old:shutil.copy2(OLD/'truth/model_truth.npz',basecopy/'model_truth.npz');shutil.copy2(OLD/'PLAN.json',basecopy/'complex_model_plan.json')
def metric(a,target,lowmask):
 e=a-target;deep=(t>=1.05)&(t<=1.6)
 return dict(RelativeError=float(np.linalg.norm(e)/np.linalg.norm(target)),LiveError=float(np.linalg.norm(e[:,live])/np.linalg.norm(target[:,live])),DeepError=float(np.linalg.norm(e[deep,:][:,live])/np.linalg.norm(target[deep,:][:,live])),LowEnergyError=float(np.linalg.norm(e[lowmask])/np.linalg.norm(target[lowmask])))
for name,ptr in [('Moderate','moderate_run.json'),('Strong','strong_run.json')]:
 run=Path(json.loads((ROOT/ptr).read_text(encoding='utf-8'))['run']);c=json.loads((run/'matlab/checks.json').read_text(encoding='utf-8'));assert c['complete'] and c['blocks']==9;checks[name]=c
 assert all(hashlib.sha256((run/'source/core_p1'/p).read_bytes()).hexdigest()==h for p,h in hashes.items())
 oldpath=(baseline_run/'matlab' if has_old else baseline_run)/(name+'.json');before=json.loads(oldpath.read_text(encoding='utf-8'));after=json.loads((run/'matlab'/(name+'.json')).read_text(encoding='utf-8'))
 if has_old:shutil.copy2(oldpath,basecopy/(name+'.json'))
 maskpath=(OLD/'analysis' if has_old else basecopy)/(name+'_evaluation_masks.npz');mask=np.load(maskpath)['low_energy']
 if has_old:shutil.copy2(maskpath,basecopy/maskpath.name)
 target=TRUTH[name+'_target'];pre=np.load(ROOT/'analysis'/(name+'_preprocess.npz'))
 assert np.linalg.norm(np.asarray(after['observed'])-pre['denoised'])<1e-12
 for method,key in [('FixedSplit','fixed'),('FixedSplitSharedCal','fixed_shared'),('FullCalRefit','unpruned'),('PrunedFullCalRefit','optimized')]:
  for stage,d in [('Before',before),('After',after)]:rows.append(dict(Dataset=name,Method=method,Stage=stage,**metric(np.asarray(d[key]),target,mask)))
 coh=TRUTH[name+'_coherent'][:,live];removed=pre['removed'][:,live]
 preprocessing.append(dict(Dataset=name,RemovedInputNormFraction=float(np.linalg.norm(pre['removed'])/np.linalg.norm(pre['original'])),CoherentModelResidualRatio=float(np.linalg.norm(coh-removed)/np.linalg.norm(coh)),RemovedCoherentCorrelation=float(np.sum(coh*removed)/(np.linalg.norm(coh)*np.linalg.norm(removed))) if np.linalg.norm(removed)>0 else 0,CleanInputControlChange=json.loads((ROOT/'analysis/clean_input_controls.json').read_text(encoding='utf-8'))[name]['relative_change']))
 arrays=[target,np.asarray(before['fixed']),np.asarray(after['fixed']),np.asarray(before['observed']),np.asarray(before['optimized']),np.asarray(after['optimized'])];limit=max(float(np.abs(a).max()) for a in arrays)
 fig,axs=plt.subplots(2,3,figsize=(16,10),sharex=True,sharey=True,layout='constrained')
 for ax,a,title in zip(axs.flat,arrays,['纯净参考','原分步法：直接补偿','原分步法：先去干扰再补偿','原始观测','优化方法：直接补偿','优化方法：先去干扰再补偿']):
  im=ax.imshow(a,extent=[x[0]/1000,x[-1]/1000,t[-1],0],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_title(title,fontsize=11);ax.set_xlabel('横向位置 / km')
 axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='六面板共用振幅色标',shrink=.85);fig.suptitle(name+'：独立去干扰前处理对补偿结果的影响',fontsize=16)
 fig.supxlabel('现有补偿代码与标定保持冻结；去噪后重新估计观测噪声水平；不使用真实干扰轨迹，不逐道归一化。',fontsize=9)
 fig.savefig(FIG/(name+'_compensation_comparison.png'),dpi=150);fig.savefig(FIG/(name+'_compensation_comparison.pdf'));plt.close(fig)
 arrays=[pre['original'],pre['denoised'],pre['removed']];limit=float(np.quantile(np.abs(np.concatenate([a.ravel() for a in arrays])),.995));fig,axs=plt.subplots(1,3,figsize=(15,7),sharex=True,sharey=True,layout='constrained')
 for ax,a,title in zip(axs,arrays,['去干扰前','去干扰后','减掉的分量']):
  im=ax.imshow(a,extent=[0,x[-1]/1000,t[-1],0],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_title(title);ax.set_xlabel('横向位置 / km')
 axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='共同振幅色标',shrink=.85);fig.suptitle(name+'：分频长孔径相干波减除',fontsize=16)
 fig.supxlabel('共同色标为三面板合并绝对振幅99.5%分位，仅显示截幅；没有逐道增益，不读取真值。',fontsize=9)
 fig.savefig(FIG/(name+'_preprocessing.png'),dpi=155);fig.savefig(FIG/(name+'_preprocessing.pdf'));plt.close(fig)

def find(name,method,stage):return next(r for r in rows if r['Dataset']==name and r['Method']==method and r['Stage']==stage)
controls=json.loads((ROOT/'analysis/clean_input_controls.json').read_text(encoding='utf-8'));gates={}
for name in ['Moderate','Strong']:
 flags=[]
 if controls[name]['relative_change']>.02:flags.append('Clean signal changed >2%')
 for method in ['FixedSplit','PrunedFullCalRefit']:
  a=find(name,method,'Before');b=find(name,method,'After')
  if b['RelativeError']>1.02*a['RelativeError']:flags.append(method+' waveform regression >2%')
  if b['LowEnergyError']>1.05*a['LowEnergyError']:flags.append(method+' low-energy regression >5%')
 gates[name]=dict(passed=not flags,reasons=flags)
field=json.loads((ROOT/'field_data/experiment_diagnostics.json').read_text(encoding='utf-8'));fk=json.loads((ROOT/'field_data/fk_diagnostics.json').read_text(encoding='utf-8'))
summary=dict(checks=checks,metrics=rows,preprocessing=preprocessing,gates=gates,field=dict(adaptive_removed_norm_fraction={k:v['removed_norm_fraction'] for k,v in field.items()},fk_removed_norm_fraction={k:v['removed_norm_fraction'] for k,v in fk.items()}),promotion='No default promotion: field adaptive detector did not trigger, no field Q calibration or clean truth. Research trial retained.',original_matlab_files_unchanged=len(hashes))
(ROOT/'analysis/summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
for filename,records in [('compensation_metrics.csv',rows),('preprocess_metrics.csv',preprocessing)]:
 with (ROOT/'analysis'/filename).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
lines=['# 相干干扰压制与真实地震资料试验','','## Material Passport','','类型：合成资料受控对照与真实炮集去噪探索。当前结果为研究尝试，不替换原补偿入口。','','## 方法与试验过程','','第一版宽频长孔径相干检测未达到0.65相似性门槛，未改变输入，完整记录在attempt_v1_broadband。第二版在按观测f90划分的频带内检测长直线相干波，预测时按道RMS预条件化，减除前恢复原幅度单位。','全孔径扫描视速度>=400m/s，1201个斜率，至少70%有效孔径，峰值相似性>=0.65，每频带最多3个波包；减除比例0.60。斜率、到时由观测估计，没有读生成器中的1850m/s、10000m/s或轨迹。','补偿与标定冻结，预处理后重新估计噪声水平。两个原复杂样本属于开发资料，不能据此宣称独立泛化。纯净输入控制只用于评分，不用于选取参数。','','## 合成数据补偿前后','','| 条件 | 方法 | 直接补偿误差 | 去干扰后补偿误差 | 低能量区：前→后 |','|---|---|---:|---:|---|']
for name in ['Moderate','Strong']:
 for method in ['FixedSplit','FixedSplitSharedCal','FullCalRefit','PrunedFullCalRefit']:
  a=find(name,method,'Before');b=find(name,method,'After');lines.append(f'| {name} | {method} | {a["RelativeError"]:.6g} | {b["RelativeError"]:.6g} | {a["LowEnergyError"]:.5g} → {b["LowEnergyError"]:.5g} |')
lines+=['','## 合成干扰建模与保真检查','','| 条件 | 减掉的输入范数比例 | 相干干扰建模残差比 | 与真实干扰相关系数 | 纯净输入变化 |','|---|---:|---:|---:|---:|']
for r in preprocessing:lines.append('| '+r['Dataset']+' | '+' | '.join(f'{r[k]:.6g}' for k in ['RemovedInputNormFraction','CoherentModelResidualRatio','RemovedCoherentCorrelation','CleanInputControlChange'])+' |')
lines+=['','相干建模指标只用于合成评分，真实干扰分量没有传给前处理。范数比例不是能量比例。相干干扰建模残差比以原干扰范数为分母，1表示未减除，不能等同于整体补偿误差。','候选检测出约1846m/s的低频波包；较缓、与反射频带重叠的线性干扰未被可靠识别。',f'保真门槛结果：{json.dumps(gates,ensure_ascii=False)}。数值检查与恢复效果分别报告。','','## 真实资料：已下载、解析并实际试验','','| 数据 | 样本与采样 | 本轮实际使用 |','|---|---|---|','| USGS Alaska NPRA 31-81 | 1981年陆地炸药源；SEG-Y IBM32；2ms、6s | L23535.SGY含3030道/30炮，选择111、120、130炮各前96道 |','| Yilmaz wz.25 | Alberta陆地炸药源；XDR浮点；2ms、4s | 97道，其中中央道是原示例插入的零道，保留并标记 |','','USGS原始文件37,090,800字节；wz.25.H为777,284字节。下载URL、时间、HTTP头、SHA256和原始文件头均保存。阿拉斯加原offset头为0，偏移距按发布者处理脚本的-5225ft、110ft间隔重建；不是原始导航坐标。wz25按文件头50m道距，头文件也保留历史几何不确定性备注。','','## 真实资料去噪结果','','| 炮集 | 长孔径分频减掉的范数比例 | 软f–k基线减掉的范数比例 |','|---|---:|---:|']
for name in field:lines.append(f'| {name} | {field[name]["removed_norm_fraction"]:.6g} | {fk[name]["removed_norm_fraction"]:.6g} |')
lines+=['','长孔径候选在四个实测炮集上均未触发，不能称为实测去噪成功。额外软f–k基线按1200—2500m/s过渡、观测f90的0.4—0.8倍频率过渡压制低频慢速谱锥，最大减除0.60；属于另一个探索基线，没有冒充主候选。','f–k减掉的能量可能包含真实直达波、折射波或反射分量，不能用去除量代替SNR或保幅改善。图中显示原始、候选、f–k和被减分量。t²只用于全部面板的同一显示增益，不参与处理；共用99.5%分位色标，显示截幅比例记录在诊断文件。','这些真实资料没有配套的独立Q标定道集和已确认的层厚/速度、真实无衰减参考。实测试验因此限于去干扰；没有填入假定Q冒充估计Q，没有宣称完成可信的实测Q保幅补偿对比。当前真实炮集还未经动校正，不能直接把横向地质/到时变化当AVO归一化。','','## 数据来源与使用说明','','- USGS NAMSS原始调查目录：https://walrus.wr.usgs.gov/namss/survey/n-pr-81-ak/','- Alaska数据处理与几何：https://ahay.org/RSF/book/data/alaska/line31-81-m8r.html','- 实际SEG-Y下载：https://ahay.org/data/alaska/L23535.SGY','- Yilmaz实测去面波论文复现代码：https://ahay.org/RSF/book/tccs/orthogroll/field.html','- 实际示例数据：https://ahay.org/data/wz/wz.25.H','- Alaska使用说明（注明USGS来源）：https://wiki.seg.org/wiki/ALASKA_2D_LAND_LINE_31-81','- wz25下载页没有单独的数据许可证文件，本报告保留来源，不把软件许可证当作数据许可证。','','解析脚本scripts/read_field.py明确IBM32和XDR格式、道选择和几何重建；处理代码coherent_filter.py、multiband_filter.py及fk_baseline.py均保留。原始文件不覆盖。']
(ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>去干扰与真实资料试验</title><style>body{font:16px system-ui;max-width:1500px;margin:24px auto;padding:0 20px}img{width:100%}p{line-height:1.7}</style><h1>去干扰与真实地震资料试验</h1><p>宽频无效果尝试、分频候选及实测软f-k基线均保留。<a href="RESULTS_ZH.md">完整报告</a> · <a href="analysis/compensation_metrics.csv">合成指标</a></p>'''
for name in ['Moderate','Strong']:page+=f'<h2>合成：{name}</h2><img src="figures/{name}_preprocessing.png" alt="前处理与减掉分量"><img src="figures/{name}_compensation_comparison.png" alt="补偿前后对照">'
for name in ['wz25','alaska_111','alaska_120','alaska_130']:page+=f'<h2>真实炮集：{name}</h2><p>候选未触发；另列软f-k探索基线。没有纯净真值，不声称保幅恢复精度。</p><img src="figures/{name}_field_comparison.png" alt="真实资料前处理尝试">'
(ROOT/'index.html').write_text(page+'</html>',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=True,indent=2),flush=True)
