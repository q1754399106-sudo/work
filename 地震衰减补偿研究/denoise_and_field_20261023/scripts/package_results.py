from pathlib import Path
import json,zipfile,hashlib
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];s=json.loads((ROOT/'analysis/summary.json').read_text(encoding='utf-8'))
for key in ['moderate_run.json','strong_run.json']:
 run=Path(json.loads((ROOT/key).read_text(encoding='utf-8'))['run']);ex=json.loads((run/'execution.json').read_text(encoding='utf-8'));assert ex['status']=='COMPLETED' and ex['exit_code']==0
manifest=json.loads((ROOT/'field_data/download_manifest.json').read_text(encoding='utf-8'))
for rec in manifest:assert rec['downloaded'] and hashlib.sha256((ROOT/'field_data'/rec['name']).read_bytes()).hexdigest()==rec['sha256']
for p in (ROOT/'figures').glob('*.png'):
 with Image.open(p) as im:im.verify()
for name in ['wz25','alaska_111','alaska_120','alaska_130']:
 for suffix,keys in [('_processed.npz',['data','denoised','removed']),('_fk.npz',['raw','filtered','removed'])]:
  d=np.load(ROOT/'field_data'/(name+suffix));assert np.linalg.norm(d[keys[0]]-d[keys[1]]-d[keys[2]])<1e-12*max(1,np.linalg.norm(d[keys[0]]))
lines=['# 去干扰与真实地震数据：从这里开始','','已完成两种复杂合成条件的前处理及冻结补偿复算，并下载解析2份真实资料、处理4个真实炮集。','','- index.html：全部合成/实测对比图。','- RESULTS_ZH.md：完整报告，包含无效果尝试、实测未触发和适用边界。','- analysis/compensation_metrics.csv：两种方法等4条分支处理前后的16行评分。','- field_data/L23535.SGY：USGS Alaska NPRA 31-81原始SEG-Y，37.09MB，3030道、30炮。','- field_data/wz.25.H：Alberta公开示例炮集，0.78MB，97道，含发布者插入的中央零道。','- field_data/*_raw.npz、*_processed.npz、*_fk.npz：解析输入、候选输出及f-k对照，保留实际幅度。','- field_data/parsed_metadata.json：采样、几何重建与数据质量备注；download_manifest.json：来源与SHA256。','- baseline/：原复杂试验真值、原结果与固定评分掩码，供对照。','- attempt_v1_broadband/：第一版未触发的完整记录。','- scripts/：检测、分频、f-k、解析、绘图与试验代码。','','结果：低频干扰部分减弱，优化方法整体误差0.56815→0.56330、0.70137→0.68462。较缓线性干扰未可靠去除。','实测候选未触发；f-k基线有减除，但可能同时减去有效波。没有纯净实测真值，不声称恢复误差或SNR改善。','当前实测缺少与原算法匹配的独立Q标定和已确认的层位/速度，也未经动校正，本轮未伪造完整实测Q补偿验证。','','实测重跑：Python环境需numpy/matplotlib/Pillow；read_field.py解析原始文件；run_field.py运行候选；field_fk_comparison.py生成基线对照。','补偿重跑：scripts/run_local.py --entry run_denoise_moderate（或run_denoise_strong） --label replay --timeout 3600 --run-pointer <新的指针路径>。现有观测/标定已随包提供。MATLAB路径按本机配置，已有run指针是本次运行的绝对路径，换目录后需更新。','注明USGS、Madagascar及原始示例来源；数据来源和使用说明详见报告。']
(ROOT/'START_HERE.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
archive=ROOT.parent/'denoise_and_real_seismic_20261023.zip';files=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='DELIVERY_CHECK.json']
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(files):z.write(p,Path(ROOT.name)/p.relative_to(ROOT))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for p in files:assert hashlib.sha256(z.read((Path(ROOT.name)/p.relative_to(ROOT)).as_posix())).digest()==hashlib.sha256(p.read_bytes()).digest()
result=dict(zip=str(archive),bytes=archive.stat().st_size,file_count=len(files),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),verified_all_members=True,verified_original_downloads=True,png_verified=len(list((ROOT/'figures').glob('*.png'))))
(ROOT/'DELIVERY_CHECK.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=True,indent=2))
