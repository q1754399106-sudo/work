from pathlib import Path
import json,hashlib,zipfile
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
summary=json.loads((ROOT/'analysis/summary.json').read_text(encoding='utf-8'));assert summary['checks']['complete']
plan=json.loads((ROOT/'PLAN.json').read_text(encoding='utf-8'));run=Path(json.loads((ROOT/'run.json').read_text(encoding='utf-8'))['run'])
execution=json.loads((run/'execution.json').read_text(encoding='utf-8'));assert execution['status']=='COMPLETED' and execution['exit_code']==0
input=json.loads((ROOT/'core_p1/inputs/complex_input.json').read_text(encoding='utf-8'))
input_errors={}
for cfg in input['datasets']:
 d=json.loads((run/'matlab'/(cfg['name']+'.json')).read_text(encoding='utf-8'));a=np.asarray(cfg['D']);b=np.asarray(d['observed']);err=float(np.linalg.norm(a-b)/np.linalg.norm(a));assert err<1e-14;input_errors[cfg['name']]=err
for p in (ROOT/'figures').glob('*.png'):
 with Image.open(p) as im:im.verify()
records=summary['metrics'];lines=['# 复杂模型对比：从这里开始','','一个2.94km宽、2.2km深、99道的复杂地质模型，两种噪声条件。模型含密集层系、褶皱、双断层、河道、薄层尖灭、低Q透镜体；观测含混合相位、非平稳衰减、多次波、相关噪声、相干干扰和坏道。','这是二维属性驱动的局部一维射线/卷积近似，不是实测或完整二维弹性波模拟。粗几何是构造并假定已知的先验。','','| 条件 | 原分步法误差 | 优化方法误差 | 优化改善的区块 |','|---|---:|---:|---:|']
for name in ['Moderate','Strong']:
 a=next(r for r in records if r['Dataset']==name and r['Method']=='FixedSplit');b=next(r for r in records if r['Dataset']==name and r['Method']=='PrunedFullCalRefit');lines.append(f'| {name} | {a["RelativeError"]:.5g} | {b["RelativeError"]:.5g} | {summary["block_wins_out_of_9"][name]}/9 |')
lines+=['','- index.html：离线完整图集。','- RESULTS_ZH.md：完整报告与局限。','- figures/complex_geology.png：原始地质属性。','- figures/Moderate_comparison.png、Strong_comparison.png：道集及误差。','- analysis/metrics.csv、block_metrics.csv：全部指标及逐区块结果。','- truth/model_truth.npz：地质属性、参考、多次波、噪声、静校正/增益等。','- core_p1/inputs/complex_input.json：求解器实际收到的观测与标定/公开几何，无真实Q。','- runs/：实际执行日志、冻结源码快照、MAT和JSON结果。','','复现：已有完整观测与求解器，可直接运行 scripts/run_local.py --entry run_complex_model --label replay --timeout 5400 --run-pointer <新的run指针文件>。MATLAB路径在run_local.py中按本机配置。','重新生成模型用scripts/build_model.py（固定种子），绘图和评分用scripts/plot_and_report.py；Python需要numpy、matplotlib、Pillow。','原有算法文件保持冻结，没有根据本轮结果重新调参。数值收敛不等于恢复准确，深部/低能量区域和失败区块均保留。']
(ROOT/'START_HERE.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
archive=ROOT.parent/'complex_fieldlike_comparison_20261023.zip'
files=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='DELIVERY_CHECK.json']
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(files):z.write(p,Path(ROOT.name)/p.relative_to(ROOT))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for p in files:assert hashlib.sha256(z.read((Path(ROOT.name)/p.relative_to(ROOT)).as_posix())).digest()==hashlib.sha256(p.read_bytes()).digest()
check=dict(archive=str(archive),bytes=archive.stat().st_size,files=len(files),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),all_zip_member_hashes_verified=True,observed_input_replay_error=input_errors,png_verified=len(list((ROOT/'figures').glob('*.png'))))
(ROOT/'DELIVERY_CHECK.json').write_text(json.dumps(check,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(check,ensure_ascii=True,indent=2))
