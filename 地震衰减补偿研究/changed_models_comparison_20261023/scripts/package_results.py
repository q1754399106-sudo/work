from pathlib import Path
import json,hashlib,zipfile
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
s=json.loads((ROOT/'summary.json').read_text(encoding='utf-8'))
manifest=json.loads((ROOT/'figures/figure_manifest.json').read_text(encoding='utf-8'))
assert len(manifest)==100 and len(list((ROOT/'figures/wiggles').glob('*.png')))==8
for rec in manifest:
 assert hashlib.sha256((ROOT/rec['png']).read_bytes()).hexdigest()==rec['image_sha256']
 assert hashlib.sha256(Path(rec['source']).read_bytes()).hexdigest()==rec['source_sha256']
for p in (ROOT/'figures').rglob('*.png'):
 with Image.open(p) as im:im.verify()
page=(ROOT/'index.html').read_text(encoding='utf-8')
page=page.replace('<p><a href="RESULTS_ZH.md">','<img src="figures/paired_errors.png" alt="所有案例的成对波形误差"><img src="figures/weak_amplitudes.png" alt="固定案例的弱反射逐道振幅"><p><a href="RESULTS_ZH.md">')
(ROOT/'index.html').write_text(page,encoding='utf-8')
normal=s['suites']['normal']['means'];colored=s['suites']['colored']['means'];wins=s['suites']['normal']['paired']['FixedSplit']['RelativeError']
text=f'''# 本轮更换模型对比

已完成4个新地质配置、80组普通噪声及20组有色噪声，原有求解代码保持冻结。

普通噪声整体误差：原分步法{normal['FixedSplit']['RelativeError']:.4f}，优化方法{normal['PrunedFullCalRefit']['RelativeError']:.4f}；优化方法在{wins['wins']}/{wins['n']}个案例更低。
普通噪声可辨识弱反射平均误差：{normal['FixedSplit']['WeakEventError']:.4f} → {normal['PrunedFullCalRefit']['WeakEventError']:.4f}；{len(s['normal_weak_regression_models'])}个模型弱反射退化。
有色噪声弱反射平均误差：{colored['FixedSplit']['WeakEventError']:.4f} → {colored['PrunedFullCalRefit']['WeakEventError']:.4f}。

- RESULTS_ZH.md：完整报告，包括失败案例和适用范围。
- index.html：离线道集图库，全部100组均可查看。
- figures/comparison.png：逐模型误差对比。
- figures/weak_amplitudes.png：4个固定案例的弱反射振幅曲线。
- figures/worst_relative_case.png：按误差比自动选出的最差案例。
- metrics.csv：原始500行评分；metrics_identifiable.csv：统一排除不可辨识弱振幅后的评分；summary.json：两种评分及汇总。
- PLAN.json / FROZEN_SOLVER_HASHES.json：运行前登记计划与原有求解代码哈希。
- runs/：MATLAB日志、源码快照、独立校准和原始MAT/JSON。

主对照是原分步法FixedSplit；FixedSplitSharedCal另使用与优化方法相同的完整校准结果。
四层、11道、Ricker子波正演体系保持不变；本轮检验新地质配置，不代表实测资料验证。
薄层WeakSection中央道有两个反射完全重合。弱反射评分在运行后复核时修正，统一排除这一不可辨识分量；原始评分保留，完整波形评分不变。详见报告的评分修正说明。
'''
(ROOT/'START_HERE.md').write_text(text,encoding='utf-8')
archive=ROOT.parent/'changed_models_comparison_20261023.zip'
files=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='DELIVERY_CHECK.json']
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(files):z.write(p,Path(ROOT.name)/p.relative_to(ROOT))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for p in files:
  assert hashlib.sha256(z.read((Path(ROOT.name)/p.relative_to(ROOT)).as_posix())).digest()==hashlib.sha256(p.read_bytes()).digest()
check=dict(png_files_verified=len(list((ROOT/'figures').rglob('*.png'))),gathers=100,wiggles=8,archive=str(archive),files=len(files),bytes=archive.stat().st_size,sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),zip_crc_and_all_members_verified=True)
(ROOT/'DELIVERY_CHECK.json').write_text(json.dumps(check,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(check,ensure_ascii=True,indent=2))
