"""Attach declared scope and all regional contrasts; primary scores unchanged."""
from common import *
p=ROOT/'RESULTS_ZH.md';s=p.read_text(encoding='utf-8').split('\n<!--final-evidence-notes-->')[0]
if '- 开发预算限制：' not in s:
    s=s.replace('- 方法：','- 开发预算限制：SC仅用F1/F2选参，当前相干规则历史上使用过更多模型；不是完全对等的历史开发预算。见[SCOPE_ZH.md](SCOPE_ZH.md)。\n- 方法：',1)
a=json.loads((ROOT/'analysis/regional_posthoc.json').read_text(encoding='utf-8'))['rows']
extras=['','<!--final-evidence-notes-->','## 区域诊断（事后，不改变主门槛）','',
'为检查整体误差是否掩盖反射结构与弱事件的取舍，按上游score_gather/score_weak原定义额外复算；完整128行见analysis/regional_posthoc.csv。此处所有有对应真值的断层、尖灭和弱反射比较均列出，不将有利区域改成主结论。','',
'| 场景 | 指标 | SC | Standard归档 |','|---|---|---:|---:|']
cases=list(dict.fromkeys(x['case'] for x in a))
for c in cases:
    for metric in ['fault','pinch','weak_projection_error']:
        b=[x['value'] for x in a if x['case']==c and x['method']=='SC' and x['metric']==metric]
        z=[x['value'] for x in a if x['case']==c and x['method']=='StandardArchived' and x['metric']==metric]
        if b and z:extras.append(f'| {c} | {metric} | {b[0]:.6f} | {z[0]:.6f} |')
extras+=['','V4上Standard的断层/尖灭区域误差较低，但弱反射投影误差反而高于SC（两种噪声均如此）。因此，本轮同样不支持“弱反射保留优于这一文献基线”的主张。区域结果只是取舍诊断，不构成新的确认门槛。','',
'## 完整复现与运行修复','',
'12份实际使用的MATLAB重拟合状态已归档在inputs/prepared，并给出SHA256；可先运行scripts/restore_prepared.py再运行主报告入口。此举不冒充重新执行稀疏支持发现阶段。',
'首次评估因空支持分组边界实现问题主动中止，原源码及日志已保留；增加与MATLAB一致的空支持分支后，选定参数和门槛未变，完整评估重跑成功，评估期间实际依赖哈希无变化。',
'如需重建本报告附注，在make_report.py后运行regional_diagnostic.py及finalize_report.py。开发预算和新颖性边界详见SCOPE_ZH.md，不把本轮描述性比较推广为全球首创或实测优越性。','']
p.write_text(s+'\n'.join(extras),encoding='utf-8')
files={}
for f in ROOT.rglob('*'):
    if not f.is_file() or any(x in ('cache','__pycache__') for x in f.relative_to(ROOT).parts):continue
    if f.name=='ARTIFACT_MANIFEST.json':continue
    files[str(f.relative_to(ROOT)).replace('\\','/')]=dict(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
writejson(ROOT/'ARTIFACT_MANIFEST.json',dict(files=files,scope='Task artifacts only; ignored cache/dependencies excluded.'))
print('Finalized report; manifest files:',len(files))
