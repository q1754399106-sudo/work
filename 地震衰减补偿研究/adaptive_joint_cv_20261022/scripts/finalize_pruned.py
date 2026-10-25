from pathlib import Path
import json,csv,hashlib
from pruning_gates import summarize,gates
ROOT=Path(__file__).resolve().parents[1]
def run(pointer):return Path(json.loads((ROOT/pointer).read_text(encoding='utf-8'))['run'])
dev=json.loads((ROOT/'analysis/pruned_development_summary.json').read_text(encoding='utf-8'));val=run('pruned_validation_run.json');checks=json.loads((val/'matlab/checks.json').read_text(encoding='utf-8'))
assert checks['complete'] and checks['cases']==75 and checks['rows']==300
rows=list(csv.DictReader((val/'matlab/metrics.csv').open(encoding='utf-8-sig')));main=summarize([r for r in rows if r['Suite']=='normal']);colored=summarize([r for r in rows if r['Suite']=='colored']);assert main['cases']==60 and colored['cases']==15
g=json.loads((ROOT/'core_p1/inputs/adaptive_cv_plan.json').read_text(encoding='utf-8'))['gates'];main['gates']=gates(main,g,checks)
colored['gates']={k:colored['means']['PrunedFullCalRefit'][k]<=colored['means']['DecoupledWeak'][k] for k in ['RelativeError','WeakEventError']};colored['gates']['solver']=checks['optimization_flags']==0;colored['gates']['all_pass']=all(colored['gates'].values())
api=run('pruned_api_run.json');api_check=json.loads((api/'matlab/checks.json').read_text(encoding='utf-8'));assert api_check['complete'] and api_check['api_replay_relative_error']<1e-8 and api_check['truth_input_rejected']
frozen=json.loads((ROOT/'core_p1/inputs/frozen_pruned_profile.json').read_text(encoding='utf-8'))
for r in [val,api]:
    source=r/'source/core_p1';manifest=json.loads((r/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    assert all(hashlib.sha256((source/p).read_bytes()).hexdigest()==v for p,v in manifest.items())
    assert all(hashlib.sha256((source/p).read_bytes()).hexdigest()==v for p,v in frozen['algorithm_sha256'].items())
assert all(hashlib.sha256((ROOT/'core_p1'/p).read_bytes()).hexdigest()==v for p,v in frozen['algorithm_sha256'].items())
ready=dev['main']['gates']['all_pass'] and dev['colored']['passes'] and main['gates']['all_pass'] and colored['gates']['all_pass']
result={'method':'PrunedFullCalRefit','ready_for_research_use':ready,'unique_cases':185,'development':dev['main'],'development_colored':dev['colored'],'validation':main,'colored_stress':colored,'api':api_check,'validation_checks':checks,'validation_run':str(val),'frozen_algorithm_unchanged':True,'scope':'8 geological models within one synthetic simulator; 100 normal and 10 colored development cases, 60 fresh normal validation cases and 15 fresh colored stress cases. Zero-phase Ricker dictionary, independent calibration, known geometry, 11 traces.'}
(ROOT/'analysis/final_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
status={'method':result['method'],'ready_for_research_use':ready,'research_only':True,'unique_cases':185,'gates':{'development':dev['main']['gates'],'validation':main['gates'],'colored_stress':colored['gates']},'api_replay_relative_error':api_check['api_replay_relative_error'],'scope':result['scope']}
(ROOT/'core_p1/inputs/validation_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 统一优化结果','','方法：完整校准＋反射系数估计＋反射组证据检验＋振幅重估（PrunedFullCalRefit）。用户入口为 `core_p1/compensate_seismic.m`，所有地质模型使用同一套规则；正则阈值与组筛选代价只随观测噪声、采样和校准主频变化。','', '| 资料 | 方法 | 整体误差 | 弱反射误差 | AVO梯度误差 | 尖灭尾部误差 |','|---|---|---:|---:|---:|---:|']
for title,s in [('开发100组',dev['main']),('全新模型60组',main),('全新有色噪声15组',colored)]:
    for m in ['FixedSplit','DecoupledWeak','PrunedFullCalRefit']:
        lines.append('| '+title+' | '+m+' | '+' | '.join(f'{s["means"][m][k]:.6g}' for k in ['RelativeError','WeakEventError','AVOGradientError','PinchoutTailError'])+' |')
lines+=['',f'统一研究版本验收：{"通过" if ready else "未通过"}。开发、新模型与有色噪声检查分别为 {dev["main"]["gates"]["all_pass"]}、{main["gates"]["all_pass"]}、{colored["gates"]["all_pass"]}。',f'整体误差低于现版：开发 {dev["main"]["case_wins"]}/100；全新模型 {main["case_wins"]}/60；全新有色噪声 {colored["case_wins"]}/15。','', '共185组不同观测：原60组开发＋第一轮50组资料已合并为开发（100组普通噪声、10组有色噪声）；第二轮为3个全新地质模型的60组普通噪声主验证及15组有色噪声压力测试。没有把第一轮用过的模型继续称为未见验证。','', '第一轮仅做振幅重估时，40组验证的尖灭尾部误差从0.001566升到0.002774，未通过门槛，完整记录见ROUND1_RESULTS_ZH.md。第二轮增加全时窗的反射组检验：删除一个组并重新拟合带来的残差增量，必须大于 σ²×(组内系数数+2)×log(nt)。该复杂度准则没有按模型调参，也没有使用尖灭位置；它不是具有已验证覆盖率的显著性检验。第二轮规则冻结后，没有根据新验证结果改动。','', '完整校准规则与现版Section域不同：同一份校准观测得到的完整估计直接用于所有采集域，现版Section仍使用几何收缩。反射系数域计算还加入零相位Ricker字典。收益来自整套机制，不能全部归于某一个正则参数或振幅重估步骤。','',f'公共API相对复现误差 {api_check["api_replay_relative_error"]:.3g}；真实Q输入拒绝测试通过。冻结算法与新验证/API运行的源码哈希一致。最大振幅重估正规方程残差 {checks["max_normal"]:.3g}；稀疏阶段未收敛标志 {checks["optimization_flags"]}。','', '适用范围：本研究仿真器内8个地质模型、给定采集几何、独立校准资料、11道输入、零相位Ricker子波假设。可作为该范围内统一调用的研究版本；尚未检验任意子波、密集反射或实测资料。弱证据组筛选有删除真实弱反射的可能，分模型和逐案例结果完整保留，不声称每一个弱事件都准确恢复。','', '道集图展示全部185组资料，振幅/误差分别使用公共色标，不做逐道归一化。普通与鲁棒CV的9组筛选、第一轮失败结果、计算中断及恢复记录均保留。']
(ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'ready':ready,'development':dev['main']['means'],'validation':main['means'],'validation_gates':main['gates'],'colored_gates':colored['gates']},ensure_ascii=False,indent=2))
