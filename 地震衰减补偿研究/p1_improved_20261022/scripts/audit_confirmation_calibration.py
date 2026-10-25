"""Summarize recorded calibrations without changing the frozen acceptance.

Three supplied-velocity conditions reuse each calibration observation. They are
not three independent calibration realizations. Shrinkage has no new optimizer.
"""
from pathlib import Path
import argparse,collections,hashlib,json,math,shutil,statistics
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def records(x):return x if isinstance(x,list) else [x] if x else []
def finite_positive(values):return all(math.isfinite(float(v)) and float(v)>0 for v in values)
def main():
 p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--output',type=Path);a=p.parse_args();run=a.run.resolve()
 assert read(run/'EXECUTION_SUMMARY.json')['status']=='COMPLETED','Wait for all numerical tasks'
 plan=read(run/'FROZEN_EXECUTION_PLAN.json');rows=[];expected=set();hashes={}
 for job in plan['jobs']:
  model=job['model']['name'];cohort='original_nine_regression' if job['reuse_original'] else 'prospective_synthetic'
  for seed in job['calibration_seeds']:
   for vg in sorted({float(s['vp_multiplier']) for s in job['suites']}):expected.add((model,int(seed),round(vg,6)))
  for file in sorted((run/model/'matlab').glob('calibration_*_v*.json')):
   d=read(file);cal=d['calibration'];old=cal['fixed_estimate'];base=d['estimate'];section=d['section_candidate_estimate'];layer=cal.get('layered',{})
   seed=int(file.stem.split('_')[1]);vg=statistics.mean(float(v)/float(t) for v,t in zip(d['provided_vp'],job['model']['vp']))
   corrected=statistics.mean(float(v)/float(t) for v,t in zip(d['corrected_vp'],job['model']['vp']))
   assert finite_positive(d['provided_vp']) and finite_positive(d['corrected_vp'])
   for est in [old,base,section]:assert finite_positive([est['effective_Q'],est['estimated_frequency'],*est.get('layer_Q',[])])
   row=dict(model=model,cohort=cohort,calibration_seed=seed,provided_vp_factor=round(vg,6),corrected_vp_factor_vs_synthetic_truth=corrected,
    velocity_correction_applied=bool(cal['applied']),proposed_velocity_boundary=bool(cal['velocity_boundary']),
    velocity_validation_gain=cal['validation_gain'],layered_applied=bool(layer.get('applied',False)),layered_validation_gain=layer.get('validation_gain'),
    joint_optimizer_flags=[int(r['exitflag']) for r in records(cal.get('optimization_runs'))],layered_optimizer_flags=[int(r['exitflag']) for r in records(layer.get('optimization_runs'))],
    original_optimizer_flag=int(old['optimizer_exitflag']),base_optimizer_flag=int(base['optimizer_exitflag']),
    original_q_boundary=bool(old['q_boundary']),original_f_boundary=bool(old['frequency_boundary']),
    base_q_boundary=bool(base['q_boundary']),base_f_boundary=bool(base['frequency_boundary']),section_q_boundary=bool(section['q_boundary']),section_f_boundary=bool(section['frequency_boundary']),
    section_optimizer_field_is_inherited_not_new_fit=True,section_uncertainty_scope=section.get('uncertainty_scope','Not a validated posterior coverage statement'),
    source_file=file.relative_to(run).as_posix())
   rows.append(row);hashes[row['source_file']]=hashlib.sha256(file.read_bytes()).hexdigest()
 keys=[(r['model'],r['calibration_seed'],r['provided_vp_factor']) for r in rows]
 assert set(keys)==expected and len(keys)==len(set(keys)),'Calibration file coverage or duplicate error'
 grouped=[]
 for cohort,vg in sorted({(r['cohort'],r['provided_vp_factor']) for r in rows}):
  rr=[r for r in rows if r['cohort']==cohort and r['provided_vp_factor']==vg];s=dict(cohort=cohort,provided_vp_factor=vg,records=len(rr))
  for field in ['velocity_correction_applied','proposed_velocity_boundary','layered_applied','original_q_boundary','original_f_boundary','base_q_boundary','base_f_boundary','section_q_boundary','section_f_boundary']:
   s[field+'_count']=sum(r[field] for r in rr)
  s['original_optimizer_nonpositive_count']=sum(r['original_optimizer_flag']<=0 for r in rr)
  s['base_optimizer_nonpositive_count']=sum(r['base_optimizer_flag']<=0 for r in rr)
  s['joint_run_flag_counts']=dict(collections.Counter(str(f) for r in rr for f in r['joint_optimizer_flags']))
  s['layered_run_flag_counts']=dict(collections.Counter(str(f) for r in rr for f in r['layered_optimizer_flags']))
  s['max_corrected_velocity_relative_error_vs_synthetic_truth']=max(abs(r['corrected_vp_factor_vs_synthetic_truth']-1) for r in rr)
  grouped.append(s)
 dest=a.output.resolve() if a.output else run/'calibration_audit';dest.mkdir(parents=True,exist_ok=False)
 result=dict(audit_complete=True,expected_records=len(expected),records=len(rows),independent_calibration_observations=len({(r['model'],r['calibration_seed']) for r in rows}),
  scope='Recorded diagnostics only; engineering acceptance is unchanged. Synthetic true velocity is used for post-run scoring only.',
  interpretation=['Nonpositive optimizer flags must be reported; a small inverse-solver residual does not establish calibration convergence.',
   'Section shrinkage is deterministic parameter combination; its optimizer flag is inherited and is not an independent fit or convergence certificate.',
   'Layerwise posterior coverage is unverified. Repeated velocity-error conditions reuse the same calibration observation.'],groups=grouped,rows=rows,source_sha256=hashes)
 save(dest/'CALIBRATION_AUDIT.json',result);shutil.copy2(Path(__file__),dest/'audit_confirmation_calibration.py')
 lines=['# 校准诊断审计','',f'覆盖 {len(rows)}/{len(expected)} 份校准记录，对应 {result["independent_calibration_observations"]} 组独立校准观测。三个速度输入条件复用同一校准观测，不能当成三次独立重复。','',
  '| 数据组 | 输入vp倍数 | 记录数 | 速度修正采用 | 分层校准采用 | 原优化标志≤0 | 新基础估计标志≤0 | 原Q/频率边界 | 新基础Q/频率边界 |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
 for s in grouped:
  lines.append(f'|{s["cohort"]}|{s["provided_vp_factor"]:.2f}|{s["records"]}|{s["velocity_correction_applied_count"]}|{s["layered_applied_count"]}|{s["original_optimizer_nonpositive_count"]}|{s["base_optimizer_nonpositive_count"]}|{s["original_q_boundary_count"]}/{s["original_f_boundary_count"]}|{s["base_q_boundary_count"]}/{s["base_f_boundary_count"]}|')
 lines+=['','非正优化停止标志与边界命中按实际记录披露；反演方程残差小不能证明校准已收敛。Section收缩没有单独的优化过程，其停止标志字段来自原估计，不应作为独立的收敛证据。','',
  '本审计不修改冻结验收门槛，也不把诊断记录齐全解释为方法合格。分层Q后验区间覆盖率尚未验证。校正速度与合成真值的比较只在运行完成后的评分阶段执行，没有传回方法或选参。',
  '', '完整逐记录结果、所有多起点优化标志、边界及实际使用的JSON文件哈希见 CALIBRATION_AUDIT.json。']
 (dest/'CALIBRATION_AUDIT_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');print(json.dumps({k:result[k] for k in ['audit_complete','records','independent_calibration_observations']}))
if __name__=='__main__':main()
