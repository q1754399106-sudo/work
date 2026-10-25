from pathlib import Path
import shutil,json,hashlib,datetime
base=Path(__file__).resolve().parent
old=base/'adaptive_joint_cv_20261022'
root=base/'changed_models_comparison_20261023'
root.mkdir(exist_ok=False)
shutil.copytree(old/'core_p1',root/'core_p1')
(root/'scripts').mkdir()
shutil.copy2(old/'scripts/run_local.py',root/'scripts/run_local.py')
models=[
 dict(name='Changed_DeepLowQ',vp=[1820,2200,2620,3100],vs=[780,1080,1400,1730],rho=[1.94,2.12,2.33,2.54],Q=[22,32,45,75],event_t=[.43,.79,1.04],fpeak=50,nominal_frequency=50),
 dict(name='Changed_ThinPolarity',vp=[2050,2400,2650,3050],vs=[920,1200,1410,1700],rho=[2.05,2.28,1.98,2.50],Q=[55,35,90,60],event_t=[.33,.66,.672],fpeak=68,nominal_frequency=70),
 dict(name='Changed_LowContrast',vp=[2100,2180,2270,2400],vs=[980,1030,1090,1180],rho=[2.10,2.12,2.14,2.16],Q=[85,65,100,130],event_t=[.32,.62,.96],fpeak=38,nominal_frequency=40),
 dict(name='Changed_HighFreqQReverse',vp=[1850,2250,2660,3150],vs=[800,1110,1420,1770],rho=[1.96,2.16,2.35,2.55],Q=[38,85,30,100],event_t=[.36,.70,.732],fpeak=88,nominal_frequency=85)
]
plan=dict(origin='Third fresh geology comparison requested 2026-10-23. Registered before execution. Frozen solver and thresholds; no tuning on these cases.',models=models,scenarios=['SmoothAVO','CoherentNoiseAVO','FaultSection','WeakSection','PinchoutSection'],snr_db=[18,12],seeds=[20262411,20262412],calibration_seed=20262401,calibration_snr_db=18,expected_normal_cases=80,colored_stress=dict(kind='colored',snr_db=12,ar1_rho=.35,seed=20262413,expected_cases=20),comparison='Primary: original reconstructed FixedSplit vs frozen PrunedFullCalRefit; extra FixedSplitSharedCal uses the same full calibration as optimized. All methods receive identical observations. No significance claim; four independent geology configurations, multiple scenes/seeds per geology.')
(root/'core_p1/inputs/pruned_validation_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
(root/'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
src=(old/'core_p1/run_refit_validation.m').read_text(encoding='utf-8')
src=src.replace('function run_refit_validation(outdir,pruned)','function run_changed_models(outdir,pruned)')
src=src.replace('if nargin<2,pruned=false;end','if nargin<2,pruned=true;end')
src=src.replace('original=fixed_v1_operators(cm,offsets);','original=fixed_v1_operators(cm,offsets);\n  sharedOriginal=fixed_v1_operators(cmfull,offsets);')
src=src.replace('[Sfixed,~,~]=fixed_v1_apply(D,original);','[Sfixed,~,~]=fixed_v1_apply(D,original);\n    [Sshared,~,~]=fixed_v1_apply(D,sharedOriginal);')
src=src.replace('for method=1:(3+pruned)','for method=1:(4+pruned)')
src=src.replace("if method==4,label=methodName;S=Snew;normal=out.normal_residual;oflag=out.optimization_flags;end","if method==4,label=methodName;S=Snew;normal=out.normal_residual;oflag=out.optimization_flags;end\n     if method==5,label='FixedSplitSharedCal';S=Sshared;normal=NaN;oflag=0;end")
src=src.replace('fixed=Sfixed;observed=D;','fixed=Sfixed;fixed_shared=Sshared;observed=D;')
src=src.replace("'baseline','fixed','observed'","'baseline','fixed','fixed_shared','observed'")
src=src.replace("'fixed',fixed,'baseline'","'fixed',fixed,'fixed_shared',fixed_shared,'baseline'")
src=src.replace('count==expected*(3+pruned)','count==expected*(4+pruned)')
(root/'core_p1/run_changed_models.m').write_text(src,encoding='utf-8')
# Record and verify that every production MATLAB solver is unchanged.
manifest={}
for p in sorted((old/'core_p1').rglob('*.m')):
 rel=p.relative_to(old/'core_p1');a=hashlib.sha256(p.read_bytes()).hexdigest();b=hashlib.sha256((root/'core_p1'/rel).read_bytes()).hexdigest();assert a==b;manifest[rel.as_posix()]=a
(root/'FROZEN_SOLVER_HASHES.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
(root/'REGISTERED_UTC.txt').write_text(datetime.datetime.now(datetime.timezone.utc).isoformat(),encoding='utf-8')
print(root)
