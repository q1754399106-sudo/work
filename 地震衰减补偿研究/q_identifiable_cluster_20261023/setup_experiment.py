from pathlib import Path
import hashlib, json, shutil, datetime

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
OLD=BASE/'complex_fieldlike_20261023'
for folder in ['core','inputs','truth','analysis','figures','vendor','runs']:
    (ROOT/folder).mkdir(exist_ok=True)
for p in (OLD/'core_p1').rglob('*.m'):
    dest=ROOT/'core'/p.relative_to(OLD/'core_p1')
    dest.parent.mkdir(parents=True,exist_ok=True)
    if not dest.exists(): shutil.copy2(p,dest)
(ROOT/'core/inputs').mkdir(exist_ok=True)
for name in ['candidate_policy.json','p1_experiment_spec.json','uq_spec.json','v2_options.json','options.json']:
    p=OLD/'core_p1/inputs'/name
    if p.exists() and not (ROOT/'core/inputs'/name).exists():shutil.copy2(p,ROOT/'core/inputs'/name)
shutil.copy2(OLD/'scripts/build_model.py',ROOT/'vendor/build_model_original.py')
hashes={p.relative_to(ROOT/'core').as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (OLD/'core_p1').rglob('*.m') for p in [ROOT/'core'/p.relative_to(OLD/'core_p1')]}
(ROOT/'FROZEN_BASELINE_HASHES.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
indices=[round(k*98/21) for k in range(22)]
spec=dict(registered_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),indices_zero_based=indices,
 development='complex_fieldlike_20261023',methods=['FullCalRefit','PrunedFullCalRefit','GlobalTSVD','ClusterModesNoQ','QAwareRidge','QAwareCluster'],
 finite_difference_log_step=.02,mode_penalty_log_nt_multiplier=1.,svd_relative_cutoff=.01,
 covariance_log_std_floor=.03,covariance_log_std_cap=.5,calibration_log_prior_std=.5,
 variants=[dict(name='FreshFoldLowQ',seed=20262701,beds=42,spacing=[18,33],qbase=28,qgradient=.025,frequency=42,phase=12,offset=600,calibration_seed=20262702,observation_seeds=[20262711,20262712]),
 dict(name='FreshDenseQReverse',seed=20262801,beds=80,spacing=[8,18],qbase=85,qgradient=-.025,frequency=48,phase=-16,offset=450,calibration_seed=20262802,observation_seeds=[20262811,20262812]),
 dict(name='FreshThinLowFrequency',seed=20262901,beds=50,spacing=[14,29],qbase=45,qgradient=.025,frequency=30,phase=35,offset=300,calibration_seed=20262902,observation_seeds=[20262911,20262912])],
 source_generator_sha256=hashlib.sha256((ROOT/'vendor/build_model_original.py').read_bytes()).hexdigest(),
 note='22 fixed spatial samples per 99-trace model; two observations are dependent conditions, not independent models.')
path=ROOT/'SPEC.json'
if path.exists():
    prev=json.loads(path.read_text(encoding='utf-8'));spec['registered_utc']=prev['registered_utc'];assert prev==spec
else:path.write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(root=str(ROOT),baseline_files=len(hashes),indices=indices),ensure_ascii=False),flush=True)
