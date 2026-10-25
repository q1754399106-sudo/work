from pathlib import Path
import json,shutil,hashlib,datetime,csv
import numpy as np

ROOT=Path(__file__).resolve().parent;OLD=ROOT.parent/'q_identifiable_cluster_20261023'
for name in ['core','inputs','truth','analysis','figures','vendor','runs']:(ROOT/name).mkdir(exist_ok=True)
for p in (OLD/'core').rglob('*'):
    if not p.is_file():continue
    dest=ROOT/'core'/p.relative_to(OLD/'core');dest.parent.mkdir(parents=True,exist_ok=True)
    if not dest.exists():shutil.copy2(p,dest)
for name in ['Development.json']:shutil.copy2(OLD/'inputs'/name,ROOT/'inputs'/name)
shutil.copy2(OLD/'truth/Development.npz',ROOT/'truth/Development.npz')
shutil.copy2(OLD/'vendor/build_model_original.py',ROOT/'vendor/build_model_original.py')
shutil.copy2(OLD/'prepare_data.py',ROOT/'prepare_data.py')
spec=json.loads((OLD/'SPEC.json').read_text(encoding='utf-8'))
spec.update(registered_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),version='output_trust_v2',kappas=[.5,1.,2.],protection_snr_threshold=2.,
 methods=['FullCalRefit','PrunedFullCalRefit','QAwareRidge','QAwareCluster','OutputTrustNoQ','QAwareOutputTrust'],
 variants=[dict(name='V2DenseLowQ',seed=20264101,beds=70,spacing=[9,22],qbase=32,qgradient=.018,frequency=40,phase=18,offset=525,calibration_seed=20264102,observation_seeds=[20264111,20264112]),
 dict(name='V2BroadContrast',seed=20264201,beds=38,spacing=[23,37],qbase=65,qgradient=.015,frequency=55,phase=-8,offset=650,calibration_seed=20264202,observation_seeds=[20264211,20264212]),
 dict(name='V2ThinQReverse',seed=20264301,beds=85,spacing=[8,16],qbase=78,qgradient=-.018,frequency=32,phase=30,offset=375,calibration_seed=20264302,observation_seeds=[20264311,20264312])])
assert not (ROOT/'SPEC.json').exists();(ROOT/'SPEC.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
oldhash={p.relative_to(ROOT/'core').as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'core').rglob('*.m')}
(ROOT/'BASELINE_HASHES.json').write_text(json.dumps(oldhash,indent=2),encoding='utf-8')
rows=[]
for stage,cases in [('development',['Development']),('validation',[v['name'] for v in json.loads((OLD/'SPEC.json').read_text())['variants']])]:
    run=Path(json.loads((OLD/f'{stage}_run.json').read_text())['run'])
    for case in cases:
        truth=np.load(OLD/'truth'/f'{case}.npz')
        for condition in ['Moderate','Strong']:
            result=json.loads((run/case/(condition+'.json')).read_text())
            ref=np.array(result['QAwareRidge']);cur=np.array(result['QAwareCluster']);target=truth[condition+'_target']
            error=ref-target;delta=cur-ref;den=np.sum(target**2)
            cross=2*np.sum(error*delta)/den;delta2=np.sum(delta**2)/den
            actual=(np.sum((cur-target)**2)-np.sum(error**2))/den
            assert abs(cross+delta2-actual)<1e-12
            rows.append(dict(SourceStage=stage,Case=case,Condition=condition,ReferenceError=float(np.linalg.norm(error)/np.linalg.norm(target)),
                ClusterError=float(np.linalg.norm(cur-target)/np.linalg.norm(target)),OutputChange=float(np.linalg.norm(delta)/np.linalg.norm(ref)),CrossTerm=float(cross),ChangeNormTerm=float(delta2),MSEChange=float(actual)))
with (ROOT/'analysis/previous_round_output_diagnosis.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(dict(root=str(ROOT),baseline_files=len(oldhash),previous_conditions=len(rows),prior_MSE_increases=sum(r['MSEChange']>0 for r in rows)),indent=2))
