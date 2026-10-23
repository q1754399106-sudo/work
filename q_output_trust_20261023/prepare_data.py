from pathlib import Path
import json, shutil, hashlib, datetime, sys
import numpy as np

ROOT=Path(__file__).resolve().parent
SPEC=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'))

def compact_case(name,payload,truth_path,calibration=None):
    idx=SPEC['indices_zero_based']
    payload=dict(payload)
    for d in payload['datasets']:
        d['D']=np.asarray(d['D'])[:,idx].tolist()
    payload['case_name']=name;payload['indices_zero_based']=idx
    payload['cached_calibration']=calibration
    (ROOT/'inputs'/f'{name}.json').write_text(json.dumps(payload,separators=(',',':')),encoding='utf-8')
    with np.load(truth_path) as t:
        arr={k:t[k] for k in ['t','x','live']+[d['name']+'_target' for d in payload['datasets']]}
        arr['x']=arr['x'][idx];arr['live']=arr['live'][idx]
        for d in payload['datasets']:arr[d['name']+'_target']=arr[d['name']+'_target'][:,idx]
    np.savez_compressed(ROOT/'truth'/f'{name}.npz',**arr)

def main():
    stage=sys.argv[1]
    if stage=='development':
        old=ROOT.parent/SPEC['development']
        payload=json.loads((old/'core_p1/inputs/complex_input.json').read_text(encoding='utf-8'))
        cal=json.loads((old/'runs/complex_frozen_comparison_20261023T032139Z/matlab/calibration.json').read_text(encoding='utf-8'))
        compact_case('Development',payload,old/'truth/model_truth.npz',cal)
        return
    assert stage=='validation'
    source=(ROOT/'vendor/build_model_original.py').read_text(encoding='utf-8')
    assert hashlib.sha256((ROOT/'vendor/build_model_original.py').read_bytes()).hexdigest()==SPEC['source_generator_sha256']
    for cfg in SPEC['variants']:
        dst=ROOT/'generated'/cfg['name']
        assert not dst.exists(),f'Refusing to overwrite generated case: {dst}'
        dst.mkdir(parents=True);(dst/'core_p1/inputs').mkdir(parents=True)
        (dst/'FROZEN_SOLVER_HASHES.json').write_text('{}',encoding='utf-8')
        text=source
        def rep(a,b):
            nonlocal text
            assert text.count(a)==1,(a,text.count(a))
            text=text.replace(a,b)
        rep("ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent",f"ROOT=Path({str(dst)!r});BASE=ROOT.parent")
        rep("OLD=BASE/'changed_models_comparison_20261023'","OLD=BASE/'unavailable_original_source'")
        rep("seed=20262501",f"seed={cfg['seed']}")
        rep("offset_m=450",f"offset_m={cfg['offset']}")
        rep("rng.uniform(12,30,59)",f"rng.uniform({cfg['spacing'][0]},{cfg['spacing'][1]},{cfg['beds']-1})")
        rep("35+.035*base_depth",f"{cfg['qbase']}+({cfg['qgradient']})*base_depth")
        rep("fpeak=36*(1+.07*np.sin(2*np.pi*xx/2940+.3));phase=24+7*np.sin(xx/560)",f"fpeak={cfg['frequency']}*(1+.07*np.sin(2*np.pi*xx/2940+.3));phase={cfg['phase']}+7*np.sin(xx/560)")
        rep("nominal_frequency=36",f"nominal_frequency={cfg['frequency']}")
        rep("default_rng(20262502)",f"default_rng({cfg['calibration_seed']})")
        rep("u=np.pi*36*(t[:,None]-ctime);s=((1-2*u*u)*np.exp(-u*u))@am;s=phase_source(s,24+7*np.sin(1470/560));",f"u=np.pi*{cfg['frequency']}*(t[:,None]-ctime);s=((1-2*u*u)*np.exp(-u*u))@am;s=phase_source(s,{cfg['phase']}+7*np.sin(1470/560));")
        rep("seed=20262511",f"seed={cfg['observation_seeds'][0]}")
        rep("seed=20262512",f"seed={cfg['observation_seeds'][1]}")
        (dst/'generated_source.py').write_text(text,encoding='utf-8')
        print('GENERATING',cfg['name'],flush=True)
        exec(compile(text,str(dst/'generated_source.py'),'exec'),{'__file__':str(dst/'generated_source.py'),'__name__':'__generated__'})
        correct_metadata(dst,cfg)
        (dst/'VARIANT.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
        payload=json.loads((dst/'core_p1/inputs/complex_input.json').read_text(encoding='utf-8'))
        compact_case(cfg['name'],payload,dst/'truth/model_truth.npz')
        print('INPUT READY',cfg['name'],flush=True)
    (ROOT/'inputs/validation_generated.json').write_text(json.dumps(dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases=[c['name'] for c in SPEC['variants']]),indent=2),encoding='utf-8')

def correct_metadata(dst,cfg):
    plan=json.loads((dst/'PLAN.json').read_text(encoding='utf-8'))
    plan.update(name=cfg['name'],geology=f"{cfg['beds']} irregular beds; inherited folding, faults, channel, wedge and low-Q lens",wavelet=f"Spatially varying {cfg['frequency']} Hz +/-7%; phase {cfg['phase']} +/-7 degrees; inherited source ghosts",calibration=f"Independent calibration seed {cfg['calibration_seed']}, nominal source {cfg['frequency']} Hz and midpoint phase; known coarse four-layer geometry")
    (dst/'PLAN.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    tp=dst/'truth/calibration_truth.json';truth=json.loads(tp.read_text(encoding='utf-8'))
    truth['source_phase_deg']=float(cfg['phase']+7*np.sin(1470/560));tp.write_text(json.dumps(truth,indent=2),encoding='utf-8')

if __name__=='__main__':main()
