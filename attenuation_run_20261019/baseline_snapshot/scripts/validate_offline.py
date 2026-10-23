"""Offline packaging/report checks; artificial acceptance fixtures are NOT research results."""
from pathlib import Path
import argparse, ast, csv, hashlib, importlib.util, json, math, re

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output=Path(output)
    if output.exists(): raise ValueError('Fresh output required')
    output.mkdir(parents=True)
    checks={}
    sources=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    checks['copied_source_hashes_match']=all(hashlib.sha256((ROOT/r['packaged']).read_bytes()).hexdigest()==r['sha256'] for r in sources)
    available=[r for r in sources if Path(r['source']).is_file()]
    checks['available_original_sources_unchanged']=all(hashlib.sha256(Path(r['source']).read_bytes()).hexdigest()==r['sha256'] for r in available)
    for p in (ROOT/'scripts').glob('*.py'): ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    checks['python_syntax']=True
    with (ROOT/'evidence/historical/metrics.csv').open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    lookup={(r['Model'],r['Scenario'],r['Seed'],r['Method'],r['Route']):r for r in rows if r['Suite']=='Main'}
    rechecks=[]
    examples=list((ROOT/'evidence/fourway/inputs/results/matlab_test').glob('example_*.json'))
    routes=[('Baseline','Direct',0,0),('Baseline','Split',1,0),('Frozen','Split',1,1),('Tuned','Split',1,2)]
    for p in examples:
        ex=json.loads(p.read_text(encoding='utf-8')); truth=ex['truth']; z=ex['Z']
        denominator=math.fsum(v*v for line in truth for v in line)
        for method,route,ir,im in routes:
            numerator=math.fsum((z[i][j][ir][im]-truth[i][j])**2 for i in range(len(truth)) for j in range(len(truth[0])))
            value=math.sqrt(numerator/denominator)
            expected=float(lookup[(ex['Model'],ex['Scenario'],'20261101',method,route)]['RelativeError'])
            rechecks.append(dict(file=p.name,method=method,route=route,difference=abs(value-expected)))
    assert len(rechecks)==24
    checks['stored_waveform_metrics_24']=max(x['difference'] for x in rechecks)<1e-12
    spec=importlib.util.spec_from_file_location('handoff_evaluator',ROOT/'scripts/evaluate_candidate.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    template=[dict(r) for r in rows if r['Method']=='Frozen']
    candidate=[]
    for r in template:
        r=dict(r); r['Method']='FixtureCandidate'
        for field in module.METRICS:
            if module.finite(r[field]): r[field]=str(float(r[field])*.9)
        candidate.append(r)
    def save(name,data):
        path=output/(name+'.csv')
        with path.open('w',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(data)
        return path
    def run(name,data):
        path=save(name,data)
        return module.evaluate(path,'FixtureCandidate','Frozen','Split',output/name)
    result=run('synthetic_pass',template+candidate)
    checks['acceptance_positive_fixture']=result['status']=='PROVISIONAL_PASS'
    result=run('missing_direct',template+[r for r in candidate if r['Route']=='Split'])
    checks['missing_direct_rejected']=result['status']=='FAIL'
    bad=[dict(r) for r in candidate]
    for r in bad:
        if r['Scenario']=='WeakSection': r['WeakEventError']='10'
    result=run('local_harm',template+bad)
    checks['global_gain_local_harm_rejected']=result['status']=='FAIL'
    # Input integrity must reject unmatched, duplicated and nonfinite cases.
    for name,kind in [('unpaired','unpaired'),('duplicate','duplicate'),('nonfinite','nonfinite')]:
        bad=[dict(r) for r in candidate]
        ix=next(i for i,r in enumerate(bad) if r['Route']=='Split')
        if kind=='unpaired': bad.pop(ix)
        elif kind=='duplicate': bad.append(dict(bad[ix]))
        else: bad[ix]['RelativeError']='NaN'
        try: run(name,template+bad)
        except ValueError: checks[name+'_rejected']=True
        else: checks[name+'_rejected']=False
    local_links=[]
    for name in ['START_HERE.md','AI_EXECUTION.md','TASKBOOK_ZH.md']:
        assert (ROOT/name).is_file()
    for p in (ROOT/'reports').glob('*.html'):
        for link in re.findall(r'(?:href|src)="([^"]+)"',p.read_text(encoding='utf-8')):
            if '://' not in link: local_links.append((p.parent/link).exists())
    checks['dashboard_links']=all(local_links) and bool(local_links)
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
            'copied_sources':len(sources),'available_originals_checked':len(available),
            'stored_waveform_count':len(rechecks),'waveform_max_absolute_delta':max(x['difference'] for x in rechecks),
            'scope':'Stored-waveform metrics and report logic only. No MATLAB inverse solve. Acceptance fixtures are artificial unit tests, not method evidence.'}
    (output/'OFFLINE_CHECKS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(result,indent=2,ensure_ascii=False))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=main(args.output)
    raise SystemExit(0 if result['status']=='PASS' else 1)
