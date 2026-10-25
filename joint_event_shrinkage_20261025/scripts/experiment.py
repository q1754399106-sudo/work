"""Stages are deliberately separated: checks -> dev selection -> frozen diagnosis."""
from common import *
from joint import prepare,solve,gradient_check,solve_trace
from scipy.ndimage import convolve1d
import argparse,time,csv,platform,scipy
DEV=['F1','F2'];DIAG=['F3','F4','V4','V5'];NOISES=['Moderate','Strong'];ETAS=[0.,.25,.5,1.]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def code_hashes():
    return {str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sorted((ROOT/'scripts').glob('*.py'))}
def input_hashes():
    return {str(p.relative_to(ROOT)).replace('\\','/'):dict(sha256=sha(p),bytes=p.stat().st_size) for base in [ROOT/'inputs',ROOT/'cache/inputs'] for p in sorted(base.rglob('*')) if p.is_file()}
def save_manifest():writejson(ROOT/'analysis/input_manifest.json',input_hashes())
def error(A,T):return float(np.linalg.norm(A-T)/max(np.linalg.norm(T),1e-30))
def case(name,noise):
    cm,D,sigma,st,tr=load_scene(name,noise)
    target=np.asarray(next(x for x in tr['tr'] if x['name']==noise)['target'])
    return cm,D,sigma,st,tr,target
def csvwrite(path,rows,fields):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def sc_output(key,target,stage):
    # Solved SC arrays were runtime cache; copy exact bytes into our immutable
    # inputs, retain converged metadata and cross-check committed score rows.
    base=REPO/'novelty_evidence_20261025/cache/solved';stem=f'{key}_0.125_10'
    for ext in ['.json','.npz']:
        dest=ROOT/'inputs/sc'/f'{stem}{ext}';src=base/f'{stem}{ext}'
        if not dest.exists():dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(src.read_bytes())
    info=json.loads((ROOT/'inputs/sc'/f'{stem}.json').read_text())
    assert info['converged'] and info['kkt_ratio']<=.02
    out=np.load(ROOT/'inputs/sc'/f'{stem}.npz')['out']
    rows=json.loads(pinned_bytes(f'novelty_evidence_20261025/analysis/{"dev" if stage=="dev" else "test"}_rows.json'))
    archived=next(x['error'] for x in rows if x['case']==key and x['method']=='SC_a0.125_m10')
    assert abs(error(out,target)-archived)<1e-12
    return out,info

def metrics(A,target,tr,cm,noise):
    ni=NOISES.index(noise);env=np.sqrt(convolve1d(target**2,np.ones(21),axis=0,mode='constant')/convolve1d(np.ones((len(target),1)),np.ones(21),axis=0,mode='constant'))
    quiet=env<.02*np.max(env)
    ret={'overall':error(A,target),'quiet_leak':float(np.linalg.norm(A[quiet])/np.linalg.norm(target))}
    for region,M in tr.get('regions',{}).items():
        M=np.asarray(M,bool)
        if M.ndim==3:M=M[:,:,ni]
        if M.shape==target.shape and np.any(M) and np.linalg.norm(target[M])>0:
            ret[region]=float(np.linalg.norm((A-target)[M])/np.linalg.norm(target[M]))
    if 'events' in tr:
        est=[];truth=[]
        for j in range(A.shape[1]):
            for e in vec(tr['weak'],int)-1:
                k=tr['events'][j,e]
                if not np.isfinite(k):continue
                k=int(k)-1
                if tr['reflectivity'][k,j]==0:continue
                ii=np.arange(max(0,k-10),min(len(A),k+11));z=np.pi*cm['fpeak']*(ii-k)*cm['dt'];w=(1-2*z*z)*np.exp(-z*z)
                est.append(float(A[ii,j]@w/(w@w)));truth.append(float(tr['reflectivity'][k,j]))
        if truth:
            est=np.asarray(est);truth=np.asarray(truth);ret['weak_projection_error']=float(np.linalg.norm(est-truth)/np.linalg.norm(truth));ret['weak_projection_gain']=float(est@truth/(truth@truth))
    return ret

def verify():
    cm,D,sigma,st,tr,target=case('F1','Moderate');p=prepare(cm,D,sigma,st)
    out,info=solve(p,0);difference=error(out,p['original']['Standard'])
    grad=max(gradient_check(t) for t in p['traces'] if t is not None)
    assert difference<1e-12 and grad<1e-7
    # All main eta values, both ablation modes, actual box KKT.
    checks=[]
    for mode in ['Joint','GOnly','NOnly']:
        for eta in ETAS:
            _,inf=solve(p,eta,mode);checks.append(dict(mode=mode,eta=eta,**inf))
    assert max(x['normal_residual'] for x in p['checks'])<1e-6
    synth=dict(E=np.array([2.,3.,1.]),V=np.array([.3,.2,.4]),lam=np.array([1.,4.,.5]))
    synth['G']=np.diag(synth['E']);synth['N']=np.diag(synth['V'])
    diagonal,_=solve_trace(synth,0.);full,_=solve_trace(synth,1.)
    orthogonal_difference=float(np.linalg.norm(diagonal-full));assert orthogonal_difference<1e-10
    output=dict(standard_eta0_relative_output_difference=difference,max_gradient_relative_error=grad,orthogonal_negative_control_gain_difference=orthogonal_difference,traces=p['checks'],solvers=checks)
    writejson(ROOT/'analysis/numerical_checks.json',output);save_manifest()
    writejson(ROOT/'analysis/environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),blas_threads=1))
    print('CHECKS PASS',difference,grad,max(x['variance_relative_error'] for x in p['checks']),flush=True)

def run_dev():
    assert (ROOT/'analysis/numerical_checks.json').exists()
    assert not (ROOT/'analysis/selection.json').exists(),'Selection exists; preserve frozen decision.'
    rows=[];allchecks={};timings=[]
    for name in DEV:
        for noise in NOISES:
            key=f'{name}_{noise}';print('DEV',key,flush=True)
            cm,D,sigma,st,tr,target=case(name,noise);p=prepare(cm,D,sigma,st);allchecks[key]=p['checks']
            rows.append(dict(case=key,method='Standard',eta=0.,error=error(p['original']['Standard'],target)))
            for eta in ETAS:
                out,info=solve(p,eta);er=error(out,target)
                rows.append(dict(case=key,method='Joint',eta=eta,error=er));timings.append(dict(case=key,method='Joint',eta=eta,preparation_seconds=p['preparation_seconds'],**info))
                print('eta',eta,'error',er,'seconds',info['seconds'],flush=True)
            writejson(ROOT/'analysis/dev_rows.json',rows);writejson(ROOT/'analysis/dev_checks.json',allchecks)
    candidates=[dict(eta=e,mean_error=float(np.mean([r['error'] for r in rows if r['method']=='Joint' and r['eta']==e]))) for e in ETAS]
    best=min(candidates,key=lambda x:(x['mean_error'],x['eta']));std=float(np.mean([r['error'] for r in rows if r['method']=='Standard']))
    save_manifest();selection=dict(**best,candidates=candidates,standard_mean=std,relative_improvement=1-best['mean_error']/std,
        selection_rule='Minimum arithmetic mean full-waveform relative error on F1/F2 x Moderate/Strong only; ties use smaller eta',
        frozen_before_diagnostic=True,plan_sha256=sha(ROOT/'PLAN_ZH.md'),code_hashes=code_hashes(),development_input_hashes=input_hashes(),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    writejson(ROOT/'analysis/dev_timings.json',timings);writejson(ROOT/'analysis/selection.json',selection)
    print('FROZEN SELECTION',best,'relative improvement',selection['relative_improvement'],flush=True)

def run_frozen():
    selection=json.loads((ROOT/'analysis/selection.json').read_text());assert selection['code_hashes']==code_hashes()
    assert selection['plan_sha256']==sha(ROOT/'PLAN_ZH.md')
    for rel,meta in selection['development_input_hashes'].items():assert sha(ROOT/rel)==meta['sha256']
    eta=selection['eta'];rows=[];regional=[];checks={};timings=[]
    # No eta grid is evaluated here, including on development ablations.
    for name in DEV+DIAG:
        for noise in NOISES:
            key=f'{name}_{noise}';stage='dev' if name in DEV else 'diagnostic';print('FROZEN',stage,key,'eta',eta,flush=True)
            cm,D,sigma,st,tr,target=case(name,noise);p=prepare(cm,D,sigma,st)
            methods={'Standard':p['original']['Standard'],'FullCalRefit':dense(st['S0'])}
            timings.append(dict(case=key,method='Standard',seconds=p['coherence_reference_seconds'],timing_scope='shared coherence plus four reference shrinkage rules, excludes dictionary/support discovery'))
            timings.append(dict(case=key,method='FullCalRefit',seconds=None,timing_scope='archived state; support discovery/refit not rerun'))
            timings.append(dict(case=key,method='SharedPreparation',seconds=p['preparation_seconds'],dictionary_seconds=p['dictionary_seconds'],coherence_reference_seconds=p['coherence_reference_seconds'],propagation_seconds=p['propagation_seconds']))
            archive=archived_mat(f'round3_20261024/results/{name}_Q+0_{noise}.mat')['gathers'];methods['StandardArchived']=archive['Standard']
            methods['SC'],scinfo=sc_output(key,target,stage)
            timings.append(dict(case=key,method='SC',seconds=scinfo['seconds'],timing_scope='archive solver, not rerun',max_kkt=scinfo['kkt_ratio']))
            for mode in ['Joint','GOnly','NOnly']:
                methods[mode],info=solve(p,eta,mode);timings.append(dict(case=key,method=mode,preparation_seconds=p['preparation_seconds'],**info))
            eta0,eta0info=solve(p,0.);diff=error(eta0,methods['Standard']);assert diff<1e-12
            checks[key]=dict(traces=p['checks'],eta0_standard_relative_output_difference=diff,
                full_reconstruction=error(p['W']@dense(st['R0']),dense(st['S0'])),soft_reconstruction=error(p['original']['Soft'],dense(st['soft'])),
                standard_archived_output_difference=error(methods['Standard'],methods['StandardArchived']))
            assert checks[key]['full_reconstruction']<1e-10 and checks[key]['soft_reconstruction']<1e-10
            path=ROOT/'outputs'/f'{key}.npz';path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**methods)
            for method,out in methods.items():
                er=error(out,target);rows.append(dict(case=key,stage=stage,method=method,error=er));print(method,er,flush=True)
                regional.extend(dict(case=key,stage=stage,method=method,metric=m,value=v) for m,v in metrics(out,target,tr,cm,noise).items())
            writejson(ROOT/'analysis/frozen_rows.json',rows);writejson(ROOT/'analysis/frozen_checks.json',checks);writejson(ROOT/'analysis/timings.json',timings)
    csvwrite(ROOT/'analysis/full_waveform.csv',rows,['case','stage','method','error']);csvwrite(ROOT/'analysis/regional.csv',regional,['case','stage','method','metric','value'])
    writejson(ROOT/'analysis/regional.json',regional);save_manifest()
    means={stage:{method:float(np.mean([r['error'] for r in rows if r['stage']==stage and r['method']==method])) for method in ['Standard','StandardArchived','FullCalRefit','SC','Joint','GOnly','NOnly']} for stage in ['dev','diagnostic']}
    dr=[r for r in rows if r['stage']=='diagnostic'];lookup={(r['case'],r['method']):r['error'] for r in dr}
    ratios={key:lookup[key,'Joint']/lookup[key,'Standard'] for key in sorted({r['case'] for r in dr})}
    weak={(r['case'],r['method']):r['value'] for r in regional if r['metric']=='weak_projection_error'}
    weakratios={f'V4_{noise}':weak[f'V4_{noise}','Joint']/weak[f'V4_{noise}','Standard'] for noise in NOISES}
    gates=dict(development_improvement_ge_2pct=selection['relative_improvement']>=.02,
        diagnostic_improvement_ge_3pct=1-means['diagnostic']['Joint']/means['diagnostic']['Standard']>=.03,
        worst_case_degradation_le_5pct=max(ratios.values())<=1.05,V4_weak_degradation_le_5pct=max(weakratios.values())<=1.05)
    summary=dict(means=means,diagnostic_relative_improvement=1-means['diagnostic']['Joint']/means['diagnostic']['Standard'],diagnostic_case_ratios=ratios,V4_weak_ratios=weakratios,
        gates=gates,continue_investment=all(gates.values()),SC_mean_ratio=means['diagnostic']['Joint']/means['diagnostic']['SC'],
        wins_over_SC=sum(lookup[key,'Joint']<lookup[key,'SC'] for key in ratios),diagnostic_cases=8,
        SC_archived_mean_check=abs(means['diagnostic']['SC']-.394409)<1e-6,selection_sha256=sha(ROOT/'analysis/selection.json'),source_hashes_unchanged=selection['code_hashes']==code_hashes())
    assert summary['SC_archived_mean_check'];writejson(ROOT/'analysis/summary.json',summary);print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['checks','dev','frozen']);args=a.parse_args()
    {'checks':verify,'dev':run_dev,'frozen':run_frozen}[args.stage]()
