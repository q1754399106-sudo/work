from pathlib import Path
import json, csv, hashlib, sys, datetime
import numpy as np

ROOT=Path(__file__).resolve().parent
SPEC=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'))

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def writecsv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def masks(target,t,live):
    power=np.apply_along_axis(lambda a:np.convolve(a*a,np.ones(15)/15,mode='same'),0,target)
    local=np.sqrt(power);region=(t[:,None]>=.3)&(t[:,None]<=1.6)&live[None,:]
    nonzero=region&(local>.001*local.max());lo,hi=np.quantile(local[nonzero],[.20,.45])
    return region&(local>=lo)&(local<=hi)

def score(s,target,t,live,low):
    e=s-target;deep=(t>=1.05)&(t<=1.6);a=target[:,live].ravel();b=s[:,live].ravel()
    ratio=lambda u,v:float(np.linalg.norm(u)/max(np.linalg.norm(v),1e-30))
    am=a-a.mean();bm=b-b.mean()
    return dict(RelativeError=ratio(e,target),LiveTraceError=ratio(e[:,live],target[:,live]),
        DeepLiveError=ratio(e[deep,:][:,live],target[deep,:][:,live]),LowEnergyRegionError=ratio(e[low],target[low]),
        CorrelationLive=float(np.dot(am,bm)/max(np.linalg.norm(am)*np.linalg.norm(bm),1e-30)),
        AmplitudeProjectionGain=float(np.dot(a,b)/max(np.dot(a,a),1e-30)),
        PreArrivalErrorRMS=float(np.sqrt(np.mean(e[t<.25,:][:,live]**2))))

def main():
    stage=sys.argv[1];pointer=json.loads((ROOT/f'{stage}_run.json').read_text(encoding='utf-8'));run=Path(pointer['run'])
    execution=json.loads((run/'execution.json').read_text(encoding='utf-8'))
    assert execution['status']=='COMPLETED' and execution['exit_code']==0,execution
    manifest=json.loads((run/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    assert all(digest(run/'source'/p)==h for p,h in manifest.items()),'Run source changed'
    rows=[];trace_rows=[];mode_rows=[];checks=[];provenance=[]
    for case in execution['cases']:
        data=json.loads((ROOT/'inputs'/f'{case}.json').read_text(encoding='utf-8'))
        truth=np.load(ROOT/'truth'/f'{case}.npz');t=truth['t'];live=truth['live'].astype(bool)
        ck=json.loads((run/case/'checks.json').read_text(encoding='utf-8'));assert ck['complete'];checks.append(ck)
        cal=json.loads((run/case/'calibration.json').read_text(encoding='utf-8'))
        for cfg in data['datasets']:
            path=run/case/(cfg['name']+'.json');result=json.loads(path.read_text(encoding='utf-8'))
            target=truth[cfg['name']+'_target'];low=masks(target,t,live)
            for method in ['Observed']+SPEC['methods']:
                s=np.array(result[method]);assert s.shape==target.shape and np.isfinite(s).all()
                row=dict(Stage=stage,Case=case,Condition=cfg['name'],Method=method,**score(s,target,t,live,low))
                row['ObservedReplayResidual']=float(np.linalg.norm(s-np.array(result['Observed']))/max(np.linalg.norm(result['Observed']),1e-30))
                rows.append(row)
                for j,idx in enumerate(SPEC['indices_zero_based']):
                    denom=np.linalg.norm(target[:,j]);trace_rows.append(dict(Stage=stage,Case=case,Condition=cfg['name'],TraceOriginal=idx+1,Method=method,RelativeError=float(np.linalg.norm(s[:,j]-target[:,j])/max(denom,1e-30))))
            for method,stats in result['diagnostics'].items():
                for j,item in enumerate(stats):
                    pathvals=np.atleast_1d(item['objective_path']).astype(float)
                    assert np.all(np.diff(pathvals)<=1e-9*max(pathvals[0],1e-30))
                    assert item['rank_after']<=item['rank_before'] and item['normal_residual']<1e-8
                    mode_rows.append(dict(Stage=stage,Case=case,Condition=cfg['name'],TraceOriginal=SPEC['indices_zero_based'][j]+1,Method=method,
                        RankBefore=item['rank_before'],RankAfter=item['rank_after'],Clusters=item['clusters'],
                        ModelErrorToNoise=item['model_error_to_noise'],NormalResidual=item['normal_residual'],
                        ObjectiveBefore=float(pathvals[0]),ObjectiveAfter=float(pathvals[-1]),Steps=len(pathvals)-1))
            provenance.append(dict(case=case,condition=cfg['name'],result=str(path),sha256=digest(path),truth_sha256=digest(ROOT/'truth'/f'{case}.npz')))
    writecsv(ROOT/'analysis'/f'{stage}_metrics.csv',rows)
    writecsv(ROOT/'analysis'/f'{stage}_trace_metrics.csv',trace_rows)
    writecsv(ROOT/'analysis'/f'{stage}_modes.csv',mode_rows)
    agg={}
    for method in ['Observed']+SPEC['methods']:
        subset=[r for r in rows if r['Method']==method]
        agg[method]={k:float(np.mean([r[k] for r in subset])) for k in ['LiveTraceError','LowEnergyRegionError','DeepLiveError','CorrelationLive','AmplitudeProjectionGain']}
    gates={}
    if stage=='validation':
        candidate=agg['QAwareCluster'];full=agg['FullCalRefit']
        gates['overall_improves_2pct']=candidate['LiveTraceError']<=.98*full['LiveTraceError']
        gates['each_model_within_2pct']=all(np.mean([r['LiveTraceError'] for r in rows if r['Case']==case and r['Method']=='QAwareCluster'])<=1.02*np.mean([r['LiveTraceError'] for r in rows if r['Case']==case and r['Method']=='FullCalRefit']) for case in execution['cases'])
        gates['better_than_noQ_1pct']=candidate['LiveTraceError']<=.99*agg['ClusterModesNoQ']['LiveTraceError']
        gates['better_than_Qridge_1pct']=candidate['LiveTraceError']<=.99*agg['QAwareRidge']['LiveTraceError']
        gates['low_energy_within_2pct']=candidate['LowEnergyRegionError']<=1.02*full['LowEnergyRegionError']
        gates['numerical']=all(c['optimization_flags']==0 and c['max_normal']<1e-8 and c['zero_Q_covariance_limit']<1e-10 and c['zero_input_exact'] for c in checks)
        gates['all_pass']=all(gates.values())
    out=dict(stage=stage,verification='ANALYZED',run=str(run),means=agg,gates=gates,checks=checks,provenance=provenance,
        scope='One development geology or three new same-family geological configurations; two dependent conditions per geology; descriptive equal-geology means, no significance claims.')
    (ROOT/'analysis'/f'{stage}_summary.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    lines=[f'# {stage} 结果','', '## Material Passport','', '- Origin Skill: academic-research-suite / experiment-agent','- Origin Mode: validate',f'- Origin Date: {datetime.date.today().isoformat()}','- Verification Status: ANALYZED',f'- Version Label: q_cluster_{stage}_v1','',
        '| 模型 | 条件 | 方法 | 正常道误差 | 深部误差 | 低能量区误差 |','|---|---|---|---:|---:|---:|']
    for r in rows:lines.append(f"| {r['Case']} | {r['Condition']} | {r['Method']} | {r['LiveTraceError']:.6f} | {r['DeepLiveError']:.6f} | {r['LowEnergyRegionError']:.6f} |")
    lines+=['','## 等模型权重汇总','','| 方法 | 正常道误差 | 低能量区误差 | 深部误差 |','|---|---:|---:|---:|']
    for method,r in agg.items():lines.append(f"| {method} | {r['LiveTraceError']:.6f} | {r['LowEnergyRegionError']:.6f} | {r['DeepLiveError']:.6f} |")
    lines+=['','## 预登记工程门槛','',json.dumps(gates,ensure_ascii=False,indent=2),'',
        '均值与胜负是描述统计；道、同模型噪声条件不视为独立地质模型。求解器通过不等于物理模型或保幅成立。低能量区评价使用真值，只在评分程序内读取。',
        '11/11类统计谬误按本实验范围检查：分模型呈现以避免聚合掩盖退化；不从区域均值推断单反射；生成地质的选择局限、支持选择和同族验证均披露；不按结果筛选/控制案例；诊断基率不适用；不将开发到验证变化归因于重复测量改善；保留全部失败；不做多重探索的显著性宣称；规则与源码冻结；归因限于配对消融；受控计算不涉及观察性反向因果。',
        '',f'原始运行：{run}；全部CSV、模式路径和源码哈希随结果保存。']
    (ROOT/'analysis'/f'{stage}_RESULTS_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(dict(stage=stage,means=agg,gates=gates),indent=2),flush=True)

if __name__=='__main__':main()
