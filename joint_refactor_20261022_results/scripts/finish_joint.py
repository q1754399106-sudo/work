"""Monitor owned run; checkpoint continuation only after its hard timeout."""
from pathlib import Path
import json,time,subprocess,sys,csv,collections,hashlib
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
    origin=Path(read(ROOT/'LATEST_RUN.json')['run']);start=time.monotonic();last=None
    while True:
        run=Path(read(ROOT/'LATEST_RUN.json')['run']);status=read(run/'execution.json');metrics=run/'matlab/metrics.csv'
        rows=list(csv.DictReader(metrics.open(encoding='utf-8-sig'))) if metrics.exists() else []
        if len(rows)!=last:
            print(json.dumps({'rows':len(rows),'status':status['status'],'run':run.name}),flush=True);last=len(rows)
        if status['status']!='RUNNING':break
        assert time.monotonic()-start<7500,'Controller exceeded observation budget; run left untouched'
        time.sleep(20)
    if status['status']=='HARD_TIMEOUT':
        # Every checkpoint is a full set of ten method/route cells.
        cells=collections.defaultdict(set)
        for r in rows:cells[tuple(r[k] for k in ['Model','Calibration','Scenario','Seed'])].add((r['Method'],r['Route']))
        assert all(len(v)==10 for v in cells.values())
        assert 0<len(rows)<640
        pointer={'metrics':str(metrics),'replay_error':read(run/'matlab/replay_check.json')['relative_difference']}
        (ROOT/'core_p1/inputs/resume_pointer.json').write_text(json.dumps(pointer,ensure_ascii=False,indent=2),encoding='utf-8')
        print('Continuing missing experimental cells after recorded hard timeout',flush=True)
        result=subprocess.run([sys.executable,str(ROOT/'scripts/run_local.py'),'--entry','run_joint_resume','--label','paired_continuation','--timeout','3600','--run-pointer',str(ROOT/'LATEST_RUN.json')])
        assert result.returncode==0,'Continuation failed; original and continuation failures retained'
        run=Path(read(ROOT/'LATEST_RUN.json')['run']);status=read(run/'execution.json')
    assert status['status']=='COMPLETED',status
    # Numerical components must match the first run in any continuation.
    for filename in ['joint_select.m','joint_physical_ensemble.m','p1_solve_branch.m']:
        a=origin/'source/core_p1/matlab_src'/filename;b=run/'source/core_p1/matlab_src'/filename
        assert a.read_bytes()==b.read_bytes(),filename
    subprocess.run([sys.executable,str(ROOT/'scripts/analyze_joint.py')],check=True)
    print('ALL EXPERIMENTAL CELLS AND ANALYSIS COMPLETE',flush=True)
if __name__=='__main__':main()
