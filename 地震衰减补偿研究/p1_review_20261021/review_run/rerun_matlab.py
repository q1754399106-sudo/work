from pathlib import Path
import csv, datetime, json, subprocess, time

here = Path(__file__).resolve().parent
stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
run = here / ('rerun_' + stamp)
run.mkdir()
quote = lambda p: p.as_posix().replace("'", "''")
command = ['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe', '-wait', '-batch',
           f"addpath('{quote(here)}'); review_smoke('{quote(run / 'results')}')"]
start = time.monotonic()
status = {'command':command, 'timeout_seconds':1200, 'output':str(run), 'status':'RUNNING'}
try:
    with (run / 'stdout.log').open('wb') as out, (run / 'stderr.log').open('wb') as err:
        proc = subprocess.Popen(command, stdout=out, stderr=err, creationflags=subprocess.CREATE_NO_WINDOW)
        status['pid'] = proc.pid
        (run / 'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
        print(str(run),flush=True)
        try:
            status['exit_code'] = proc.wait(timeout=1200)
            status['status'] = 'COMPLETED' if status['exit_code']==0 else 'FAILED'
        except subprocess.TimeoutExpired:
            subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True)
            status['status'] = 'HARD_TIMEOUT'
    if status['status']=='COMPLETED':
        def read(p):
            with p.open(encoding='utf-8-sig',newline='') as f:
                return list(csv.DictReader(f))
        old = read(here / 'matlab_smoke_corrected/metrics.csv')
        new = read(run / 'results/metrics.csv')
        key = lambda r: tuple(r[k] for k in ('Model','Suite','Scenario','Seed','CalibrationSeed','Method','Route'))
        oldmap = {key(r):r for r in old}
        differences = []
        for r in new:
            previous = oldmap.get(key(r))
            if previous != r:
                differences.append({'key':key(r),'changed_fields':[k for k in r if previous is None or previous.get(k)!=r[k]]})
        paired = {}
        for r in new:
            if r['Route']=='Split' and r['Method'] in ('Frozen','CombinedUncertainty'):
                paired.setdefault(r['Scenario'],{})[r['Method']]=float(r['RelativeError'])
        b = sum(r['Frozen'] for r in paired.values())/len(paired)
        c = sum(r['CombinedUncertainty'] for r in paired.values())/len(paired)
        comparison = {'old_rows':len(old),'new_rows':len(new),'identical_rows':not differences and len(old)==len(new),
                      'differences':differences,'scenarios':paired,'frozen_mean':b,'combined_mean':c,
                      'relative_improvement':1-c/b,'solver_violations':[r for r in new if float(r['NormalResidual'])>2e-7 or float(r['FinalFlag'])!=0]}
        (run / 'comparison.json').write_text(json.dumps(comparison,ensure_ascii=False,indent=2),encoding='utf-8')
finally:
    status['elapsed_seconds'] = time.monotonic()-start
    (run / 'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
    print(json.dumps(status,indent=2),flush=True)
