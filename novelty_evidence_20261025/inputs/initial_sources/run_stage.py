"""Logged, bounded process; kills only its own process on declared timeout."""
from pathlib import Path
import sys,subprocess,time,json,datetime,hashlib
root=Path(__file__).resolve().parents[1]
stage=sys.argv[1]
assert stage in ('dev','test','theory')
cmd=[sys.executable,str(root/'scripts'/('check_theory.py' if stage=='theory' else 'benchmark.py'))]
if stage!='theory':cmd.append(stage)
paths=[root/'PLAN_ZH.md']+sorted((root/'scripts').glob('*.py'))
freeze={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
if stage=='test':freeze['analysis/selection.json']=hashlib.sha256((root/'analysis/selection.json').read_bytes()).hexdigest()
stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
meta=dict(command=cmd,start_utc=stamp,timeout_seconds=1800,input_hashes=freeze)
with (root/f'logs/{stage}.log').open('w',encoding='utf-8') as log:
    beg=time.monotonic();p=subprocess.Popen(cmd,cwd=root.parent,stdout=log,stderr=subprocess.STDOUT)
    meta['pid']=p.pid
    (root/f'logs/{stage}_execution.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    try:code=p.wait(timeout=1800)
    except subprocess.TimeoutExpired:
        p.kill();p.wait();code=124;meta['timed_out']=True
    meta.update(exit_code=code,seconds=time.monotonic()-beg,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    meta['source_unchanged']=all(hashlib.sha256((root/k).read_bytes()).hexdigest()==v for k,v in freeze.items())
    (root/f'logs/{stage}_execution.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in meta.items() if k!='input_hashes'}))
sys.exit(code)
