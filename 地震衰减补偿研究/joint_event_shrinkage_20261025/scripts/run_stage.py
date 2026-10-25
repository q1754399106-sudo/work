"""Run a stage with a permanent command/source/input audit and stdout log."""
from pathlib import Path
import sys,subprocess,time,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
stage=sys.argv[1];stamp=time.strftime('%Y%m%dT%H%M%SZ',time.gmtime());stem=f'{stage}_{stamp}'
cmd=[sys.executable,str(ROOT/'scripts/experiment.py'),stage]
audit=dict(command=cmd,cwd=str(ROOT),start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'scripts').glob('*.py')})
with (ROOT/'logs'/f'{stem}.log').open('w',encoding='utf-8') as f:
    proc=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
audit.update(exit_code=proc.returncode,end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
(ROOT/'logs'/f'{stem}.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print('Stage:',stage,'exit:',proc.returncode,'log:',str(ROOT/'logs'/f'{stem}.log'))
sys.exit(proc.returncode)
