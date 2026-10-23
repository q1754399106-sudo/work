from pathlib import Path
import subprocess,json,time,sys
ROOT=Path(__file__).resolve().parent;q=lambda p:str(p).replace('\\','/').replace("'","''")
label=sys.argv[1] if len(sys.argv)>1 else 'conditional_covariance_check'
dest=ROOT/'analysis'/f'{label}.json';assert not dest.exists()
cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',f"addpath('{q(ROOT/'core')}');addpath('{q(ROOT/'core/matlab_src')}');check_output_covariance_mc('{q(dest)}');"]
start=time.monotonic()
with (ROOT/'analysis'/f'{label}_stdout.log').open('wb') as out,(ROOT/'analysis'/f'{label}_stderr.log').open('wb') as err:
 p=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));code=p.wait(timeout=180)
(ROOT/'analysis'/f'{label}_execution.json').write_text(json.dumps(dict(pid=p.pid,exit_code=code,elapsed_seconds=time.monotonic()-start),indent=2));print('Covariance check exit:',code);raise SystemExit(code)
