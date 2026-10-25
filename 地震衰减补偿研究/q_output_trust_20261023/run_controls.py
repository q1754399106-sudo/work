from pathlib import Path
import subprocess,json,time
ROOT=Path(__file__).resolve().parent;q=lambda p:str(p).replace('\\','/').replace("'","''")
dest=ROOT/'analysis/standalone_controls.json';assert not dest.exists()
cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',f"addpath('{q(ROOT/'core')}');addpath('{q(ROOT/'core/matlab_src')}');a=qtrust_controls();uq_writejson('{q(dest)}',a);"]
started=time.monotonic()
with (ROOT/'analysis/control_stdout.log').open('wb') as out,(ROOT/'analysis/control_stderr.log').open('wb') as err:
 p=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));code=p.wait(timeout=120)
(ROOT/'analysis/control_execution.json').write_text(json.dumps(dict(pid=p.pid,exit_code=code,elapsed_seconds=time.monotonic()-started),indent=2));print('Control exit:',code);raise SystemExit(code)
