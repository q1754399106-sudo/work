from pathlib import Path
import subprocess,json,datetime,time
ROOT=Path(__file__).resolve().parent
dest=ROOT/'analysis/numerical_controls.json';assert not dest.exists()
q=lambda p:str(p).replace('\\','/').replace("'","''")
statement=f"addpath('{q(ROOT/'core')}');addpath('{q(ROOT/'core/matlab_src')}');c=qcluster_controls();uq_writejson('{q(dest)}',c);"
cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',statement]
start=time.monotonic()
with (ROOT/'analysis/controls_stdout.log').open('wb') as out,(ROOT/'analysis/controls_stderr.log').open('wb') as err:
    proc=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    status=dict(pid=proc.pid,command=cmd,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    code=proc.wait(timeout=120);status.update(exit_code=code,elapsed_seconds=time.monotonic()-start)
(ROOT/'analysis/controls_execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
print(json.dumps(status,indent=2));raise SystemExit(code)
