from pathlib import Path
import json,subprocess,time
ROOT=Path(__file__).resolve().parent;run=Path(json.loads((ROOT/'development_run.json').read_text())['run']);q=lambda p:str(p).replace('\\','/').replace("'","''")
dest=ROOT/'analysis/cache_replay_check.json';assert not dest.exists()
expr=f"addpath('{q(ROOT/'core')}');addpath('{q(ROOT/'core/matlab_src')}');verify_prepared_replay('{q(run/'Development/Moderate_prepared.mat')}','{q(run/'Development/Moderate.json')}','{q(dest)}');"
cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',expr];start=time.monotonic()
with (ROOT/'analysis/cache_replay_stdout.log').open('wb') as out,(ROOT/'analysis/cache_replay_stderr.log').open('wb') as err:
 p=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));code=p.wait(timeout=180)
(ROOT/'analysis/cache_replay_execution.json').write_text(json.dumps(dict(exit_code=code,pid=p.pid,elapsed_seconds=time.monotonic()-start),indent=2));print('Cache replay exit:',code);raise SystemExit(code)
