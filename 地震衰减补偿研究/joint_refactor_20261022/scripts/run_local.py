"""Run a fresh hidden MATLAB batch, snapshot exact source, retain every failure."""
from pathlib import Path
import argparse,datetime,hashlib,json,shutil,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--entry',default='run_calibration_study');p.add_argument('--label',default='calibration_development');p.add_argument('--timeout',type=float,default=3600);p.add_argument('--run-pointer',type=Path);a=p.parse_args()
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');run=ROOT/'runs'/(a.label+'_'+stamp);run.mkdir(parents=True)
 source=run/'source'/'core_p1';shutil.copytree(ROOT/'core_p1',source)
 if a.run_pointer:a.run_pointer.write_text(json.dumps({'run':str(run)},ensure_ascii=False,indent=2),encoding='utf-8')
 manifest={x.relative_to(source).as_posix():hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(source.rglob('*')) if x.is_file()}
 (run/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 q=lambda x:x.as_posix().replace("'","''")
 cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',f"addpath('{q(source)}');{a.entry}('{q(run/'matlab')}');"]
 status=dict(status='RUNNING',command=cmd,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),timeout_seconds=a.timeout)
 start=time.monotonic();print(str(run),flush=True)
 try:
  with (run/'stdout.log').open('wb') as out,(run/'stderr.log').open('wb') as err:
   proc=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));status['pid']=proc.pid
   (run/'execution.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
   try:
    status['exit_code']=proc.wait(timeout=a.timeout);status['status']='COMPLETED' if proc.returncode==0 else 'FAILED'
   except subprocess.TimeoutExpired:
    print(f'Hard timeout; terminating only owned MATLAB process tree {proc.pid}',flush=True)
    subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True);status['status']='HARD_TIMEOUT'
 finally:
  status['elapsed_seconds']=time.monotonic()-start
  (run/'execution.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
  print(json.dumps(status,ensure_ascii=True,indent=2),flush=True)
 return 0 if status['status']=='COMPLETED' else 1
if __name__=='__main__':raise SystemExit(main())
