from pathlib import Path
import datetime, json, subprocess, time

here = Path(__file__).resolve().parent
command = ['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe', '-wait', '-batch',
           "addpath('" + here.as_posix().replace("'", "''") + "'); review_smoke"]
start = time.monotonic()
status = {'command': command, 'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'timeout_seconds': 1200, 'execution': 'new independent hidden batch process'}
try:
    with (here / 'matlab_run2_stdout.log').open('w', encoding='utf-8') as stdout, (here / 'matlab_run2_stderr.log').open('w', encoding='utf-8') as stderr:
        proc = subprocess.Popen(command, stdout=stdout, stderr=stderr, creationflags=subprocess.CREATE_NO_WINDOW)
        status['pid'] = proc.pid
        (here / 'launch_status2.json').write_text(json.dumps(status, indent=2), encoding='utf-8')
        try:
            status['exit_code'] = proc.wait(timeout=1200)
            status['status'] = 'COMPLETED' if status['exit_code'] == 0 else 'FAILED'
        except subprocess.TimeoutExpired:
            subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'], capture_output=True)
            status['status'] = 'HARD_TIMEOUT'
finally:
    status['elapsed_seconds'] = time.monotonic() - start
    (here / 'launch_status2.json').write_text(json.dumps(status, indent=2), encoding='utf-8')
    print(json.dumps(status, indent=2))
