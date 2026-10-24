"""Run one MATLAB stage in a hidden process with logs, provenance and a hard timeout.

Usage: run_stage.py <stage> [numeric args...], e.g. run_stage.py run_reps 1 1 50
"""
from pathlib import Path
import datetime, json, subprocess, sys, time

sys.stdout.reconfigure(errors="replace")   # MATLAB error logs can hold characters the console codec rejects
ROOT = Path(__file__).resolve().parents[1]
MATLAB = r"D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe"
STAGES = {"run_stage0", "run_reps", "summarize"}
stage, args = sys.argv[1], sys.argv[2:]
assert stage in STAGES and all(a.replace("_", "").isalnum() for a in args)
tag = "_".join([stage, *args])
call = "%s(%s);" % (stage, ",".join(a if a.isdigit() else "'%s'" % a for a in args))
cmd = [MATLAB, "-wait", "-nojvm", "-batch",
       "addpath('%s');%s" % ((ROOT / "scripts").as_posix().replace("'", "''"), call)]
status = dict(stage=stage, args=args, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              command=cmd, timeout_seconds=6 * 3600, status="RUNNING")
path = ROOT / "logs" / (tag + "_execution.json")
start = time.monotonic()
with (ROOT / "logs" / (tag + ".log")).open("wb") as out, (ROOT / "logs" / (tag + "_err.log")).open("wb") as err:
    proc = subprocess.Popen(cmd, stdout=out, stderr=err, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    status["pid"] = proc.pid
    path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    try:
        status["exit_code"] = proc.wait(timeout=status["timeout_seconds"])
        status["status"] = "COMPLETED" if proc.returncode == 0 else "FAILED"
    except subprocess.TimeoutExpired:
        print(f"Hard timeout: terminating owned MATLAB process tree {proc.pid}", flush=True)
        subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
        status["status"] = "HARD_TIMEOUT"
status["elapsed_seconds"] = time.monotonic() - start
path.write_text(json.dumps(status, indent=2), encoding="utf-8")
print(json.dumps(status, ensure_ascii=True), flush=True)
print((ROOT / "logs" / (tag + ".log")).read_text(errors="replace"))
if status["status"] != "COMPLETED":
    print((ROOT / "logs" / (tag + "_err.log")).read_text(errors="replace"))
    raise SystemExit(1)
