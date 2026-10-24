"""Run a hidden MATLAB process with provenance, logs, and a hard timeout."""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]
stage = sys.argv[1]
assert stage in {"export_inputs", "predict", "score_outputs"}
cmd = [r"D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe", "-wait", "-nojvm", "-batch",
       "addpath('%s');%s();" % ((ROOT / "scripts").as_posix().replace("'", "''"), stage)]
status = dict(stage=stage, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              command=cmd, timeout_seconds=2700, status="RUNNING")
path = ROOT / "logs" / (stage + "_execution.json")
start = time.monotonic()
with (ROOT / "logs" / (stage + ".log")).open("wb") as out, (ROOT / "logs" / (stage + "_err.log")).open("wb") as err:
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
if status["status"] != "COMPLETED":
    print((ROOT / "logs" / (stage + "_err.log")).read_text(errors="replace"))
    raise SystemExit(1)
if stage == "predict":
    paths = [ROOT / "PLAN_ZH.md", ROOT / "scripts" / "predict.m", ROOT / "scripts" / "export_inputs.m",
             ROOT / "scripts" / "score_outputs.m", ROOT / "scripts" / "run_stage.py",
             *sorted((ROOT / "inputs").glob("*.mat")), ROOT / "predictions" / "risk.csv",
             *sorted((ROOT / "analysis").glob("*.csv"))]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
    receipt = dict(sealed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   historical_labels_already_known=True, truth_used_for_predictions=False, hashes=hashes)
    (ROOT / "predictions" / "SEALED_SHA256.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
