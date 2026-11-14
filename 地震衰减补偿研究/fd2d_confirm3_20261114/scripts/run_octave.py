"""Run one FD stage (make_shots or assemble_models) in Octave with logs and provenance (run_stage.py of fd2d_confirm_20261109
with MATLAB replaced by Octave + octave_compat_v2_20261025 (MATLAB rand / randn / rng) + scripts/octave_shim (discretize)).
Usage: run_octave.py make_shots M FIRST LAST | run_octave.py assemble_models"""
from pathlib import Path
import datetime, json, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]; COMPAT = ROOT.parent / "octave_compat_v2_20261025"
stage, args = sys.argv[1], sys.argv[2:]
assert stage in {"make_shots", "assemble_models"} and all(a.isdigit() for a in args)
tag = "_".join([stage, *args]); call = "%s(%s);" % (stage, ",".join(args))
cmd = ["nice", "-n", "10", "octave", "--no-gui", "--quiet", "--eval",
       "addpath('%s');addpath('%s');addpath('%s');%s" % (COMPAT.as_posix(), (ROOT / "scripts" / "octave_shim").as_posix(), (ROOT / "scripts").as_posix(), call)]
(ROOT / "logs").mkdir(exist_ok=True); path = ROOT / "logs" / (tag + "_execution.json")
status = dict(stage=stage, args=args, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), command=cmd, status="RUNNING")
start = time.monotonic()
with (ROOT / "logs" / (tag + ".log")).open("wb") as out, (ROOT / "logs" / (tag + "_err.log")).open("wb") as err:
    proc = subprocess.Popen(cmd, stdout=out, stderr=err); status["pid"] = proc.pid; path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    status["exit_code"] = proc.wait(); status["status"] = "COMPLETED" if proc.returncode == 0 else "FAILED"
status["elapsed_seconds"] = time.monotonic() - start; path.write_text(json.dumps(status, indent=2), encoding="utf-8")
print(json.dumps(status, ensure_ascii=True), flush=True)
raise SystemExit(0 if status["status"] == "COMPLETED" else 1)
