"""Run one stage with logs and a hard timeout (after qscan_offset_20261025/scripts/run_stage.py).

  run_matlab.py make_scenes            MATLAB: scripts/make_scenes.m -> models/, truth/ (W1-W6, X1-X3)
  run_matlab.py stage_dev              Python: record model/truth SHA-256 (logs/MODEL_HASHES.json),
                                       copy models/*.mat into cache/runner/dev
  run_matlab.py run_combo3 W1 W2 ...  MATLAB: cache/runner/run_combo3.m on these models
"""
from pathlib import Path
import datetime, hashlib, json, shutil, subprocess, sys, time

sys.stdout.reconfigure(errors="replace")
ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "cache" / "runner"
MATLAB = r"D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe"
stage, args = sys.argv[1], sys.argv[2:]
assert stage in {"make_scenes", "stage_dev", "run_combo3"} and all(a.replace("_", "").isalnum() for a in args)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

if stage == "stage_dev":
    hashes = {p.relative_to(ROOT).as_posix(): sha(p) for d in ("models", "truth") for p in sorted((ROOT / d).glob("*.mat"))}
    assert len(hashes) == 18, len(hashes)
    out = ROOT / "logs" / "MODEL_HASHES.json"
    assert not out.exists(), "already staged"
    out.write_text(json.dumps(hashes, indent=2), encoding="utf-8")
    for d in ("dev", "cache", "analysis"):
        (RUNNER / d).mkdir(parents=True, exist_ok=True)
    for p in sorted((ROOT / "models").glob("*.mat")):
        shutil.copyfile(p, RUNNER / "dev" / p.name)
        assert sha(RUNNER / "dev" / p.name) == hashes[f"models/{p.name}"]
    print(f"staged {len(hashes) // 2} models")
    raise SystemExit(0)

if stage == "make_scenes":
    path, call = ROOT / "scripts", "make_scenes();"
else:
    hashes = json.loads((ROOT / "logs" / "MODEL_HASHES.json").read_text(encoding="utf-8"))
    for a in args:
        assert sha(RUNNER / "dev" / f"{a}.mat") == hashes[f"models/{a}.mat"], a
    path, call = RUNNER, "run_combo3({%s});" % ",".join("'%s'" % a for a in args)
tag = "_".join([stage, *args]) if len(args) < 3 else f"{stage}_{args[0]}_to_{args[-1]}"
cmd = [MATLAB, "-wait", "-nojvm", "-batch", "addpath('%s');%s" % (path.as_posix().replace("'", "''"), call)]
status = dict(stage=stage, args=args, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              command=cmd, timeout_seconds=8 * 3600, status="RUNNING")
log = ROOT / "logs"
log.mkdir(exist_ok=True)
spath = log / (tag + "_execution.json")
start = time.monotonic()
with (log / (tag + ".log")).open("wb") as out, (log / (tag + "_err.log")).open("wb") as err:
    proc = subprocess.Popen(cmd, stdout=out, stderr=err, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    status["pid"] = proc.pid
    spath.write_text(json.dumps(status, indent=2), encoding="utf-8")
    try:
        status["exit_code"] = proc.wait(timeout=status["timeout_seconds"])
        status["status"] = "COMPLETED" if proc.returncode == 0 else "FAILED"
    except subprocess.TimeoutExpired:
        subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
        status["status"] = "HARD_TIMEOUT"
status["elapsed_seconds"] = time.monotonic() - start
spath.write_text(json.dumps(status, indent=2), encoding="utf-8")
print(json.dumps(status, ensure_ascii=True), flush=True)
print((log / (tag + ".log")).read_text(errors="replace"))
if status["status"] != "COMPLETED":
    print((log / (tag + "_err.log")).read_text(errors="replace"))
    raise SystemExit(1)
