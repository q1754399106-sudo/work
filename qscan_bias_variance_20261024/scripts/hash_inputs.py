"""Record SHA-256 of upstream code/data and of this run's frozen plan and scripts."""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
BAE = REPO / "bae_q_20261024"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rel(p):
    return p.relative_to(REPO).as_posix()


head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
upstream = [*sorted((BAE / "core").glob("*.m")), *sorted((BAE / "source").glob("*.m")),
            *sorted((BAE / "inputs" / "independent").glob("*")), *sorted((BAE / "truth" / "independent").glob("*.mat")),
            *sorted((BAE / "results" / "ind").glob("F?_Q+0_*.mat")), *sorted((BAE / "results" / "ind_qscan").glob("F?_Q+0_*.mat")),
            BAE / "make_independent.m", BAE / "run_qscan.m", BAE / "SPEC.json",
            REPO / "theory_predict_diagnostic_20261024" / "predictions" / "risk.csv",
            REPO / "theory_predict_diagnostic_20261024" / "analysis" / "case_scores.csv",
            REPO / "theory_predict_diagnostic_20261024" / "analysis" / "scored_windows.csv"]
if sys.argv[1:] == ["upstream"]:
    out = dict(git_head=head, hashes={rel(p): sha(p) for p in upstream})
    (ROOT / "UPSTREAM_HASHES.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
elif sys.argv[1:] == ["prerun"]:
    frozen = [ROOT / "PLAN_ZH.md", ROOT / "UPSTREAM_HASHES.json", *sorted((ROOT / "scripts").glob("*.m")),
              *sorted((ROOT / "scripts").glob("*.py")), ROOT / "analysis" / "reproduction_check.csv"]
    out = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
               note="Recorded after stage 0 and before stages 1-2; noise-free and replicate results did not exist yet.",
               hashes={rel(p): sha(p) for p in frozen})
    (ROOT / "logs" / "PRE_RUN_HASHES.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
elif sys.argv[1:] == ["addendum"]:
    frozen = [ROOT / "PLAN_ZH.md", ROOT / "scripts" / "run_snr_probe.m", ROOT / "scripts" / "summarize_snr.m",
              ROOT / "scripts" / "run_stage.py"]
    out = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
               note="Addendum A registered while stage 2 was running and before it was summarized; no probe result existed.",
               hashes={rel(p): sha(p) for p in frozen})
    (ROOT / "logs" / "ADDENDUM_A_HASHES.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
else:
    raise SystemExit("usage: hash_inputs.py upstream|prerun|addendum")
print("written")
