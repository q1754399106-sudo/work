# Soft refit trial execution report

## Material Passport

- Origin Mode: run
- Verification Status: UNVERIFIED
- Experiment: frozen SoftGroupRefit 99-trace trial, expected 180 blocks
- Scope: `soft_refit_20261023`

## Attempt 1

- Exact command: `& 'C:/Users/qq175/Documents/ChatGPT/研究/soft_refit_20261023/launch_matlab.ps1'`
- Working directory: `C:/Users/qq175/Documents/ChatGPT/研究`
- Session: 1672
- MATLAB launcher PID: 7428
- Start: 2026-10-23 22:01:23 +08:00
- Process completion recorded: 2026-10-23 22:01:46 +08:00
- Launcher and terminal exit code: 0
- Standard output: 0 bytes
- Standard error: 0 bytes
- Result files: 0
- Checkpoint files: 0
- `progress.json`: absent
- `checks.json`: absent
- `RUN_LEDGER.json`: remains `prepared`, 0 / 180 blocks
- Baseline-match diagnostic: unavailable; no block executed
- Optimization flags: unavailable; no block executed
- Initial RSS: unavailable; launcher exited before first process inventory

## Execution status

The launcher returned exit code 0 without creating experiment outputs. This is an anomalous no-output completion and does **not** confirm mathematical tests or experiment execution. A likely cause is flattening the unquoted MATLAB `-batch` expression through PowerShell `Start-Process -ArgumentList`; the expression contains a space between `cd(...);` and `run_soft_trial`. This cause has not been verified by rerun.

No retry, source modification, parameter change, or process termination was performed. No scientific quality judgment is made in this report.
