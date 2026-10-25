# Soft refit trial execution report

## Material Passport

- Origin Mode: run
- Verification Status: UNVERIFIED
- Experiment: frozen SoftGroupRefit 99-trace trial, 180 blocks
- Scope: `soft_refit_20261023`
- Scientific quality assessment: not performed by the execution agent

## Execution

- Exact command: `& 'C:/Users/qq175/Documents/ChatGPT/研究/soft_refit_20261023/launch_matlab.ps1'`
- Working directory: `C:/Users/qq175/Documents/ChatGPT/研究`
- Attempt: 2, explicitly instructed after the primary agent corrected launcher argument quoting
- Session: 87532
- MATLAB launcher PID: 8572
- Start: 2026-10-23 22:05:09 +08:00
- Exit recorded: 2026-10-23 22:07:45 +08:00
- Launcher/terminal exit code: 0
- Experiment-reported elapsed time: 136.7617492 seconds
- First observed launcher RSS: 9,560,064 bytes; this is the launcher process and does not measure MATLAB worker memory
- Final stdout: 12,368 bytes
- Final stderr: 0 bytes
- Timeout: 3,600 seconds; not reached

## Observed completion evidence

- Explicit stdout: “Soft variance propagation, Monte Carlo, zero-noise and dead-trace checks passed.”
- Completed blocks: 180 / 180
- Result JSON files: 12, consisting of 10 Q-scenario outputs and 2 references
- MAT checkpoint files: 180
- `checks.json`: `complete: true`
- `RUN_LEDGER.json`: `status: computed`
- Optimization flags: 0
- Maximum recomputed FullCalRefit versus frozen baseline relative difference: 0
- Maximum initial normal-equation residual: 9.528636571967633e-15
- Final progress: complex model, Q +30%, Strong noise, block 9 / 9
- No output stall, crash, timeout, or process termination observed

## Output inventory

| Path relative to experiment root | Bytes |
|---|---:|
| results/primary/reference.json | 4,099,250 |
| results/primary/Q_-30.json | 8,192,361 |
| results/primary/Q_-10.json | 8,234,226 |
| results/primary/Q_+0.json | 8,253,306 |
| results/primary/Q_+10.json | 8,260,496 |
| results/primary/Q_+30.json | 8,265,172 |
| results/complex/reference.json | 7,516,513 |
| results/complex/Q_-30.json | 14,956,244 |
| results/complex/Q_-10.json | 15,046,583 |
| results/complex/Q_+0.json | 15,083,627 |
| results/complex/Q_+10.json | 15,115,703 |
| results/complex/Q_+30.json | 15,142,569 |

Each Q-scenario JSON contains `error_percent`, `used_Q`, and two condition records. Each condition record contains `name`, `fixed`, `unpruned`, `hard`, `soft`, and nine block diagnostics. No method-quality comparison is made here.

## Prior launch attempt

Attempt 1 returned exit code 0 without executing the trial. Its report and empty logs were preserved by the primary agent as `EXECUTION_REPORT_ATTEMPT1.md`, `matlab_stdout_attempt1.log`, `matlab_stderr_attempt1.log`, and `matlab_exit_code_attempt1.txt`. The primary agent corrected only launcher argument quoting and explicitly instructed attempt 2. The execution agent did not modify scripts or parameters and did not retry without that instruction.

Completion and internal assertions establish execution evidence. Independent scientific/statistical validation remains for the reviewer; the Material Passport therefore remains UNVERIFIED.
