# Launch the step-1 chain (PLAN_ZH.md section 6) as detached hidden workers (run_step1_worker.ps1), one per
# argument; each worker runs va_run(<models>) in MATLAB at high priority and restarts it after a crash.
# Workers are created through Win32_Process.Create, so they do not belong to the calling shell's process
# tree or job and keep running after it ends.
#   powershell -NoProfile -ExecutionPolicy Bypass -File run_step1.ps1 D2 D4      # two processes in parallel
#   powershell -NoProfile -ExecutionPolicy Bypass -File run_step1.ps1 "D2,D4"    # one process, D2 then D4
# Logs: logs/run_<models>.log (MATLAB diary), logs/worker_<models>.log (attempts). Resumable.
$worker = Join-Path $PSScriptRoot 'run_step1_worker.ps1'
foreach ($a in $args) {
  $cl = 'powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File "' + $worker + '" "' + $a + '"'
  $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cl; CurrentDirectory = $PSScriptRoot }
  Write-Output ("worker for {0}: pid {1} (return {2})  {3}" -f $a, $r.ProcessId, $r.ReturnValue, (Get-Date -Format 'HH:mm:ss'))
}
