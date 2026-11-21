# Launch the field chain (PLAN_ZH.md section 5) as one detached hidden worker (run_fe2_worker.ps1), which runs
# fe2_main in MATLAB at high priority and restarts it after a crash (resumable). The worker is created through
# Win32_Process.Create, so it does not belong to the calling shell's process tree and keeps running after it.
#   powershell -NoProfile -ExecutionPolicy Bypass -File run_fe2.ps1
# Logs: logs/run.log (MATLAB diary), logs/worker.log (attempts).
$worker = Join-Path $PSScriptRoot 'run_fe2_worker.ps1'
$cl = 'powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File "' + $worker + '"'
$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cl; CurrentDirectory = $PSScriptRoot }
Write-Output ("worker pid {0} (return {1})  {2}" -f $r.ProcessId, $r.ReturnValue, (Get-Date -Format 'HH:mm:ss'))
