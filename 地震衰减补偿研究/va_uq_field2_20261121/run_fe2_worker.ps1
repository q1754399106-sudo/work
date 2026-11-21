# Worker for run_fe2.ps1: runs fe2_main in a hidden MATLAB at high priority and restarts it (resumable, state
# every 5 frequencies) if it exits before analysis/fe2_judge.json exists, at most 6 attempts.
$matlab = 'D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe'
$root = $PSScriptRoot; $scripts = Join-Path $root 'scripts'; $wlog = Join-Path $root 'logs\worker.log'
function Test-Done { return (Test-Path (Join-Path $root 'analysis\fe2_judge.json')) }
for ($try = 1; $try -le 6 -and -not (Test-Done); $try++) {
  $p = Start-Process -FilePath $matlab -ArgumentList @('-batch', 'fe2_main') -WorkingDirectory $scripts -WindowStyle Hidden -PassThru
  Add-Content $wlog ("{0} attempt {1}: fe2_main, launcher pid {2}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $try, $p.Id)
  Start-Sleep -Seconds 30
  Get-CimInstance Win32_Process -Filter ("ParentProcessId={0}" -f $p.Id) | ForEach-Object {
    try { (Get-Process -Id $_.ProcessId).PriorityClass = 'High'; Add-Content $wlog ("  MATLAB pid {0} priority High" -f $_.ProcessId) } catch {} }
  $p.WaitForExit()
  Add-Content $wlog ("{0} attempt {1} exited, code {2}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $try, $p.ExitCode)
  if (-not (Test-Done)) { Start-Sleep -Seconds 60 }
}
Add-Content $wlog ("{0} worker end, done = {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), (Test-Done))
