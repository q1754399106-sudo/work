# Worker for run_step1.ps1: runs va_run(<models>) in a hidden MATLAB at high priority and restarts it
# (resumable, state every 5 frequencies) if it exits before analysis/va_step1_<model>.json exists for every
# model, at most 6 attempts. Usage: run_step1_worker.ps1 "D2"  or  "D2,D4".
param([string]$Models)
$matlab = 'D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe'
$root = $PSScriptRoot; $scripts = Join-Path $root 'scripts'
$names = ($Models -split ',') | ForEach-Object { $_.Trim() }
$wlog = Join-Path $root ("logs\worker_" + ($names -join '_') + ".log")
$cmd = "va_run({" + (($names | ForEach-Object { "'" + $_ + "'" }) -join ',') + "})"
function Test-Done { foreach ($n in $names) { if (-not (Test-Path (Join-Path $root "analysis\va_step1_$n.json"))) { return $false } }; return $true }
for ($try = 1; $try -le 6 -and -not (Test-Done); $try++) {
  $p = Start-Process -FilePath $matlab -ArgumentList @('-batch', "`"$cmd`"") -WorkingDirectory $scripts -WindowStyle Hidden -PassThru
  Add-Content $wlog ("{0} attempt {1}: {2}, launcher pid {3}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $try, $cmd, $p.Id)
  Start-Sleep -Seconds 30
  Get-CimInstance Win32_Process -Filter ("ParentProcessId={0}" -f $p.Id) | ForEach-Object {
    try { (Get-Process -Id $_.ProcessId).PriorityClass = 'High'; Add-Content $wlog ("  MATLAB pid {0} priority High" -f $_.ProcessId) } catch {} }
  $p.WaitForExit()
  Add-Content $wlog ("{0} attempt {1} exited, code {2}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $try, $p.ExitCode)
  if (-not (Test-Done)) { Start-Sleep -Seconds 60 }
}
Add-Content $wlog ("{0} worker end, done = {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), (Test-Done))
