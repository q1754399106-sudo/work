# Launch the step-1 chain (PLAN_ZH.md section 6) as detached, hidden MATLAB processes at high priority.
#   powershell -NoProfile -ExecutionPolicy Bypass -File run_step1.ps1 D2 D4      # two processes in parallel
#   powershell -NoProfile -ExecutionPolicy Bypass -File run_step1.ps1 "D2,D4"    # one process, D2 then D4
# Each process: va_run(<models>) in scripts/; log logs/run_<models>.log (MATLAB diary). Resumable.
$matlab = 'D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe'
$scripts = Join-Path $PSScriptRoot 'scripts'
$t0 = Get-Date
foreach ($a in $args) {
  $names = ($a -split ',') | ForEach-Object { "'" + $_.Trim() + "'" }
  $cmd = "va_run({" + ($names -join ',') + "})"
  $p = Start-Process -FilePath $matlab -ArgumentList @('-batch', "`"$cmd`"") -WorkingDirectory $scripts -WindowStyle Hidden -PassThru
  Write-Output ("started {0}  launcher pid {1}  {2}" -f $cmd, $p.Id, (Get-Date -Format 'HH:mm:ss'))
}
Start-Sleep -Seconds 30
Get-Process MATLAB -ErrorAction SilentlyContinue | Where-Object { $_.StartTime -ge $t0 } | ForEach-Object {
  $_.PriorityClass = 'High'; Write-Output ("MATLAB pid {0} priority High" -f $_.Id) }
