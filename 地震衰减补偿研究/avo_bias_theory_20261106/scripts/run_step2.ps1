# Runs step 2 for the three models in parallel, then the velocity stress group.
$py = 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$root = Split-Path $PSScriptRoot; $logs = Join-Path $root 'logs'; New-Item -ItemType Directory -Force $logs | Out-Null
$marker = Join-Path $root 'RUNNING.txt'
"Step 2 running, started $(Get-Date -Format s); do not clean this folder" | Out-File -Encoding ascii $marker
$env:OMP_NUM_THREADS = '4'; $env:OPENBLAS_NUM_THREADS = '4'; $env:MKL_NUM_THREADS = '4'
$p = foreach ($m in 'M1', 'M10', 'M11') {
  Start-Process -FilePath $py -ArgumentList 'step2.py', 'run', $m -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru `
    -RedirectStandardOutput (Join-Path $logs "step2_$m.log") -RedirectStandardError (Join-Path $logs "step2_$m.err")
}
$p | Wait-Process
$p | ForEach-Object { "exit $($_.Id) $($_.ExitCode)" }
$env:OMP_NUM_THREADS = '12'; $env:OPENBLAS_NUM_THREADS = '12'; $env:MKL_NUM_THREADS = '12'
$s = Start-Process -FilePath $py -ArgumentList 'step2.py', 'stress' -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru -Wait `
    -RedirectStandardOutput (Join-Path $logs 'step2_stress.log') -RedirectStandardError (Join-Path $logs 'step2_stress.err')
"stress exit $($s.ExitCode)"
Remove-Item $marker -Force -Confirm:$false
