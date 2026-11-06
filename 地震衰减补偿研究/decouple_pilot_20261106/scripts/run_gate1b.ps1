# Runs gate 1b for the three models in parallel, then the mismatch stress group.
$py = 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$root = Split-Path $PSScriptRoot; $logs = Join-Path $root 'logs'; New-Item -ItemType Directory -Force $logs | Out-Null
$marker = Join-Path $root 'RUNNING.txt'
"Gate 1b running, started $(Get-Date -Format s); do not clean this folder" | Out-File -Encoding ascii $marker
$env:OMP_NUM_THREADS = '4'; $env:OPENBLAS_NUM_THREADS = '4'; $env:MKL_NUM_THREADS = '4'
$p = foreach ($m in 'M1', 'M2', 'M3') {
  Start-Process -FilePath $py -ArgumentList 'gate1b.py', 'run', $m -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru `
    -RedirectStandardOutput (Join-Path $logs "gate1b_$m.log") -RedirectStandardError (Join-Path $logs "gate1b_$m.err")
}
$p | Wait-Process
$p | ForEach-Object { "exit $($_.Id) $($_.ExitCode)" }
$env:OMP_NUM_THREADS = '8'; $env:OPENBLAS_NUM_THREADS = '8'; $env:MKL_NUM_THREADS = '8'
$s = Start-Process -FilePath $py -ArgumentList 'gate1b.py', 'stress' -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru -Wait `
    -RedirectStandardOutput (Join-Path $logs 'gate1b_stress.log') -RedirectStandardError (Join-Path $logs 'gate1b_stress.err')
"stress exit $($s.ExitCode)"
Remove-Item $marker -Force -Confirm:$false
