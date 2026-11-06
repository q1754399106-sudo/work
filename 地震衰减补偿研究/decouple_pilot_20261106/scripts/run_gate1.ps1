# Runs gate 1 for the three models in parallel; run with this folder as the current directory.
$py = 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$marker = Join-Path (Split-Path $PSScriptRoot) 'RUNNING.txt'
"Gate 1 running, started $(Get-Date -Format s); do not clean this folder" | Out-File -Encoding ascii $marker
$env:OMP_NUM_THREADS = '4'; $env:OPENBLAS_NUM_THREADS = '4'; $env:MKL_NUM_THREADS = '4'
$p = foreach ($m in 'M1', 'M2', 'M3') {
  Start-Process -FilePath $py -ArgumentList 'gate1.py', 'run', $m -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru `
    -RedirectStandardOutput (Join-Path $PSScriptRoot "gate1_$m.log") -RedirectStandardError (Join-Path $PSScriptRoot "gate1_$m.err")
}
$p | Wait-Process
$p | ForEach-Object { "exit $($_.Id) $($_.ExitCode)" }
Remove-Item $marker -Force -Confirm:$false
