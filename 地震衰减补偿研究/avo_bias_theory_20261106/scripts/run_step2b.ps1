# Runs step 2b for the four new models in parallel.
$py = 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$root = Split-Path $PSScriptRoot; $logs = Join-Path $root 'logs'
$env:OMP_NUM_THREADS = '3'; $env:OPENBLAS_NUM_THREADS = '3'; $env:MKL_NUM_THREADS = '3'
$p = foreach ($m in 'M12', 'M13', 'M14', 'M15') {
  Start-Process -FilePath $py -ArgumentList 'step2b.py', 'run', $m -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru `
    -RedirectStandardOutput (Join-Path $logs "step2b_$m.log") -RedirectStandardError (Join-Path $logs "step2b_$m.err")
}
$p | Wait-Process
$p | ForEach-Object { "exit $($_.ExitCode)" }
