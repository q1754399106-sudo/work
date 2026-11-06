$py = 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$root = Split-Path $PSScriptRoot; $logs = Join-Path $root 'logs'
$env:OMP_NUM_THREADS = '4'; $env:OPENBLAS_NUM_THREADS = '4'; $env:MKL_NUM_THREADS = '4'
$p = foreach ($m in 'M1', 'M10', 'M11') {
  Start-Process -FilePath $py -ArgumentList $args[0], $m -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru `
    -RedirectStandardOutput (Join-Path $logs "$($args[0])_$m.log") -RedirectStandardError (Join-Path $logs "$($args[0])_$m.err")
}
$p | Wait-Process
$p | ForEach-Object { "exit $($_.ExitCode)" }
