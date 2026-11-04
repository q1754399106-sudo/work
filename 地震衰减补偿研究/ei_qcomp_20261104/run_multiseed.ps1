# Runs the frozen multi-seed plan one MATLAB process at a time; finished runs are skipped.
$ErrorActionPreference = 'Continue'
Set-Location $PSScriptRoot
$matlab = 'D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe'
$marker = Join-Path $PSScriptRoot 'RUNNING.txt'
New-Item -ItemType Directory -Force outputs | Out-Null
"multi-seed runs in progress, started $(Get-Date -Format s); do not clean this folder" | Out-File -Encoding ascii $marker
$runs = @(@{iters=4000; cfg='1 2 3'}, @{iters=12000; cfg='1 2'})
foreach ($r in $runs) {
  foreach ($seed in 1,2,3) {
    $tag = "Qerr+0.00_snr30_seed${seed}_it$($r.iters)_cosine"
    if (Test-Path "outputs\history_$tag.csv") { "skip $tag"; continue }
    $env:EI_SEED = "$seed"; $env:EI_ITERS = "$($r.iters)"; $env:EI_CFG = $r.cfg; $env:EI_SCHED = 'cosine'
    & $matlab -batch "tic; ei_qcomp_net; fprintf('elapsed %.1f s\n', toc)" -logfile "outputs\run_$tag.log" | Out-Null
    "$tag exit=$LASTEXITCODE $(Get-Date -Format s)"
  }
}
Remove-Item $marker -Force -Confirm:$false
