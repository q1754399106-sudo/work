# Runs the frozen non-stationary reflectivity plan one MATLAB process at a time; finished runs are skipped.
$ErrorActionPreference = 'Continue'
Set-Location $PSScriptRoot
$matlab = 'D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe'
$marker = Join-Path $PSScriptRoot 'RUNNING.txt'
"Non-stationary reflectivity runs in progress, started $(Get-Date -Format s); do not clean this folder" | Out-File -Encoding ascii $marker
Remove-Item Env:EI_SMOKE -ErrorAction SilentlyContinue
foreach ($seed in 1,2,3) {
  foreach ($scen in 'density','amp','broad2narrow','narrow2broad') {
    $tag = "Qerr+0.00_snr30_seed${seed}_it4000_cosine_${scen}"
    if (Test-Path "outputs\history_$tag.csv") { "skip $tag"; continue }
    $env:EI_SCEN = $scen; $env:EI_QERR = '0'; $env:EI_SEED = "$seed"; $env:EI_ITERS = '4000'; $env:EI_CFG = '1 2'; $env:EI_SCHED = 'cosine'
    & $matlab -batch "tic; ei_qcomp_net; fprintf('elapsed %.1f s\n', toc)" -logfile "outputs\run_$tag.log" | Out-Null
    "$tag exit=$LASTEXITCODE $(Get-Date -Format s)"
  }
}
Remove-Item $marker -Force -Confirm:$false
