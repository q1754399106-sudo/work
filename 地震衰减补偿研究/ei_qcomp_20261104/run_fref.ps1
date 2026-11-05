# Runs the frozen reference-frequency plan one MATLAB process at a time; finished runs are skipped.
$ErrorActionPreference = 'Continue'
Set-Location $PSScriptRoot
$matlab = 'D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe'
$marker = Join-Path $PSScriptRoot 'RUNNING.txt'
"Reference-frequency runs in progress, started $(Get-Date -Format s); do not clean this folder" | Out-File -Encoding ascii $marker
Remove-Item Env:EI_SMOKE -ErrorAction SilentlyContinue
foreach ($seed in 1,2,3) {
  foreach ($q in '+0.00','+0.10','-0.10') {
    $tag = "Qerr${q}_snr30_seed${seed}_it4000_cosine_fref25"
    if (Test-Path "outputs\history_$tag.csv") { "skip $tag"; continue }
    $env:EI_FREF = '25'; $env:EI_SCEN = 'stationary'; $env:EI_QERR = $q; $env:EI_SEED = "$seed"; $env:EI_ITERS = '4000'; $env:EI_CFG = '1 2'; $env:EI_SCHED = 'cosine'
    & $matlab -batch "tic; ei_qcomp_net; fprintf('elapsed %.1f s\n', toc)" -logfile "outputs\run_$tag.log" | Out-Null
    "$tag exit=$LASTEXITCODE $(Get-Date -Format s)"
  }
}
Remove-Item $marker -Force -Confirm:$false
