$d = Split-Path (Split-Path $PSScriptRoot); Set-Location $PSScriptRoot
$py = 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'; $m = 'D:\Program Files (x86)\MATLAB2024B\bin\matlab.exe'
$ex = "$d\abt\ex"; $src = "$d\base\compensation_package_20261024\source"; $src2 = "$d\base\attenuation_research_handoff_20261019\core\matlab_src"; $inp = "$d\base\attenuation_research_handoff_20261019\core\inputs"
& $py dev_base2.py export
& $m -batch "addpath('$PSScriptRoot');base_fcr('$ex','$src','dev2_M1_*')"
& $m -batch "addpath('$PSScriptRoot');base_frozen('$ex','$src2','$inp','dev2_M1_*')"
$env:OMP_NUM_THREADS = '4'; $env:OPENBLAS_NUM_THREADS = '4'; $env:MKL_NUM_THREADS = '4'
$p = foreach ($i in 0, 1) { Start-Process -FilePath $py -ArgumentList 'dev_base2.py', 'sc', $i -WorkingDirectory $PSScriptRoot -NoNewWindow -PassThru -RedirectStandardOutput "$d\abt\logs\dev2_sc$i.log" -RedirectStandardError "$d\abt\logs\dev2_sc$i.err" }
$p | Wait-Process
Get-Content "$d\abt\logs\dev2_sc0.log", "$d\abt\logs\dev2_sc1.log", "$d\abt\logs\dev2_sc0.err", "$d\abt\logs\dev2_sc1.err"
& $py dev_base2.py score
