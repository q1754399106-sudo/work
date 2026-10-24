$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$env:UV_PYTHON_INSTALL_DIR = Join-Path $PSScriptRoot '.tools\python'
$env:UV_PYTHON_BIN_DIR = Join-Path $PSScriptRoot '.tools\bin'
$env:UV_CACHE_DIR = Join-Path $PSScriptRoot '.cache\uv'
& py -3 (Join-Path $PSScriptRoot 'bootstrap_uv.py')
if ($LASTEXITCODE -ne 0) { throw 'uv bootstrap failed' }
$uvexe = (Get-ChildItem -LiteralPath '.tools\uv' -Filter uv.exe -Recurse | Select-Object -First 1).FullName
& $uvexe python install 3.12.15
if ($LASTEXITCODE -ne 0) { throw 'Python installation failed' }
if (-not (Test-Path -LiteralPath '.venv\Scripts\python.exe')) {
    & $uvexe venv --python 3.12.15 .venv
    if ($LASTEXITCODE -ne 0) { throw 'Virtual environment creation failed' }
}
& $uvexe pip install --python '.venv\Scripts\python.exe' -r requirements-lock.txt
if ($LASTEXITCODE -ne 0) { throw 'DeepTutor installation failed' }
& $uvexe pip check --python '.venv\Scripts\python.exe'
if ($LASTEXITCODE -ne 0) { throw 'Dependency check failed' }
Write-Host 'Ready. Double-click Start-DeepTutor.cmd to launch.'
