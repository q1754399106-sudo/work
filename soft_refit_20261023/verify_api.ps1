$ErrorActionPreference = 'Stop'
$trialRoot = $PSScriptRoot
$verifyExpression = "cd('$($trialRoot.Replace('\','/'))');verify_public_entry"
$verifyArgs = @('-wait','-nojvm','-batch',('"' + $verifyExpression + '"'))
$verifyProcess = Start-Process -FilePath 'D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe' -ArgumentList $verifyArgs -WindowStyle Hidden -PassThru -Wait -RedirectStandardOutput (Join-Path $trialRoot 'api_stdout.log') -RedirectStandardError (Join-Path $trialRoot 'api_stderr.log')
$verifyProcess.ExitCode | Set-Content -LiteralPath (Join-Path $trialRoot 'api_exit_code.txt') -Encoding utf8
exit $verifyProcess.ExitCode
