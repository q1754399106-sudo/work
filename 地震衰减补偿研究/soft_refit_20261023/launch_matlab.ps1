$ErrorActionPreference = 'Stop'
$trialRoot = $PSScriptRoot
$matlabExe = 'D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe'
$batchExpression = "cd('$($trialRoot.Replace('\','/'))');run_soft_trial"
$matlabArgs = @('-wait','-nojvm','-batch',('"' + $batchExpression + '"'))
$trialProcess = Start-Process -FilePath $matlabExe -ArgumentList $matlabArgs -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $trialRoot 'matlab_stdout.log') -RedirectStandardError (Join-Path $trialRoot 'matlab_stderr.log')
$trialProcess.Id | Set-Content -LiteralPath (Join-Path $trialRoot 'matlab_pid.txt') -Encoding utf8
$trialProcess.WaitForExit()
$trialProcess.Refresh()
$trialExit = $trialProcess.ExitCode
$trialExit | Set-Content -LiteralPath (Join-Path $trialRoot 'matlab_exit_code.txt') -Encoding utf8
Write-Output "soft MATLAB exit: $trialExit"
exit $trialExit
