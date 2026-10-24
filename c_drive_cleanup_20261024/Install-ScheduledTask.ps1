$ErrorActionPreference = 'Stop'
$taskName = 'Codex-CDriveSafeCleanup'
$runtime = 'D:\ChatGPT\CDriveMaintenance'
if (-not (Test-Path -LiteralPath (Join-Path $runtime 'Cleanup-CDrive.ps1') -PathType Leaf)) {
    throw 'Deploy Cleanup-CDrive.ps1 to D:\ChatGPT\CDriveMaintenance first.'
}
if (Get-ScheduledTask -TaskName $taskName -ErrorAction SilentlyContinue) {
    throw 'Task already exists. Inspect it before replacing or changing it.'
}
if ((Get-TimeZone).Id -ne 'Taipei Standard Time') { throw 'Review schedule: expected Taipei time zone.' }
$powerShellExe = Join-Path $PSHOME 'pwsh.exe'
if (-not (Test-Path -LiteralPath $powerShellExe -PathType Leaf)) { throw 'Run this installer with PowerShell 7.' }
$action = New-ScheduledTaskAction -Execute $powerShellExe -Argument '-NoProfile -NonInteractive -WindowStyle Hidden -File "D:\ChatGPT\CDriveMaintenance\Cleanup-CDrive.ps1" -Apply' -WorkingDirectory $runtime
$trigger = New-ScheduledTaskTrigger -Daily -At '09:00'
$settings = New-ScheduledTaskSettingsSet -StartWhenAvailable -MultipleInstances IgnoreNew -ExecutionTimeLimit (New-TimeSpan -Minutes 10) -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries
$settings.Priority = 4
$principal = New-ScheduledTaskPrincipal -UserId ([Security.Principal.WindowsIdentity]::GetCurrent().Name) -LogonType Interactive -RunLevel Limited
Register-ScheduledTask -TaskName $taskName -Action $action -Trigger $trigger -Settings $settings -Principal $principal -Description 'Daily safe C-drive temp cleanup with verified D-drive quarantine retained for 30 days.'
