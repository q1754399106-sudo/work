[CmdletBinding()]
param([switch]$Apply)

# Fixed allowlist. No recursive deletion, registry changes, or system-folder cleanup.
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$stateRoot = 'D:\ChatGPT\CDriveMaintenance\state'
$quarantine = Join-Path $stateRoot 'quarantine'
$logs = Join-Path $stateRoot 'logs'
$now = [DateTime]::UtcNow
$cutoff = $now.AddDays(-7)
$roots = @(
    (Join-Path $env:LOCALAPPDATA 'Temp'),
    (Join-Path $env:LOCALAPPDATA 'CrashDumps')
)
$extensions = @('.tmp', '.temp', '.log', '.etl', '.dmp', '.cache')
$runLimit = 5GB
$storeLimit = 10GB
$timer = [Diagnostics.Stopwatch]::StartNew()

function Assert-PlainPath([string]$Path, [string]$Root) {
    $full = [IO.Path]::GetFullPath($Path)
    $base = [IO.Path]::GetFullPath($Root).TrimEnd('\')
    if ($full -ne $base -and -not $full.StartsWith($base + '\', [StringComparison]::OrdinalIgnoreCase)) {
        throw "Path outside allowed root: $full"
    }
    # Inspect every existing ancestor, including ancestors above Root.
    $current = $full
    while ($current) {
        if (Test-Path -LiteralPath $current) {
            $item = Get-Item -LiteralPath $current -Force
            if ($item.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw "Reparse point: $current" }
        }
        $parent = [IO.Directory]::GetParent($current)
        if ($null -eq $parent) { break }
        $current = $parent.FullName
    }
    return $full
}

function Get-FreeBytes([string]$Drive) {
    return ([IO.DriveInfo]::new($Drive)).AvailableFreeSpace
}

function Save-Json($Value, [string]$Path) {
    $Value | ConvertTo-Json -Depth 7 | Set-Content -LiteralPath $Path -Encoding UTF8
}

Assert-PlainPath $stateRoot 'D:\ChatGPT\CDriveMaintenance' | Out-Null
New-Item -ItemType Directory -Path $quarantine, $logs -Force | Out-Null
$lock = $null
try {
    $lock = [IO.File]::Open((Join-Path $stateRoot 'run.lock'), 'OpenOrCreate', 'ReadWrite', 'None')
} catch { throw 'Another cleanup is running, or the state directory is unavailable.' }

$report = [ordered]@{
    startedUtc = $now.ToString('o'); mode = $(if ($Apply) { 'apply' } else { 'preview' })
    powerShellVersion = $PSVersionTable.PSVersion.ToString()
    cFreeBefore = (Get-FreeBytes 'C:\'); cFreeAfter = 0
    eligibleCount = 0; eligibleBytes = [long]0; movedCount = 0; movedBytes = [long]0
    purgedCount = 0; purgedBytes = [long]0; skipped = 0
    warnings = [Collections.Generic.List[string]]::new()
    sample = [Collections.Generic.List[object]]::new()
}
try {
    Save-Json @{startedUtc=$now.ToString('o'); processId=$PID; powerShellVersion=$PSVersionTable.PSVersion.ToString(); mode=$report.mode} (Join-Path $stateRoot 'started.json')
    # Only purge files created by this script, after 30 days and an unchanged SHA-256.
    foreach ($meta in @(Get-ChildItem -LiteralPath $quarantine -File -Filter '*.json' -Force)) {
        try {
            Assert-PlainPath $meta.FullName $quarantine | Out-Null
            $entry = Get-Content -LiteralPath $meta.FullName -Raw | ConvertFrom-Json
            if ($entry.status -ne 'quarantined' -or [DateTime]::Parse($entry.quarantinedUtc).ToUniversalTime() -ge $now.AddDays(-30)) { continue }
            $expected = Join-Path $quarantine ($meta.BaseName + '.data')
            if ($entry.storedPath -ne $expected) { throw 'Invalid quarantine record path.' }
            Assert-PlainPath $expected $quarantine | Out-Null
            $stored = Get-Item -LiteralPath $expected -Force
            if ($stored.PSIsContainer -or $stored.Length -ne $entry.length) { throw 'Quarantined file size changed.' }
            if ((Get-FileHash -LiteralPath $expected -Algorithm SHA256).Hash -ne $entry.sha256) { throw 'Quarantined file hash changed.' }
            if ($Apply) {
                Remove-Item -LiteralPath $expected -Force
                $entry.status = 'purged'
                Save-Json $entry $meta.FullName
                $report.purgedCount++
                $report.purgedBytes += $stored.Length
            }
        } catch { $report.warnings.Add("Quarantine skipped: $($_.Exception.Message)") }
    }
    $usage = [long]0
    foreach ($data in @(Get-ChildItem -LiteralPath $quarantine -File -Filter '*.data' -Force)) { $usage += $data.Length }
    foreach ($root in $roots) {
        if (-not (Test-Path -LiteralPath $root -PathType Container)) { continue }
        Assert-PlainPath $root $root | Out-Null
        if (-not ([IO.Path]::GetFullPath($root)).StartsWith('C:\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Expected a C: source root.' }
        $pending = [Collections.Generic.Stack[string]]::new()
        $pending.Push($root)
        while ($pending.Count -gt 0 -and $timer.Elapsed.TotalSeconds -lt 180) {
            $directory = $pending.Pop()
            try {
                Assert-PlainPath $directory $root | Out-Null
                $children = @(Get-ChildItem -LiteralPath $directory -Force -ErrorAction Stop)
            } catch { $report.skipped++; continue }
            foreach ($file in $children) {
                if ($timer.Elapsed.TotalSeconds -ge 180) { break }
                if ($file.Attributes -band ([IO.FileAttributes]::ReparsePoint -bor [IO.FileAttributes]::System -bor [IO.FileAttributes]::ReadOnly)) { continue }
                if ($file.PSIsContainer) { $pending.Push($file.FullName); continue }
                if ($file.LastWriteTimeUtc -ge $cutoff -or $file.CreationTimeUtc -ge $cutoff -or $file.Length -eq 0) { continue }
                if ($root.EndsWith('CrashDumps')) { if ($file.Extension.ToLowerInvariant() -ne '.dmp') { continue } }
                elseif ($extensions -notcontains $file.Extension.ToLowerInvariant()) { continue }
                $report.eligibleCount++
                $report.eligibleBytes += $file.Length
                if ($report.sample.Count -lt 20) { $report.sample.Add(@{path=$file.FullName; bytes=$file.Length}) }
                if (-not $Apply) { continue }
                if ($report.movedCount -ge 500 -or $report.movedBytes + $file.Length -gt $runLimit -or $usage + $file.Length -gt $storeLimit -or (Get-FreeBytes 'D:\') - $file.Length -lt 20GB) { $report.skipped++; continue }
                $source = $null; $target = $null; $sha = $null; $destination = $null
                try {
                    Assert-PlainPath $file.FullName $root | Out-Null
                    # Exclusive read access skips files held by applications and blocks new writers during copying.
                    $source = [IO.File]::Open($file.FullName, 'Open', 'Read', 'None')
                    $fresh = Get-Item -LiteralPath $file.FullName -Force
                    if ($fresh.LastWriteTimeUtc -ge $cutoff -or $fresh.Length -ne $file.Length) { throw 'File changed during scan.' }
                    $sha = [Security.Cryptography.SHA256]::Create()
                    $hash = [BitConverter]::ToString($sha.ComputeHash($source)).Replace('-', '')
                    $source.Position = 0
                    $id = [Guid]::NewGuid().ToString('N')
                    $destination = Join-Path $quarantine ($id + '.data')
                    $metaPath = Join-Path $quarantine ($id + '.json')
                    Assert-PlainPath $destination $quarantine | Out-Null
                    $entry = [ordered]@{status='copying'; originalPath=$file.FullName; storedPath=$destination; length=$fresh.Length; sha256=$hash; lastWriteUtc=$fresh.LastWriteTimeUtc.ToString('o'); quarantinedUtc=$now.ToString('o')}
                    Save-Json $entry $metaPath
                    $target = [IO.File]::Open($destination, 'CreateNew', 'ReadWrite', 'None')
                    $source.CopyTo($target)
                    $target.Flush($true)
                    $target.Position = 0
                    $copiedHash = [BitConverter]::ToString($sha.ComputeHash($target)).Replace('-', '')
                    if ($copiedHash -ne $hash -or $target.Length -ne $source.Length) { throw 'Copy verification failed; source preserved.' }
                    $target.Dispose(); $target = $null
                    $entry.status = 'copied_verified'
                    Save-Json $entry $metaPath
                    $source.Dispose(); $source = $null
                    Assert-PlainPath $file.FullName $root | Out-Null
                    $final = Get-Item -LiteralPath $file.FullName -Force
                    if ($final.PSIsContainer -or $final.LastWriteTimeUtc -ne $fresh.LastWriteTimeUtc -or $final.Length -ne $fresh.Length) { throw 'Source changed; source preserved.' }
                    Remove-Item -LiteralPath $file.FullName -Force
                    $entry.status = 'quarantined'
                    Save-Json $entry $metaPath
                    $report.movedCount++; $report.movedBytes += $fresh.Length
                } catch {
                    $report.skipped++
                    if ($report.warnings.Count -lt 20) { $report.warnings.Add("File skipped: $($file.FullName): $($_.Exception.Message)") }
                } finally {
                    if ($source) { $source.Dispose() }
                    if ($target) { $target.Dispose() }
                    if ($sha) { $sha.Dispose() }
                    # Include incomplete copies in the storage budget as well.
                    if ($destination -and (Test-Path -LiteralPath $destination -PathType Leaf)) {
                        $usage += (Get-Item -LiteralPath $destination -Force).Length
                    }
                }
            }
        }
    }
    if ($timer.Elapsed.TotalSeconds -ge 180) { $report.warnings.Add('Time budget reached; remaining files deferred.') }
    if ($usage -ge $storeLimit) { $report.warnings.Add('Quarantine storage limit reached.') }
    $report.cFreeAfter = Get-FreeBytes 'C:\'
    $report['quarantineBytes'] = $usage
    $report['finishedUtc'] = [DateTime]::UtcNow.ToString('o')
    $logPath = Join-Path $logs ($now.ToString('yyyyMMddTHHmmssfffZ') + '-' + $report.mode + '.json')
    Save-Json $report $logPath
    Save-Json $report (Join-Path $stateRoot 'latest.json')
    $report | ConvertTo-Json -Depth 7
} finally { if ($lock) { $lock.Dispose() } }
