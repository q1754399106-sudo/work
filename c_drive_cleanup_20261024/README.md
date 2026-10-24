# C 盘自动维护

由 Windows 计划任务 `Codex-CDriveSafeCleanup` 每天台北时间 09:00 执行，电脑开机且当前用户登录时运行，错过计划时间后在条件满足时补跑，不依赖 Codex 保持打开。部署副本位于 `D:\ChatGPT\CDriveMaintenance\Cleanup-CDrive.ps1`，该副本是持久运行依赖，不属于研究仓库工作副本清理范围。

任务使用本机已验证的 PowerShell 7，以当前用户普通权限隐藏运行，防止同一任务重叠，最长运行 10 分钟。脚本自带独占锁。`Install-ScheduledTask.ps1` 在 PowerShell 7 中用于首次安装任务；已有同名任务时停止，避免覆盖。状态和上次运行结果可在 Windows 任务计划程序中查看，执行详情记录在下面的本地日志中。暂停可执行 `Disable-ScheduledTask -TaskName Codex-CDriveSafeCleanup`，恢复使用 `Enable-ScheduledTask -TaskName Codex-CDriveSafeCleanup`。

任务使用普通进程优先级，使其在研究计算持续运行时仍能获得执行时间。

只处理当前用户 `AppData\Local\Temp` 内的 `.tmp/.temp/.log/.etl/.dmp/.cache`，以及 `AppData\Local\CrashDumps` 内的 `.dmp`。创建时间和最后修改时间均须超过 7 天。跳过占用、只读、系统文件和所有重解析点。不会递归删除目录，不处理研究资料、下载、回收站、浏览器登录数据和 Windows 系统目录。

先独占打开文件并复制到 D 盘，验证 SHA-256 和大小，再重新确认源文件未变化后移除原文件。隔离文件保留 30 天，仅清除有成功隔离记录且哈希未变化的文件。每轮最多 5 GiB、500 个文件；隔离区上限 10 GiB，D 盘至少保留 20 GiB。目录遍历预算 180 秒，单文件复制可能超出此时间。复制失败时保留源文件；中断或未完成的隔离记录不会自动清除。

预览（默认，不移动或删除文件，但写入日志）：

```powershell
& 'D:\ChatGPT\CDriveMaintenance\Cleanup-CDrive.ps1'
```

执行：

```powershell
& 'D:\ChatGPT\CDriveMaintenance\Cleanup-CDrive.ps1' -Apply
```

日志：`D:\ChatGPT\CDriveMaintenance\state\logs`；最近结果：`state\latest.json`；隔离区：`state\quarantine`。每个 `.json` 记录原路径、隔离数据路径和 SHA-256。恢复时暂停自动任务，核对哈希，确认原位置没有新文件后，将对应 `.data` 复制回记录中的 `originalPath`；不要覆盖新文件。

Windows 自带存储感知已开启，保留其当前配置。此脚本补充定向清理，不能保证一次解决低空间。C 盘迁移旧副本仍需单独核查，不属于自动清理范围。

参考项目：[WinDirStat](https://github.com/windirstat/windirstat)、[BleachBit](https://github.com/bleachbit/bleachbit)。本脚本没有下载或执行这些项目的代码，也不依赖第三方清理器。
