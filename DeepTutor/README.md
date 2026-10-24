# DeepTutor 本机安装

安装日期：2026-10-24。项目：https://github.com/HKUDS/DeepTutor 。

安装目录：`D:\ChatGPT\研究\DeepTutor`。

## 打开与启动

- 网页：<http://127.0.0.1:3782/>。
- 双击 `Start-DeepTutor.cmd`。已运行时打开网页；未运行时启动官方 Web 服务，并在就绪后打开浏览器。保留启动控制台，按 Ctrl+C 停止。
- `Open-DeepTutor.url` 只打开网页，不启动服务。
- `Stop-DeepTutor.cmd` 仅适用于用官方 `--detach` 方式启动的实例；控制台启动请用 Ctrl+C。

本次安装验证后服务保持运行。首次验证使用 Windows 隐藏后台进程，日志在 `data/user/runtime/server.stdout.log`、`server.stderr.log`。今后使用双击入口启动时会显示控制台。

## 模型配置

界面已设为中文。打开设置中的“模型与服务”，配置自己的提供商、模型与认证。知识库检索还需要嵌入模型，联网搜索需要相应搜索服务。当前未配置模型，也未进行付费在线模型请求。

## 独立环境与复原

- DeepTutor：1.6.15（官方 PyPI 完整 Web + CLI 包）。
- Python：3.12.15，位于 `.tools/python/`；不替换系统 Python 3.10。
- 虚拟环境：`.venv/`；使用本机 Node.js 24.20.0。
- `requirements-lock.txt` 固定本次安装的 175 个依赖版本。
- uv：0.13.5；下载自 PyPI，并核对其公布的 SHA-256。
- 在 Windows PowerShell 中运行 `powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\Install-DeepTutor.ps1` 可以复原环境；引导阶段需要已有的 Python Launcher（`py -3`）及网络。

启动入口将后端等待上限设为 300 秒、前端设为 180 秒。本机首次启动超过官方默认的 60 秒；网页缓存准备完毕后已验证成功。

## 已完成验证

- 官方 CLI 可运行。
- `uv pip check`：175 个包依赖兼容；复原安装脚本重复运行成功。
- `http://127.0.0.1:8001/health/live`：HTTP 200。
- `http://127.0.0.1:8001/health/ready`：HTTP 200。
- `http://127.0.0.1:3782/`：HTTP 200。
- 浏览器实际显示 DeepTutor 中文聊天界面；模型选择和发送按钮因未配置模型而不可用。
- `deeptutor doctor`：运行存储通过；模型和凭证待配置，在线调用未测试。

## 数据与归档

`data/` 保存设置、模型凭证、学习资料和会话；`.tools/`、`.venv/`、`.cache/` 保存本地运行环境和缓存。它们全部被本目录的 `.gitignore` 排除。仅归档本说明、启动入口、复原脚本和依赖版本。运行中的程序与本地数据保留，不执行归档后清理。
