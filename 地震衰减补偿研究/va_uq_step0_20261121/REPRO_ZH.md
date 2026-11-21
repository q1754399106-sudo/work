# 复现记录（2026-11-21）

- 日期：2026-11-21，Windows 11，本机（16 逻辑核，同时有别的 MATLAB/Python 长任务在跑）。
- 解释器：codex 自带 Python 3.12.14（MSC v.1944 64 bit），NumPy 2.3.5；均以 `-I` 运行。README 写的是 Python 3.11。
- SciPy：codex Python 里没有 SciPy，没有安装。只有 `t_r0band.py` 要用（`scipy.optimize.brentq`），用一个包装脚本在 `-I` 下把仓库里已有的 `novelty_evidence_20261025/cache/deps`（SciPy 1.18.1）追加到 `sys.path`，再 `runpy.run_path` 运行。其余三个脚本只用 NumPy。
- 运行前读过全部 6 个脚本：只做 NumPy/SciPy 计算，不联网，只写命令行给定的输出文件（`t_convention_noise.py` 按路径导入 `t_convention.py`）。
- README 的 4 条命令并行运行，输出写到会话临时目录，没有覆盖本目录 `analysis/`。
- 比较方法：`.log` 忽略行尾 CR 后逐字节比较；`.json` 递归比较全部数值（NaN 须同为 NaN），键和长度必须一致。

| 脚本 | 用时（4 个并行） | 日志 | JSON |
|---|---|---|---|
| `t_r0band.py` | 35 s | 一致（逐字） | 无 JSON 输出 |
| `t_hfocus3.py` | 27 min 28 s | 一致（逐字） | 106 个数值，最大绝对差 1.3e-14，结构一致 |
| `t_convention.py` | 16 min 9 s | 一致（逐字） | 540 个数值，最大绝对差 5.0e-14，结构一致 |
| `t_convention_noise.py` | 15 min 41 s | 一致（逐字） | 120 个数值，最大绝对差 3.3e-16，结构一致 |

结论：4 条命令全部复现，差异都在浮点舍入量级（≤ 5e-14，远小于 1e-6 的容差）。用时比 README 写的“2 核不到 10 分钟”长，原因是与别的任务共用机器。

没有重跑：`t_hfocus.py`、`t_hfocus2.py`（README 的运行命令里没有列出；`t_hfocus2` 的单频结论已由 `t_hfocus3` 取代）。
