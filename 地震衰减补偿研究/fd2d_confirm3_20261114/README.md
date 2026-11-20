# 二维黏声有限差分数据第三批（UQW-S 确认用，2026-11-14）

`fd2d_confirm_20261109` 的正演与处理脚本原样复制，只换了模型表（`scripts/fd_setup.m`）：6 个新模型 H13–H18，Q 形态与 H7–H12 一一对应（递增、递减、反转、中等、整体低、高且递减），参数见 `uqw_snr_20261114/PLAN_CONFIRM_S_ZH.md`。

- **用 Octave 运行**：
  - 本环境没有 MATLAB，`scripts/run_octave.py` 用 Octave 加 `octave_compat_v2_20261025` 运行；
  - `scripts/octave_shim/discretize.m` 补上 Octave 缺少的 `discretize`；
  - 快照脚本 `scripts/snapshot_inputs.py` 只加了目录前缀。
- **等价性核对**：`scripts/check/check_octave.m` 用归档的 H7 做核对。地质的上尺度 Q 相差 2.8×10⁻¹⁴；处理后的第 50 道与归档 MATLAB 道相对差 5.1×10⁻⁶。
- **出炮**：`scripts/run_all_shots.sh` 起 4 个进程，每炮约 135 秒。
- `models/H*.mat`、`truth/H*_truth.mat`：汇总后的无噪声道集与弹性目标。
- `cache/`（不上传）：输入快照与每炮原始记录。

汇总结果（`logs/assemble_models.log`）：上尺度 Q 分别为
- H13：45 53.7 89.6 122.4
- H14：110 61.2 56.5 46.1
- H15：55 71.3 54.6 107.1
- H16：75 55.7 79.0 104.2
- H17：36 43.9 60.7 78.6
- H18：140 84.0 90.8 83.4

衰减/弹性能量比 0.364、0.687、0.373、0.625、0.283、0.744。600 炮 4 个 Octave 进程约 4.5 小时，没有失败。

使用记录：2026-11-14 UQW-S 确认检验（`uqw_snr_20261114/RESULTS_CONFIRM_S_ZH.md`）。从此不再是未见数据。

注：`models/`、`truth/` 的 12 个 .mat（约 4 MB）暂按普通 Git 文件保存，没有用 LFS。原因是本会话的网络策略拦截了 lfs.github.com，LFS 上传无法校验。
