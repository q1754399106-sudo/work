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
