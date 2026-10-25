# 二维黏声有限差分正演下的复验

**状态：已完成。** 结果见 [RESULTS_ZH.md](RESULTS_ZH.md)：E1 判为“偏移不稳定”，E2–E4 成立。冻结方案见 [PLAN_ZH.md](PLAN_ZH.md)，哈希见 `logs/PRE_RUN_HASHES.json`、`logs/PRE_SWEEP_HASHES.json`。

上游：`qscan_offset_20261025`（QScan 偏移随模型变化；无噪声最优在上尺度 Q0）。这两条结论只在 Kennett 一维细层正演上得到。本轮换成二维黏声有限差分（GSLS 衰减、断层绕射）检验它们能否复现。不使用 U1–U4。

## 内容

- `scripts/fd2d_visco.m`：二维黏声有限差分核心（8 阶交错网格、4 机制 GSLS、C-PML、单精度）。
- `scripts/gsls_design.m`：GSLS 弛豫时间与每个 Q 的 τ（5–100 Hz 极小极大拟合）。
- `scripts/fd_setup.m`：路径、常数、H1–H6 的种子与宏观模型、断层参数、场景、噪声种子规则。
- `scripts/fd_geology.m`：细层生成器（与 `qscan_offset_20261025/scripts/make_dev.m` 逐字相同）铺成二维网格，加断层。
- `scripts/fd_shot.m`、`scripts/process_zo.m`：单炮正演，以及零偏移距处理（去直达波、去子波与线源滤波、二维扩散校正、重采样）。
- `scripts/test_fd2d.m`、`scripts/diag_fd_accuracy.m`、`scripts/bench_fd.m`：正演验证、GSLS 与网格收敛诊断、测速。只用测试地层。
- `scripts/make_shots.m`、`scripts/assemble_models.m`：逐炮生成（可续跑），汇成公开配置、无噪声道集（`models/`）与评分目标（`truth/`）。
- `scripts/run_fd.m`：扫描 worker，与 `qscan_offset_20261025/scripts/run_offset.m` 逐句相同。
- `scripts/summarize_fd.m`：误差曲线、D/S1/S2/B/A/最优、偏移、无噪声最优，以及 E1–E4 判定。
- `scripts/summarize_fault.m`、`scripts/fault_compare.py`：断层附近误差的补算（冻结后加入，只作描述），见 [FAULT_ZH.md](FAULT_ZH.md)。
- `scripts/qscan_select.m`、`noise_replica.m`、`estimate_sigmas.m`：与 `qscan_offset_20261025/scripts` 中的文件逐字相同。
- `scripts/snapshot_inputs.py`、`hash_prerun.py`、`run_stage.py`：输入快照、两个冻结点的哈希、带日志与超时的 MATLAB 运行器。
- `cache/`（不上传）：输入快照、每炮原始记录（约 40 KB/炮）、每个实现的扫描输出。

## 运行

需要 MATLAB R2024b 和 Python。每个 MATLAB `-nojvm` 进程约需 2 GB 内存，最多同时跑 2 个。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_stage.py test_fd2d
& $py scripts\hash_prerun.py plan
& $py scripts\run_stage.py make_shots 1 0 99     # 模型 1、炮 0–99；模型 1–6
& $py scripts\run_stage.py assemble_models
& $py scripts\hash_prerun.py models
& $py scripts\run_stage.py run_fd 1 1 0 6        # 模型 1、场景 1、实现 0–6；模型 1–6 × 场景 1–3
& $py scripts\run_stage.py summarize_fd
```
