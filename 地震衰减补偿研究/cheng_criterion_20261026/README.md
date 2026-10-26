# Cheng 2025 / Ma 2018 选 Q 准则的偏移检验

**状态：方案已冻结，尚未运行。** 方案见 [PLAN_ZH.md](PLAN_ZH.md)，运行前哈希见 `PRE_RUN_HASHES.txt`。

把 Cheng et al. 2025 的浅深一致性准则和 Ma et al. 2018 型 lp 准则接到与 QScan 相同的 lasso 扫描上，在 G1–G8、H1–H6 上检验选 Q 偏移是否因模型而异。框架沿用 `yan_criterion_20261026`。

## 内容

- `scripts/snapshot_inputs.py`：与 `yan_criterion_20261026/scripts/snapshot_inputs.py` 逐字相同。把 Git 输入复制到 `cache/inputs`，并记录本地扫描缓存的 SHA-256。
- `scripts/cheng_stats.m`：逐道集重建含噪数据（σ 须逐位一致），对 13 个因子做冻结的 lasso（核对 QScan 选择），保存补偿数据 C = W·(X/scale) 的浅深窗平均振幅谱和 lp 和。
- `scripts/score.py`：各准则的选择、偏移、C1 判定和描述量。
- `cache/`（不上传）：输入快照与逐道集统计量。

## 运行

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py -I scripts\snapshot_inputs.py
matlab -batch "addpath('scripts');cheng_stats('G','1','2','3','4','5','6','7','8')"
matlab -batch "addpath('scripts');cheng_stats('H','1','2','3','4','5','6')"
& $py -I scripts\score.py
```
