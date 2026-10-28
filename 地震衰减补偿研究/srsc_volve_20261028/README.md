# 实验 11：改进 QScan＋SC——SR 相位校正＋保守振幅

问题：把 Q 系数改由 SR 估计，SC-a 的振幅补偿用 c_SR×1.3（偏保守），再把频散相位校正到 c_SR。这个改进方案（SRSCa_1.3+ph）是否优于现流程 QScan＋SC-a？方案在 `PLAN_ZH.md`，冻结时的哈希在 `PRE_RUN_HASHES.txt`。

- `scripts/snapshot_inputs.py`：把输入恢复到被忽略的 `cache/inputs/` 并校验，包括实验 7、8、9、10 的程序，以及实验 10 的 H 组输出（本地文件，记录哈希）。
- `scripts/sc11.py`：公共代码，包括 SC-a、QScan、相位校正 `phase_fix` 和 M1。
- `scripts/dev11.py`、`devsum11.py`：开发阶段（实验 10 的 H 组，已见数据）。结果在 `analysis/dev/`、`analysis/dev_summary.json`。
- `scripts/run11.py`：检验集 T/L 和实测 inline 10141。输出 `results/`，输出道存在本地 `cache/`。
- `scripts/analyze11.py`：判定，输出 `analysis/summary.json`。
