# 实验 12：改进方案对"同一个 Q、别处保守"的基线

问题：实验 11 的改进方案（SR 估 Q，SC-a 振幅用 c_SR×1.3，再把相位校正到 c_SR）"振幅、相位用不同 Q"，它是否胜过振幅、相位用同一个 Q、靠别处保守的经典做法？比较对象：
- 调好稀疏权重的 SC-a；
- 调好增益的 Wang 稳定化反 Q；
- 纯相位反 Q＋时变 Wiener（van der Baan 型）。

方案在 `PLAN_ZH.md`，冻结时的哈希在 `PRE_RUN_HASHES.txt`。

- `scripts/snapshot_inputs.py`：把输入恢复到被忽略的 `cache/inputs/` 并校验，包括实验 7–11 的程序、实验 10/11 的本地 H 组输出（记录哈希）。
- `scripts/lib12.py`：新代码，包括稳定化反 Q `iq`（振幅、相位可用不同系数）、时变 Wiener `tvw`、综合分 S。
- `scripts/dev12.py`、`dev12_ext.py`、`dev12_ext2.py`、`dev12_ext3.py`、`devsum12.py`：开发阶段（实验 10 的 H 组，已见数据）。
  - 结果在 `analysis/dev/`、`analysis/dev_ext2/`、`analysis/dev_ext3/`；
  - 汇总在 `analysis/dev_summary.json`，选参在 `analysis/dev_selection.json`。
- `scripts/run12.py`：检验集 T（20 个，全新）。输出 `results/synth/`，输出道存在本地 `cache/out/`。
- `scripts/analyze12.py`：判定，输出 `analysis/summary.json`。
- `scripts/figures12.py`：出图，输出 `figures/`。
