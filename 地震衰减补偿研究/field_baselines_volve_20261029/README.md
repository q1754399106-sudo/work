# 实验 13：实验 12 的方法用到 Volve 实测（多井标定）

问题：实验 12 在半合成数据上的结论，在实测井上方向是否还成立？比较的方法：
- 调好增益的稳定化反 Q（70 dB）；
- 稀疏项极小的 SC-a；
- 实验 11 的改进方案；
- 用户的流程 QScan＋SC-a，原样一版、调小 a 一版。

参数全部沿用实验 12 冻结的值，不在实测上调。方案在 `PLAN_ZH.md`，冻结时的哈希在 `PRE_RUN_HASHES.txt`。

- `scripts/snapshot_inputs.py`：把输入恢复到被忽略的 `cache/inputs/` 并校验，包括实验 7–12 的程序、实验 9 的井资料和 MIG_RAW 哈希。
- `scripts/lib13.py`：公共代码，包括整条 inline 读取、SR/QScan 系数按 inline 取值、方法集。
- `scripts/run13.py`：
  - `lines` 计算要处理的 inline，写 `analysis/lines.json`；
  - `qscan` 在 10 条线上算 QScan 系数，写 `results/qscan/`；
  - `a`、`b` 处理各条 inline，写 `results/lines/`，输出道存在本地 `cache/out/`。
- `scripts/tie13.py`：井震相关，写 `results/tie.json`。
- `scripts/analyze13.py`：判定，写 `analysis/summary.json`。
- `scripts/figures13.py`：出图（冻结后写，只读结果），写 `figures/`。

**结果**：见 `RESULTS_ZH.md`。四项判定三项支持（F1、F2、F3）、F4 无法判定，没有相反；实测方向与实验 12 半合成一致。
