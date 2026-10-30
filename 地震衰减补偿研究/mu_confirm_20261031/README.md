# 实验 16：SC-a 去掉横向相干项的确认检验

- 方案：[PLAN_ZH.md](PLAN_ZH.md)（开跑前提交）
- 数据：[第二批断层＋复杂 Q 半合成数据](../semisynth_fault_q2_20261031/DATA_ZH.md)，构造 D、E、F 共 12 个，封存

检验 SC-a 的 mu_relative 从 10 改为 0 是否更好，覆盖三条路线：SR 估 Q、用户的 QScan 选 Q（a 调小）、真值 Q。同时报告"QScan＋调小 a＋mu 0"相对用户原流程的总收益，以及 mu 0 的 SC-a 与 70 dB 反 Q 的比较。

运行：`scripts/hash_prerun.py` → `scripts/run16.py test K 4` → `scripts/analyze16.py`。
