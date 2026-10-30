# 实验 15：断层与复杂 Q 下的确认检验

- 方案：[PLAN_ZH.md](PLAN_ZH.md)（开跑前提交 b79fc73）
- 结果：[RESULTS_ZH.md](RESULTS_ZH.md)
- 数据：[semisynth_fault_q_20261030](../semisynth_fault_q_20261030/DATA_ZH.md) 的构造 B、C（Z 组 12＋孪生 Y 组 12）

判定（Z 组）：
- H1 用户流程调小 a：支持，+0.123，12/12；
- H2 SR 型补偿胜过原流程：支持，70 dB 反 Q +0.237，SC-a +0.128；
- H3 补偿有用：支持；
- H4 锚定理论压缩斜率外推到非比例 Q：支持，A = 0.36，理论 0.33；
- H5 SR 比 QScan 更准：无法判定。

另外：本数据上 70 dB 反 Q 比 SC-a 高 0.11，用真值 Q 时 SC-a 也只有 0.77；Q 复杂度本身代价 ≤ 0.035；气藏下方 IQ70 比孪生 Y 低 0.15。

运行：`scripts/run15.py test K N` → `scripts/analyze15.py` → `scripts/figures15.py`。Python 与依赖路径见 `scripts/lib15.py`。
