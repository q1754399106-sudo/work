# 实验 8：Volve 半合成检验（真实噪声＋测井反射系数，Q 已知）

**结果见 `RESULTS_ZH.md`**；运行后的改动见 `PLAN_ADDENDUM_ZH.md`。

- `PLAN_ZH.md`：预登记方案（2026-10-28）；`PRE_RUN_HASHES.txt`：运行前的代码、方案和构建结果哈希。
- `scripts/semi.py`：读测线、估构造时移、提取实测噪声、1 ms 正演、常规稳定化反 Q 滤波、评价指标；实验 7 的 `lib.py`（b1d76f8）和求解器（ff95a69）原样导入。
- `scripts/snapshot_inputs.py`：从登记提交恢复输入到被忽略的 `cache/inputs/`，并核对 MIG_RAW 数据体。
- `scripts/test_semi.py`：不用实测数据的自检；`scripts/smoke.py`：登记前的冒烟测试（不看分数）。
- `scripts/run.py build|solve K N`：构建半合成数据与真实性检查（`analysis/build.json`）；并行运行各方法。
- `scripts/analyze.py`：H1–H5 判定与汇总表。
- `scripts/posthoc_lag.py`：事后描述，允许整体时移的 M1；`scripts/figures.py`：作图。
