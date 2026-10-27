# 实验 7：SC / QScan＋SC 的实测适配与应用

- `PLAN_ZH.md`：预登记方案（2026-10-28）；`PRE_RUN_HASHES.txt`：运行前的代码与方案哈希。
- `scripts/lib.py`：从数据估计子波、按相邻道相干信噪比限带输出、新字典上的 QScan、评价函数；求解器取自 `novelty_evidence_20261025`（ff95a69），未改。
- `scripts/snapshot_inputs.py`：从登记提交恢复实验 5/6 的幕帘、参照和对照输出到被忽略的 `cache/inputs/`。
- `scripts/test_lib.py`：合成道集自检（不用实测数据）。
- `scripts/run.py wells|inline`：两口井的运行、评价和门槛；inline 10141 的应用。
