# QScan＋SC 在一维数据 F3/F4、V4/V5 上的正交对照

**状态：已完成。** 先读 [RESULTS_ZH.md](RESULTS_ZH.md)：Kennett 组 QScan＋SC 稳定更好（24/24），V 组反而是 QScan＋标准版稳定更好（24/24），结果按生成器分化。

- `PLAN_ZH.md`：预登记方案；`PRE_RUN_HASHES.txt`：运行前的代码、方案与输入哈希。
- `RESULTS_ZH.md`：结果与判定（运行后写入）。
- `scripts/`：
  - `snapshot_inputs.py`：复制 novelty 的 SC 代码与参数（提交 ff95a69，逐字节核对与 HEAD 一致）、round3 的 F3/F4/V4/V5 道集、真值、归档结果与日志，记录 novelty 本地 Q+0 SC 输出的 SHA-256 → `cache/SNAPSHOT.json`；
  - `run_sc.py`：按所用 Q（SC）和 round3 归档的 QScan 选择（QScan＋SC）原样运行 SC → `cache/solved/`；先算 Q+0 的 SC 并与 novelty 输出核对；
  - `score.py` → `analysis/rows.csv`、`analysis/decisions.json`。

运行顺序（Python 3.12.14、numpy 2.3.5、scipy 1.18.1，BLAS 单线程）：

```
python scripts/snapshot_inputs.py
python scripts/run_sc.py F3 V4     # 三个进程：F3 V4 / F4 / V5
python scripts/score.py
```

需要同级目录 `../novelty_evidence_20261025/cache/`（scipy 依赖与 Q+0 SC 输出，本地，不上传）。
