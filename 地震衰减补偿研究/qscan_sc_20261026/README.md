# QScan＋SC 组合在 H1–H6 上的检验

- `PLAN_ZH.md`：预登记方案（提交 8fdea24）；`PRE_RUN_HASHES.txt`：运行前的代码、方案与输入哈希。
- `RESULTS_ZH.md`：结果与判定。
- `scripts/`：
  - `snapshot_inputs.py`：从登记的提交复制 Git 输入（fd2d 模型与真值、novelty 的 SC 代码与参数），逐字节核对与 HEAD 一致；并记录本地输入（fd2d 扫描缓存、sc_neutral 的重建道集与 SC 输出，共 462 个文件）的 SHA-256；
  - `run_qsc.py`：以 fd2d 归档的 QScan 选择缩放 Q，再原样运行 SC → `cache/solved/*_QSC.npz/json`；
  - `score.py` → `analysis/replicates.csv`、`analysis/scenes.csv`、`analysis/decisions.json`。

运行顺序（Python 3.12.14、numpy 2.3.5、scipy 1.18.1，BLAS 单线程）：

```
python scripts/snapshot_inputs.py
python scripts/run_qsc.py H1 H2     # 三个进程：H1 H2 / H3 H4 / H5 H6
python scripts/score.py
```

需要同级目录 `../novelty_evidence_20261025/cache/deps`（scipy，不上传）、`../sc_neutral_20261026/cache/` 与 `../fd2d_20261025/cache/sweep/`（本地，不上传）。
