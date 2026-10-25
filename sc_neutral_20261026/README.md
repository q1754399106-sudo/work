# 中立数据 H1–H6 上的相干标准版对 SC 文献基线

- `PLAN_ZH.md`：预登记方案（提交 18c85cb）；`PRE_RUN_HASHES.txt`：运行前的代码、方案与输入哈希。
- `RESULTS_ZH.md`：结果与判定。
- `scripts/`：
  - `snapshot_inputs.py`：把输入从登记的提交复制到被忽略的 `cache/inputs/`，逐字节核对与 HEAD 一致，并记录 fd2d 扫描缓存（本地，未上传）126 个文件的 SHA-256；
  - `prep_noise.m`（由会话“研究后续方向”编写并运行）：用 fd2d 的种子重建 H1–H6 的含噪道集，要求重估的 σ 与扫描缓存逐位相同（需要 4 线程）→ `cache/noisy/`、`logs/prep_noise.*`；
  - `dict_check.m`、`check_dictionary.py`：运行前核对 SC 用的 Python 字典与 Standard 用的冻结 MATLAB 字典 → `analysis/dictionary_check.json`；
  - `run_sc.py`（由会话“研究后续方向”编写，本会话运行）：逐字调用 novelty 的 `benchmark.solve`、`predictor` 和 `common.dictionary`，参数冻结为 SC（0.125, 10）和 SC0（0.5, 0）→ `cache/solved/`；
  - `score.py` → `analysis/replicates.csv`、`analysis/scenes.csv`、`analysis/decisions.json`。

运行顺序（Python 3.12.14、numpy 2.3.5、scipy 1.18.1，BLAS 单线程；MATLAB R2024B）：

```
python scripts/snapshot_inputs.py
matlab -batch "addpath('scripts');prep_noise"
matlab -batch "addpath('scripts');dict_check"
python scripts/check_dictionary.py
python scripts/run_sc.py H1 H2      # 三个进程：H1 H2 / H3 H4 / H5 H6
python scripts/score.py
```

`run_sc.py` 和 `score.py` 需要 scipy：放在 `novelty_evidence_20261025/cache/deps` 下（不上传）。
