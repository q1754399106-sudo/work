# 为什么 V 组偏向相干标准版：生成器因子实验

**状态：方案已冻结，尚未运行。** 方案见 [PLAN_ZH.md](PLAN_ZH.md)，运行前哈希见 `PRE_RUN_HASHES.txt`。

上游：`qscan_sc_round3_20261026`。那一轮在 Kennett 的 F3/F4 上 QScan＋SC 稳定更好，在 V4/V5 上 QScan＋相干标准版稳定更好。本轮在 24 个全新模型上，用 2×2×2 配对设计拆开三项生成器差别：A 反射系数稀疏/密集，B 横向有无构造，C 正演是 Kennett 还是求解器自身的算子。

## 内容

- `scripts/snapshot_inputs.py`：把 round3 的 `core/`、`source/`、`run_combo.m` 逐字节复制到 `cache/runner/`，生成只改误差档的 `run_combo3.m`；复制 novelty 的 SC 代码（ff95a69）→ `cache/SNAPSHOT.json`。
- `scripts/make_models.m`：生成 Z1–Z3 × 8 个格子 → `models/`（求解端）、`truth/`（评分端）。
- `scripts/run_matlab.py`：`make_models`、`stage_dev`（记录模型哈希，复制到 runner）、`run_combo3`（标准版、QScan＋标准版等）。
- `scripts/run_sc.py`：SC、QScan＋SC → `cache/solved/`。
- `scripts/score.py` → `analysis/rows.csv`、`models.csv`、`decisions.json`。

## 运行

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_matlab.py make_models
& $py scripts\run_matlab.py stage_dev
& $py scripts\run_matlab.py run_combo3 Z1_000 Z1_001 ...   # 两个 MATLAB 进程分担 24 个模型
& $py scripts\run_sc.py Z1_000 Z1_001 ...                  # 最多三个进程
& $py scripts\score.py
```

需要 MATLAB R2024b 和同级目录 `../novelty_evidence_20261025/cache/deps`（scipy，本地）。`cache/` 不上传。
