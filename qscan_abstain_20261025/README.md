# “有把握才修正”的 QScan 规则：新模型确认

**状态：已完成，不确认（A1 通过，A2 未通过）。** 先读 [RESULTS_ZH.md](RESULTS_ZH.md)。方案见 [PLAN_ZH.md](PLAN_ZH.md)，冻结哈希见 `logs/PRE_RUN_HASHES.json`，封存见 `decisions/SEALED_SHA256.json`。

规则来自 `qscan_predictor_20261025` 的事后观察（用户选择的选项 (a)）：若 QScan 所选因子处的 lasso 支撑数 Z1 > 72.657，就把选择下移 2 个网格步长，否则保持 QScan 的选择。

## 内容

- `scripts/as_setup.m`：路径、K1–K8 的种子与宏观模型、4 个场景、冻结的规则。
- `scripts/make_models.m`：生成器，与 `qscan_offset_20261025/scripts/make_dev.m` 逐句相同。`make_models('G')` 重新生成 G1–G8 作检验，`make_models('K')` 生成 K1–K8。
- `scripts/run_confirm.m`：每个道集的决定（不读真值），以及标准版在 f=1、b、b−2 三处的输出。
- `scripts/seal_decisions.py`、`score.m`、`analyze.py`、`make_figure.py`。
- `models/`、`truth/`：K1–K8 的公开配置与无噪声道集，以及评分目标。
- `decisions/`：封存的逐道集决定。
- `cache/`（不上传）：输入快照与输出道集。

## 运行（Octave 8.4，需要 `octave_compat_v2_20261025`、NumPy、pandas、Matplotlib）

```bash
python3 scripts/snapshot_inputs.py
octave --eval "addpath('scripts'); make_models('G')"   # 生成器检验
python3 scripts/hash_prerun.py
octave --eval "addpath('scripts'); make_models('K')"
octave --eval "addpath('scripts'); run_confirm(m,s)"   # m=1..8, s=1..4
python3 scripts/seal_decisions.py
octave --eval "addpath('scripts'); score()"
python3 scripts/analyze.py && python3 scripts/make_figure.py
```
