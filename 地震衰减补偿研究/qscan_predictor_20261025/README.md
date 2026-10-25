# 用不依赖真值的特征预测 QScan 的偏移（开发检验）

**状态：已完成，开发门槛 P1 未通过。** 先读 [RESULTS_ZH.md](RESULTS_ZH.md)。方案见 [PLAN_ZH.md](PLAN_ZH.md)，冻结哈希见 `logs/PRE_RUN_HASHES.json`。

上游：`qscan_offset_20261025`。在 QScan＋相干标准版中，QScan 相对误差最优 Q 的偏移是模型的确定属性，范围 −4 到 +4 步，所以常数去偏不可行。本轮检验 6 个不依赖真值的特征能否预测这个偏移。

## 内容

- `scripts/qp_setup.m`、`load_model.m`：路径、10 个开发模型（G1–G8、F1、F2）、3 个场景与种子规则，均与原研究相同。
- `scripts/features.m`：逐实现的 QScan 选择与特征 Z1–Z6。
- `scripts/qscan_select.m`、`full_refit.m`、`estimate_sigmas.m`、`noise_replica.m`：与上游研究中的同名文件逐字相同。
- `scripts/analyze.py`：复现检查、留一模型交叉验证、判定 P1。
- `scripts/make_figure.py`：描述图（评分之后新增）。
- `cache/`（不上传）：从 Git 对象取出的输入快照。

## 运行（Octave 8.4，Linux）

```bash
python3 scripts/snapshot_inputs.py
python3 scripts/hash_prerun.py                    # 冻结，拒绝覆盖
octave --eval "addpath('scripts'); for k=[1 2 3 4 9], features(k); end"   # 另一进程跑 [5 6 7 8 10]
python3 scripts/analyze.py
python3 scripts/make_figure.py
```

需要 `octave_compat_20261025`（由 `qp_setup.m` 自动加入路径）以及 NumPy、pandas、Matplotlib。
