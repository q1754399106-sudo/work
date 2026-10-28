# 实验 9：Volve 多井标定

问题：在 8 组互相独立的井上，用 SR 的 Q 做常规反 Q 补偿后，井震相关是否高于不补偿。方案在 `PLAN_ZH.md`，冻结时的哈希在 `PRE_RUN_HASHES.txt`。

- `INVENTORY_ZH.md`：资料核查（13e2fa7）。`fetch_volve.sh`、`fetch_list.txt`、`DATA_LICENCE_ZH.md`、`DOWNLOAD_MANIFEST_SHA256.txt`、`scripts/inventory.py`：下载、许可和核查程序。
- `scripts/snapshot_inputs.py`：把输入从登记的提交恢复到被忽略的 `cache/inputs/`，并做校验。
- `scripts/wq.py`：测井时深关系、Gardner 反射系数、井轨迹，以及开发阶段的 WR 估计器。
- `scripts/synth.py`、`dev.py`、`devsum.py`：开发阶段的半合成井。结果在 `analysis/dev/`、`analysis/dev_summary.json`，日志在 `logs/dev_*.log`。
- `scripts/tie.py`：多井标定（主检验 H1），输出 `results/tie.json`。
- `scripts/analyze.py`：判定，输出 `analysis/summary.json`。
- `scripts/colour.py`：只用测井的 SR 颜色偏差检查，输出 `results/colour.json`。
