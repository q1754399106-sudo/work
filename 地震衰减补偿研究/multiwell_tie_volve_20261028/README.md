# 实验 9（准备中）：Volve 多井标定

目前只完成资料核查，见 `INVENTORY_ZH.md`；方案要等实验 8 的结果后再写并冻结。

- `fetch_volve.sh`、`fetch_list.txt`：从 Volve Data Village 下载 20 口 F 井的组合测井和全工区分层表到 `raw/`；`DATA_LICENCE_ZH.md`：来源与许可；`DOWNLOAD_MANIFEST_SHA256.txt`：校验值。
- `scripts/inventory.py`：读声波、密度覆盖段，用分层点换算双程时间和井位，查离数据体最近道的距离；输出 `inventory/wells.csv` 和 `inventory/logs/*.csv`。
