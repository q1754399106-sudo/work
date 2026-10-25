# 多炮真实资料试验

先看 `field_pilot_results.pdf` 或 `index.html`。详细结论和限制见 `RESULTS_ZH.md`。

实验共八炮，开发记录 111-114，留出记录 115-118；均来自 Alaska 同一测线。Q 稳定性检查未通过，没有运行实测 Q 补偿。输出为研究候选，不替换原方法。

## 重现

完整 ZIP 包已含原始 SEG-Y `raw/L23535.SGY`、原全局基线 `scripts/fk_baseline.py`。

安装 Python、NumPy、Matplotlib，然后在本目录运行：

```powershell
python run_field_pilot.py
```

程序校验原始文件哈希；首次计划与冻结参数已保存。已存在的冻结参数必须逐字符一致，不允许根据留出结果重新选择参数。重复运行是同参数复现，不增加独立验证样本。

`build_delivery.py` 使用 ReportLab、pypdf、pypdfium2、Pillow 和 Windows 微软雅黑字体生成 PDF、图集与压缩包。

## 主要文件

- `data/record_*.npz`：八炮原始及处理后数据，时间、偏移距、采样间隔。
- `analysis/field_metrics.csv`：所有炮集指标。
- `analysis/reflection_probe_controls.csv`：27 个有限人工反射控制，不是实测真值。
- `analysis/spectral_ratio_measurements.csv`：包括未通过拟合检查的所有逐道表观 Q 测量。
- `analysis/q_summary.csv`、`analysis/results.json`：汇总及拒绝原因。
- `PLAN.json`、`FROZEN_PARAMETERS.json`：开发 / 留出划分、方法、门槛及固定事件假设。
- `VALIDATION_NOTE.md`：首次运行末尾的数值容差修正说明；原失败日志保留。
- `MANIFEST.json`：交付文件哈希；`DELIVERY_CHECK.json`：本地打包检查结果。

图中的时间增益仅用于显示，谱比使用原幅度。减去范数和慢速谱锥功率降低均不能解释为实测信噪比或保幅改善。
