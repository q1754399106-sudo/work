# 结果图 / Result figures

`render_results.py` 只读取保存的CSV及JSON，使用本次运行目录 `dependencies/` 内的NumPy、pandas、Matplotlib。它不连接MATLAB，也不调用MATLAB图形。输出目录必须不存在。

```text
python candidate/visualization/render_results.py --metrics test/metrics.csv --examples test/examples --output test_figures --worst-cases test/worst_cases
```

输出14幅中英双语PNG/PDF、可浏览的 `index.html`、`FIGURE_MANIFEST.json`、`ALT_TEXT.json` 与图表源CSV。`--worst-cases` 可省略，默认读取 `examples` 同级的 `worst_cases/`。脚本使用Windows的Microsoft YaHei字体，若不可用则使用SimHei；两者均缺少时停止，并保留明确错误。字体路径和哈希写入清单。

所有原始指标行保留在 `source_data/metrics_all_rows.csv`。统计图保留求解失败但数值有限的得分，无法数值展示的非有限值另行计数，不以零替换。噪声次数不被解释为独立地质样本；图中分位范围仅描述分布，不是置信区间。Main和LowerSNR逐模型展示；全部suite用于完整误差分布、失败比例和最差案例。局部指标的NaN可能表示场景不适用，因此单独标注missing or N/A，原值仍在源表。

道集案例使用预先指定的suite/scenario及最先保存的模型、校准与观测种子，不根据效果选图。未保存的指定案例显示缺失占位，不用其他结果补位。道集和局部窗在同一比较行中使用全部显示数组的共享绝对最大值，禁止振幅裁剪；方法缺失或失败时保留占位面板。频谱是完整Nyquist范围的 `abs(rfft)/N`，不作平滑、滤波或每道归一化。局部时间窗固定为0.45–0.85秒，尖灭尾部固定9–11道。

额外两幅失败案例图按明确的事后评分规则，选取三个已声明候选分支Split中**最大有限RelativeError**，包括得分有限的求解失败。平分时按固定分支顺序及模型、suite、场景和种子排序。它展示实际最差候选、同观测Frozen Split、真值与观测，并绘制固定1/6/11道波形及频谱。案例键、状态、双方误差、完整选择规则和缺失失败数写入清单；无法数值重建的非有限失败继续在CSV和失败比例中保留。配对键在CSV提供时包含 `SNRdB`，避免不同信噪比误配。

This renderer preserves every input row and retains finite scores from failed solver runs. Nonfinite values are explicitly counted and remain missing in the source tables. The figures are descriptive; repeated observation noise is not treated as independent geology. Example selection uses only predesignated suite/scenario and IDs/seeds. Shared amplitude limits include every displayed value, with no clipping or independent panel normalization. Every figure includes labels and a legend/colorbar or an explicit missing-artifact notice.

Syntax and dependency imports are checked separately from rendering. `RENDERED_REQUIRES_VISUAL_INSPECTION` means that PNG/PDF export completed; it does not establish scientific validity or verify the rendered layout. Inspect the actual files before delivery. The scientific-visualization skill informed the workflow; attribution and its bibliographic verification belong in the main research report.
