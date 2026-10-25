# 从这里开始

这是一个统一调用的研究版本：同一算法与固定规则用于全部模型，资料相关参数自动从观测与校准估计。验收状态见 `core_p1/inputs/validation_status.json`，完整结果见 `RESULTS_ZH.md`，算法说明见 `METHOD_ZH.md`。

## 运行附带示例

在MATLAB中切换到本目录，执行：

```matlab
addpath('core_p1');
result = demo_optimized();
imagesc(result.gather); colorbar;
```

不接收输出参数执行 `demo_optimized`，可显示观测与优化结果的共同色标对比。示例输入在 `core_p1/examples/example_inputs.mat`，只包含观测、校准资料、公开几何和采集方式，不包含真实Q或评分真值。

## 处理自己的资料

```matlab
result = compensate_seismic(D, offsets, public, calibrationD, calibrationOffsets, 'Section');
S = result.gather;
```

- `D`：nt×11的观测道集；`offsets`：11道实际偏移距。同偏移距剖面用Section，偏移距变化的道集用Offset。
- `public`：采样间隔dt、样点数nt、时间向量t、公开offset列表、四层速度vp、前三层厚度h、名义频率nominal_frequency及传播算子数值设置；附带示例展示完整字段。
- `calibrationD`：同一采样轴上的独立nt×11校准道集；`calibrationOffsets`为对应偏移距。当前校准需要已知层位几何。
- 不要向public传入Q、真实主频fpeak、真实波形、密度或真实反射系数；入口会拒绝这些字段。

当前接口基于四层传播模型、11道输入和零相位Ricker子波近似。不同地质模型共用规则，但尚无任意子波、任意层数或实测资料的泛化结论。

## 查看对比

结果包中的 `index.html` 可离线选择每一个道集图；`figures/overall_comparison.png` 汇总整体与弱反射误差。`runs/` 保存原始数据、MATLAB输出、源码快照和运行日志；第一轮失败结果见 `ROUND1_RESULTS_ZH.md`。

代码包可独立运行示例和统一入口。完整实验包还附带开发输入及失败/恢复记录；跨目录复现实验前可运行 `scripts/relocate_inputs.py` 更新当前运行入口的本机路径，历史源码快照与原始记录保持原样。历史复现实验使用MATLAB R2024b；绘图脚本使用NumPy/Matplotlib。
