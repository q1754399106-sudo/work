# 最初方法与新方法：相同观测、多场景道集对比

本包是真实执行后的配对试验：3模型×5场景×2噪声档位×2观测种子，共60组观测、120条方法结果。

## 查看结果

1. 双击 **index.html**，离线切换全部60组道集；不需要联网。
2. 每组可切换色彩道集、黑白波形道集、第6道波形，并下载对应PDF。
3. **figures/overview_OriginalReference_18dB.png** 是预先指定的常规模型五场景总览；其他两个模型也各有一张总览。
4. **RESULTS_ZH.md** 为完整比较报告；**analysis/** 包含逐观测、逐模型、逐场景及局部指标CSV。

所有图同时保留真值、含吸收观测、最初FixedSplit和新DecoupledWeak。完整60组和两个种子均展示，没有按恢复效果挑选样例。

## 比的是哪两个方法

- **FixedSplit**：最初归档的重构算法，先炮检距附加吸收、再深度吸收；两步正则参数均固定为1e-3。不是缺失的毕设原代码。
- **DecoupledWeak**：上一轮选定用于本轮对比的研究配置，结构系数0.3；保留新方法的噪声估计、预处理与结构约束。不是Legacy，也不是Physical。

每组观测只合成一次，两方法共用同一个D及同一份校准所得Q、速度、主频。最初固定分步公式保持不变，但本轮接收共享的校准估计模型；没有让旧方法使用真实Q。

本轮是两套恢复流程的比较，不将差异只归因于新正则或分步形式。两个Dev模型此前用于开发；本轮采用新校准与新观测种子，不是全新独立地质或实测验证。

## 道集图的显示规则

- 色彩图上四幅共用同组振幅色标；下两幅误差图另用两方法共用的误差色标。色标均明确画出，不能跨两种色标直接比较颜色深浅。
- 黑白道集四面板采用同组相同波形增益，不逐道、逐方法归一化；正振幅填黑，保留完整0–1.2s记录。
- 前两个AVO场景横轴为炮检距；Section场景横轴为横向位置，炮检距固定750m，不是实际迁移成像剖面。
- 第6道波形用原始绝对振幅；右侧细节图固定取0.45–1.0s。
- 不滤波、不裁振幅、不平滑或插值。图件说明、尺度和SHA256见figures/manifest.json。

## 原始结果和源码

- runs/paired_gathers_20261022T140103Z/matlab/metrics.csv：全部120行评分。
- 同目录 cases/：全部60组绘图原始数组，JSON和MAT各一份。
- 同目录 calibration_input_*.mat、calibration_*.json：独立校准输入与估计。
- 同目录 checks.json：同观测断言、原方法矩阵复合重放和两方法求解残差。
- runs/下的source/core_p1/：实际运行时的冻结源码；SOURCE_MANIFEST.json记录其哈希。
- core_p1/original_v1/：最初v1归档的13个原文件；三个核心数值函数与当前调用版本逐字节相同，核验见analysis/fairness_checks.json。

## 复现

在MATLAB中切换到解压目录，执行：

```matlab
addpath(fullfile(pwd,'core_p1'));
run_fixed_vs_new(fullfile(pwd,'new_matlab_results'));
```

输出文件夹必须尚不存在。MATLAB入口会生成原始数据和评分；图件由scripts/plot_and_report.py生成。Python自动驱动、快照和日志入口为scripts/run_local.py，其MATLAB可执行文件路径需匹配本机安装；Python绘图依赖版本见scripts/figure_requirements.lock.txt。

本次已完成试验，查看结果不需要重复计算。运行记录中的绝对路径是原机器的来源记录；移动本包后，历史数据和图件仍可按本说明中的相对路径查看。
