# 研究模型 MATLAB 图册

手机请打开 [GitHub 在线图册](README.md)，可直接浏览全部图片；本页保留模型来源与复现说明。

直接打开 `index.html` 浏览图册，`figures/00_research_models_overview.png` 为四类物理模型的代表总览。所有图由 MATLAB R2024b 绘制，同时导出 PNG（200 dpi）、PDF（矢量文字/线条与原始矩阵图像）及可编辑 FIG。

本图册覆盖 29 个主要模型配置：原始四层对照与 P1 九个地质配置，V1–V5，F1–F4，U1–U4，以及 H1–H6。它覆盖主要建模家族和近期核心验证模型，不声称穷尽历史中所有 G/K/C 等参数变体；实测资料未被解释成已知地下真值。

## 模型与坐标

| 模型 | 展示内容 | 原始依据 |
|---|---|---|
| OriginalControl + P1 九配置 | P 波速度、Q 深度剖面 | `config_model.m`、`p1_experiment_spec.json`；前三层厚度 `h=diff([0,event_t]).*vp(1:3)/2`；底部半空间仅画到 1800 m |
| V1/V2/V4 | 时间域反射系数、无衰减目标、Moderate 观测 | 原始 reflectivity、tr.target、obs.D；没有将反射事件强行解释为深度速度界面 |
| V3/V5 | 速度/Q 剖面、河道/楔状体/含气异常轮廓、目标及观测 | 归档 g 参数及原 `gen_complex_model.m` 地质生成段；未重跑正演、标定或噪声 |
| F1–F4/U1–U4 | 细层速度/Q 剖面及第 50 道测井 | 归档 base 日志 + 原 lateral 规则；2 m 显示网格，按所属细层直接取值，无平滑；只有道号，没有假设横向米制距离 |
| H1–H6 | 二维黏声有限差分速度/Q 剖面 | 原 fd_setup、fd_geology、固定随机种子；4 m 网格和 0.25 m 层内子采样的 Backus 等效参数；仅画物理区，不画 PML |

深度向下为正。速度单位 m/s，Q 无量纲。F/U 和 H 的同变量总览使用共同色标，V3/V5 也共用各自色标。目标与观测成对共用对称振幅色标；没有逐道归一化、自动增益或分位数裁幅。绘图程序遇到超出色标范围的数据即报错，避免静默截断。各图标明坐标，不采用等比例几何，横纵跨度应由坐标轴读取。

F/U 的密集细层以原层属性在 2 m 网格上显示；这属于模型显示采样，不增加原始地质细节。H 图展示的是原二维数值网格等效属性，不能视为未上采样细层的逐点物性。

## 复现

MATLAB 命令窗口：

```matlab
addpath('D:\ChatGPT\研究\地震衰减补偿研究\model_visualization_20261026');
plot_research_models
```

打开并修改图：

```matlab
openfig('D:\ChatGPT\研究\地震衰减补偿研究\model_visualization_20261026\figures\00_research_models_overview.fig','visible');
```

运行会更新本目录生成的图和数据。无需额外工具箱。`snapshot_inputs.py` 从 Git 提交提取输入并验证 LFS 对象，不改变共享稀疏检出；已保存的 `inputs/` 可直接复跑，无需联网。

`SOURCE_MANIFEST.json` 记录源提交、源路径、文件长度和 SHA-256；`model_display_data.mat` 保存绘图数据；`FIGURE_INDEX.json` 列出全部图；`VALIDATION.json` 与 `MATLAB_run.log` 保存运行验证。V3 重建与归档 vp/Q 逐元素一致；H1–H6 的 base 逐字段与归档真值一致。V5 的归档仅有生成参数，按相同代码/种子重建，其几何生成流程由 V3 对照检验支持。

图件采用 scientific-visualization 技能中的数据溯源、共同色标、单位和导出后目视检查原则；没有指定期刊版式要求。
