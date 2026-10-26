# MATLAB 补偿道集图册（手机可直接看）

覆盖上一册的 **29 个主要模型配置**，另给出 H1–H6 的 Q 低估 30% 场景；共 35 个场景、38 张图。

每张主图：**无衰减目标 → 含衰减观测 → 补偿后 → 误差道集**，并附两道原始幅度波形叠加。点击图片或“原图”可放大，无需 MATLAB。

使用项目当前推荐规则：给定参考 Q 用 **SC**，四层标定 Q 不确定或 Q 低估时用 **QScan＋SC**。这些是当前推荐方案，不声称每个模型上都是所有方法中的最低误差；V3/V4/V5 另给出自研方法对照。

[先看 F3](figures/F3.png) · [先看复杂模型 V5](figures/V5.png) · [先看二维模型 H1](figures/H1.png) · [原模型图册](../model_visualization_20261026/README.md)

[模型、数据与复现说明](README_ZH.md) · [全部指标 CSV](metrics.csv) · [MATLAB 绘图源码](plot_gathers.m)

全部场景采用一个固定的名义 14 dB 噪声实现。V3/V5 还含附加相干噪声、多次波等干扰，14 dB 指原始随机噪声设定，并非总干扰信噪比。四层十配置是新生成的 FaultSection 展示观测；其余沿用原观测。F1/F2 为 SC 选参数据。图中相对波形误差为 ‖补偿−目标‖F / ‖目标‖F，越低越好；复杂模型死道保留，全道指标包含死道，另在 CSV 给出仅活道指标。

统一原始幅度色标，无自动增益（AGC）、无逐道归一化、无分位数裁幅；残差也使用同一色标。SC 参数全局固定，没有对本册新模型调参。

## 快速索引

| 模型 / 场景 | 方法 | 观测误差 | 补偿误差 | 图 |
|---|---|---:|---:|---|
| Dev_DeepContrast | SC | 0.981 | 0.428 | [查看](figures/Dev_DeepContrast.png) |
| Dev_ModerateQ | SC | 0.997 | 0.459 | [查看](figures/Dev_ModerateQ.png) |
| F1 | SC | 0.857 | 0.310 | [查看](figures/F1.png) |
| F2 | SC | 0.724 | 0.207 | [查看](figures/F2.png) |
| F3 | SC | 0.848 | 0.413 | [查看](figures/F3.png) |
| F4 | SC | 0.643 | 0.194 | [查看](figures/F4.png) |
| H1 | SC | 0.802 | 0.287 | [查看](figures/H1.png) |
| H2 | SC | 0.494 | 0.311 | [查看](figures/H2.png) |
| H3 | SC | 0.878 | 0.347 | [查看](figures/H3.png) |
| H4 | SC | 0.600 | 0.268 | [查看](figures/H4.png) |
| H5 | SC | 0.897 | 0.422 | [查看](figures/H5.png) |
| H6 | SC | 0.443 | 0.218 | [查看](figures/H6.png) |
| Holdout_DeepLowQ | SC | 0.962 | 0.300 | [查看](figures/Holdout_DeepLowQ.png) |
| Holdout_HighContrastReversal | SC | 0.967 | 0.332 | [查看](figures/Holdout_HighContrastReversal.png) |
| Holdout_StrongAtten | SC | 1.010 | 0.487 | [查看](figures/Holdout_StrongAtten.png) |
| Holdout_ThinInterbed | SC | 0.986 | 0.491 | [查看](figures/Holdout_ThinInterbed.png) |
| OriginalControl | SC | 0.972 | 0.360 | [查看](figures/OriginalControl.png) |
| SoftLowQ | SC | 1.013 | 0.485 | [查看](figures/SoftLowQ.png) |
| ThinHighFrequency | SC | 0.989 | 0.531 | [查看](figures/ThinHighFrequency.png) |
| U1 | SC | 0.823 | 0.307 | [查看](figures/U1.png) |
| U2 | SC | 0.535 | 0.336 | [查看](figures/U2.png) |
| U3 | SC | 0.879 | 0.324 | [查看](figures/U3.png) |
| U4 | SC | 0.605 | 0.249 | [查看](figures/U4.png) |
| V1 | SC | 0.818 | 0.126 | [查看](figures/V1.png) |
| V2 | SC | 0.770 | 0.151 | [查看](figures/V2.png) |
| V3 | QScan+SC | 0.998 | 0.545 | [查看](figures/V3.png) |
| V4 | SC | 0.811 | 0.148 | [查看](figures/V4.png) |
| V5 | QScan+SC | 1.004 | 0.605 | [查看](figures/V5.png) |
| VelocityReversal | SC | 0.960 | 0.338 | [查看](figures/VelocityReversal.png) |
| H1_Qminus30 | QScan+SC | 0.802 | 0.407 | [查看](figures/H1_Qminus30.png) |
| H2_Qminus30 | QScan+SC | 0.494 | 0.304 | [查看](figures/H2_Qminus30.png) |
| H3_Qminus30 | QScan+SC | 0.878 | 0.492 | [查看](figures/H3_Qminus30.png) |
| H4_Qminus30 | QScan+SC | 0.600 | 0.312 | [查看](figures/H4_Qminus30.png) |
| H5_Qminus30 | QScan+SC | 0.897 | 0.481 | [查看](figures/H5_Qminus30.png) |
| H6_Qminus30 | QScan+SC | 0.443 | 0.257 | [查看](figures/H6_Qminus30.png) |

## 一、四层模型

### Dev_DeepContrast · SC

相对波形误差：**0.981 → 0.428**；新生成展示观测。

![Dev_DeepContrast 目标、观测、补偿及误差道集](figures/Dev_DeepContrast.png)

[原图](figures/Dev_DeepContrast.png) · [PDF](figures/Dev_DeepContrast.pdf) · [可编辑 FIG](figures/Dev_DeepContrast.fig)

### Dev_ModerateQ · SC

相对波形误差：**0.997 → 0.459**；新生成展示观测。

![Dev_ModerateQ 目标、观测、补偿及误差道集](figures/Dev_ModerateQ.png)

[原图](figures/Dev_ModerateQ.png) · [PDF](figures/Dev_ModerateQ.pdf) · [可编辑 FIG](figures/Dev_ModerateQ.fig)

### Holdout_DeepLowQ · SC

相对波形误差：**0.962 → 0.300**；新生成展示观测。

![Holdout_DeepLowQ 目标、观测、补偿及误差道集](figures/Holdout_DeepLowQ.png)

[原图](figures/Holdout_DeepLowQ.png) · [PDF](figures/Holdout_DeepLowQ.pdf) · [可编辑 FIG](figures/Holdout_DeepLowQ.fig)

### Holdout_HighContrastReversal · SC

相对波形误差：**0.967 → 0.332**；新生成展示观测。

![Holdout_HighContrastReversal 目标、观测、补偿及误差道集](figures/Holdout_HighContrastReversal.png)

[原图](figures/Holdout_HighContrastReversal.png) · [PDF](figures/Holdout_HighContrastReversal.pdf) · [可编辑 FIG](figures/Holdout_HighContrastReversal.fig)

### Holdout_StrongAtten · SC

相对波形误差：**1.010 → 0.487**；新生成展示观测。

![Holdout_StrongAtten 目标、观测、补偿及误差道集](figures/Holdout_StrongAtten.png)

[原图](figures/Holdout_StrongAtten.png) · [PDF](figures/Holdout_StrongAtten.pdf) · [可编辑 FIG](figures/Holdout_StrongAtten.fig)

### Holdout_ThinInterbed · SC

相对波形误差：**0.986 → 0.491**；新生成展示观测。

![Holdout_ThinInterbed 目标、观测、补偿及误差道集](figures/Holdout_ThinInterbed.png)

[原图](figures/Holdout_ThinInterbed.png) · [PDF](figures/Holdout_ThinInterbed.pdf) · [可编辑 FIG](figures/Holdout_ThinInterbed.fig)

### OriginalControl · SC

相对波形误差：**0.972 → 0.360**；新生成展示观测。

![OriginalControl 目标、观测、补偿及误差道集](figures/OriginalControl.png)

[原图](figures/OriginalControl.png) · [PDF](figures/OriginalControl.pdf) · [可编辑 FIG](figures/OriginalControl.fig)

### SoftLowQ · SC

相对波形误差：**1.013 → 0.485**；新生成展示观测。

![SoftLowQ 目标、观测、补偿及误差道集](figures/SoftLowQ.png)

[原图](figures/SoftLowQ.png) · [PDF](figures/SoftLowQ.pdf) · [可编辑 FIG](figures/SoftLowQ.fig)

### ThinHighFrequency · SC

相对波形误差：**0.989 → 0.531**；新生成展示观测。

![ThinHighFrequency 目标、观测、补偿及误差道集](figures/ThinHighFrequency.png)

[原图](figures/ThinHighFrequency.png) · [PDF](figures/ThinHighFrequency.pdf) · [可编辑 FIG](figures/ThinHighFrequency.fig)

### VelocityReversal · SC

相对波形误差：**0.960 → 0.338**；新生成展示观测。

![VelocityReversal 目标、观测、补偿及误差道集](figures/VelocityReversal.png)

[原图](figures/VelocityReversal.png) · [PDF](figures/VelocityReversal.pdf) · [可编辑 FIG](figures/VelocityReversal.fig)


## 二、Kennett 细层 F1–F4 / U1–U4

### F1 · SC

相对波形误差：**0.857 → 0.310**。

![F1 目标、观测、补偿及误差道集](figures/F1.png)

[原图](figures/F1.png) · [PDF](figures/F1.pdf) · [可编辑 FIG](figures/F1.fig)

### F2 · SC

相对波形误差：**0.724 → 0.207**。

![F2 目标、观测、补偿及误差道集](figures/F2.png)

[原图](figures/F2.png) · [PDF](figures/F2.pdf) · [可编辑 FIG](figures/F2.fig)

### F3 · SC

相对波形误差：**0.848 → 0.413**。

![F3 目标、观测、补偿及误差道集](figures/F3.png)

[原图](figures/F3.png) · [PDF](figures/F3.pdf) · [可编辑 FIG](figures/F3.fig)

### F4 · SC

相对波形误差：**0.643 → 0.194**。

![F4 目标、观测、补偿及误差道集](figures/F4.png)

[原图](figures/F4.png) · [PDF](figures/F4.pdf) · [可编辑 FIG](figures/F4.fig)

### U1 · SC

相对波形误差：**0.823 → 0.307**。

![U1 目标、观测、补偿及误差道集](figures/U1.png)

[原图](figures/U1.png) · [PDF](figures/U1.pdf) · [可编辑 FIG](figures/U1.fig)

### U2 · SC

相对波形误差：**0.535 → 0.336**。

![U2 目标、观测、补偿及误差道集](figures/U2.png)

[原图](figures/U2.png) · [PDF](figures/U2.pdf) · [可编辑 FIG](figures/U2.fig)

### U3 · SC

相对波形误差：**0.879 → 0.324**。

![U3 目标、观测、补偿及误差道集](figures/U3.png)

[原图](figures/U3.png) · [PDF](figures/U3.pdf) · [可编辑 FIG](figures/U3.fig)

### U4 · SC

相对波形误差：**0.605 → 0.249**。

![U4 目标、观测、补偿及误差道集](figures/U4.png)

[原图](figures/U4.png) · [PDF](figures/U4.pdf) · [可编辑 FIG](figures/U4.fig)


## 三、V1–V5 反射与复杂地层

### V1 · SC

相对波形误差：**0.818 → 0.126**。

![V1 目标、观测、补偿及误差道集](figures/V1.png)

[原图](figures/V1.png) · [PDF](figures/V1.pdf) · [可编辑 FIG](figures/V1.fig)

### V2 · SC

相对波形误差：**0.770 → 0.151**。

![V2 目标、观测、补偿及误差道集](figures/V2.png)

[原图](figures/V2.png) · [PDF](figures/V2.pdf) · [可编辑 FIG](figures/V2.fig)

### V3 · QScan+SC

相对波形误差：**0.998 → 0.545**。

![V3 目标、观测、补偿及误差道集](figures/V3.png)

[原图](figures/V3.png) · [PDF](figures/V3.pdf) · [可编辑 FIG](figures/V3.fig)

**同一观测上的自研方法对照（包括表现更好的候选）**

![V3 固定候选方法对照](figures/compare_V3.png)

[对照原图](figures/compare_V3.png) · [PDF](figures/compare_V3.pdf)

### V4 · SC

相对波形误差：**0.811 → 0.148**。

![V4 目标、观测、补偿及误差道集](figures/V4.png)

[原图](figures/V4.png) · [PDF](figures/V4.pdf) · [可编辑 FIG](figures/V4.fig)

**同一观测上的自研方法对照（包括表现更好的候选）**

![V4 固定候选方法对照](figures/compare_V4.png)

[对照原图](figures/compare_V4.png) · [PDF](figures/compare_V4.pdf)

### V5 · QScan+SC

相对波形误差：**1.004 → 0.605**。

![V5 目标、观测、补偿及误差道集](figures/V5.png)

[原图](figures/V5.png) · [PDF](figures/V5.pdf) · [可编辑 FIG](figures/V5.fig)

**同一观测上的自研方法对照（包括表现更好的候选）**

![V5 固定候选方法对照](figures/compare_V5.png)

[对照原图](figures/compare_V5.png) · [PDF](figures/compare_V5.pdf)

### VelocityReversal · SC

相对波形误差：**0.960 → 0.338**；新生成展示观测。

![VelocityReversal 目标、观测、补偿及误差道集](figures/VelocityReversal.png)

[原图](figures/VelocityReversal.png) · [PDF](figures/VelocityReversal.pdf) · [可编辑 FIG](figures/VelocityReversal.fig)


## 四、二维模型 H1–H6：参考 Q

### H1 · SC

相对波形误差：**0.802 → 0.287**。

![H1 目标、观测、补偿及误差道集](figures/H1.png)

[原图](figures/H1.png) · [PDF](figures/H1.pdf) · [可编辑 FIG](figures/H1.fig)

### H2 · SC

相对波形误差：**0.494 → 0.311**。

![H2 目标、观测、补偿及误差道集](figures/H2.png)

[原图](figures/H2.png) · [PDF](figures/H2.pdf) · [可编辑 FIG](figures/H2.fig)

### H3 · SC

相对波形误差：**0.878 → 0.347**。

![H3 目标、观测、补偿及误差道集](figures/H3.png)

[原图](figures/H3.png) · [PDF](figures/H3.pdf) · [可编辑 FIG](figures/H3.fig)

### H4 · SC

相对波形误差：**0.600 → 0.268**。

![H4 目标、观测、补偿及误差道集](figures/H4.png)

[原图](figures/H4.png) · [PDF](figures/H4.pdf) · [可编辑 FIG](figures/H4.fig)

### H5 · SC

相对波形误差：**0.897 → 0.422**。

![H5 目标、观测、补偿及误差道集](figures/H5.png)

[原图](figures/H5.png) · [PDF](figures/H5.pdf) · [可编辑 FIG](figures/H5.fig)

### H6 · SC

相对波形误差：**0.443 → 0.218**。

![H6 目标、观测、补偿及误差道集](figures/H6.png)

[原图](figures/H6.png) · [PDF](figures/H6.pdf) · [可编辑 FIG](figures/H6.fig)


## 五、H1–H6：Q 低估 30%，使用 QScan＋SC

### H1_Qminus30 · QScan+SC

相对波形误差：**0.802 → 0.407**。

![H1_Qminus30 目标、观测、补偿及误差道集](figures/H1_Qminus30.png)

[原图](figures/H1_Qminus30.png) · [PDF](figures/H1_Qminus30.pdf) · [可编辑 FIG](figures/H1_Qminus30.fig)

### H2_Qminus30 · QScan+SC

相对波形误差：**0.494 → 0.304**。

![H2_Qminus30 目标、观测、补偿及误差道集](figures/H2_Qminus30.png)

[原图](figures/H2_Qminus30.png) · [PDF](figures/H2_Qminus30.pdf) · [可编辑 FIG](figures/H2_Qminus30.fig)

### H3_Qminus30 · QScan+SC

相对波形误差：**0.878 → 0.492**。

![H3_Qminus30 目标、观测、补偿及误差道集](figures/H3_Qminus30.png)

[原图](figures/H3_Qminus30.png) · [PDF](figures/H3_Qminus30.pdf) · [可编辑 FIG](figures/H3_Qminus30.fig)

### H4_Qminus30 · QScan+SC

相对波形误差：**0.600 → 0.312**。

![H4_Qminus30 目标、观测、补偿及误差道集](figures/H4_Qminus30.png)

[原图](figures/H4_Qminus30.png) · [PDF](figures/H4_Qminus30.pdf) · [可编辑 FIG](figures/H4_Qminus30.fig)

### H5_Qminus30 · QScan+SC

相对波形误差：**0.897 → 0.481**。

![H5_Qminus30 目标、观测、补偿及误差道集](figures/H5_Qminus30.png)

[原图](figures/H5_Qminus30.png) · [PDF](figures/H5_Qminus30.pdf) · [可编辑 FIG](figures/H5_Qminus30.fig)

### H6_Qminus30 · QScan+SC

相对波形误差：**0.443 → 0.257**。

![H6_Qminus30 目标、观测、补偿及误差道集](figures/H6_Qminus30.png)

[原图](figures/H6_Qminus30.png) · [PDF](figures/H6_Qminus30.pdf) · [可编辑 FIG](figures/H6_Qminus30.fig)
