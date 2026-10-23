# 新方法完整包：先读这里

包含本轮联合正则重构代码、640行配对实验、旧/新校准交叉对照、图件、求解与入口检查，以及实际运行时的冻结源码。

## 先看结果

- RESULTS_ZH.md：完整结论、各方法误差、局部指标和工程筛选。
- figures/：整体对比与局部配对图，均有PNG、PDF和说明。
- analysis/：逐模型汇总、局部失败、物理惩罚诊断及检查点完整性。
- evidence/paired_continuation_20261022T113517Z/matlab/metrics.csv：合并后的全部640行，含首批保留的500行与续跑140行。
- evidence/：两批实验及代数/入口检查的命令、退出状态、日志、源码清单和冻结源码。

本轮新校准下，DecoupledWeak 相对同校准 Legacy 的平均误差降低0.75%；Physical 增加12.10%。新候选均未通过预定工程门槛。它们是研究选项，不能把示例配置当作推荐默认。

## 方法和入口

core_p1/joint_compensate.m 是新增研究入口，需要观测D、炮检距、公共参数public、独立校准资料，以及明确的采集域和研究配置。

```matlab
% 先将MATLAB当前文件夹切换到本包根目录。
addpath(fullfile(pwd,'core_p1'));
options=struct('mode','DecoupledWeak','beta_factor',0.3);
out=joint_compensate(D,offsets,public,calibrationD,calibrationOffsets,'Offset',options);
S=out.gather;
```

D等变量需要由调用者提供。当前实现要求601点模型配置下经验证的11道布局；完整输入检查见入口源码。public使用白名单字段，不得加入真实Q等评分信息；Q从校准资料估计。Section域需要明确传入'Section'。

本轮配置：

| 名称 | 定义 | beta_factor |
|---|---|---:|
| Legacy | CalShrink025Light配置下的原绑定对照，保留combined_uncertainty_v2入口 | 原绑定规则 |
| DecoupledWeak | H整体归一化，结构系数固定、与模态强度分开 | 0.3 |
| Decoupled | 与上一配置同一构造，固定结构系数较大 | 1 |
| SharedIsotropic | 共同物理坐标正则与均值算子，不加RQ | 0.3 |
| Physical | 与SharedIsotropic相同，再加入八个分层Q情景构造的RQ | 0.3 |

Legacy不是最初固定参数分步方法；本包与最初方法包应分开使用。旧/新校准使用相同正则构造规则，数据与校准变化后，数值权重依规则更新。

## 自带合成试验入口

无需准备外部D，可以在MATLAB中运行：

```matlab
addpath(fullfile(pwd,'core_p1'));
run_joint_trial(fullfile(pwd,'new_trial_results'));
```

输出文件夹必须尚不存在。完整640行试验在本机曾累计耗时约79分钟；本次打包未重新启动实验。

独立小矩阵核对和研究入口核对分别为test_joint_algebra、test_joint_api；均接受尚不存在的输出目录。本包保留了本轮已完成的检查证据，无需为查看结果重复运行。

## 关于路径和历史记录

README_ZH.md记录原工作区的示例路径；原始日志、RESULTS_ZH.md和LATEST_RUN.json中也保留了实际执行时的绝对路径，它们是来源记录。解压到其他位置后，查看历史结果请使用本包的evidence/和analysis/，不要依赖旧路径。

如使用Python驱动脚本，run_local.py中的MATLAB可执行文件路径需与本机安装位置一致；绘图依赖版本见scripts/figure_requirements.lock.txt。数值实验本身使用MATLAB，Python执行与汇总脚本使用标准库。

打包完整性与文件哈希见MANIFEST.json。本包增加本说明，算法源码和已完成的实验数值均保持原样。
