# development 结果

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-10-23
- Verification Status: ANALYZED
- Version Label: q_cluster_development_v1

| 模型 | 条件 | 方法 | 正常道误差 | 深部误差 | 低能量区误差 |
|---|---|---|---:|---:|---:|
| Development | Moderate | Observed | 0.998957 | 1.020581 | 0.977991 |
| Development | Moderate | FullCalRefit | 0.532683 | 0.628869 | 0.691159 |
| Development | Moderate | PrunedFullCalRefit | 0.550456 | 0.650337 | 0.735334 |
| Development | Moderate | GlobalTSVD | 0.534874 | 0.632265 | 0.694067 |
| Development | Moderate | ClusterModesNoQ | 0.546751 | 0.643691 | 0.694748 |
| Development | Moderate | QAwareRidge | 0.530441 | 0.628247 | 0.691654 |
| Development | Moderate | QAwareCluster | 0.553532 | 0.657131 | 0.696344 |
| Development | Strong | Observed | 1.000259 | 1.022308 | 0.988200 |
| Development | Strong | FullCalRefit | 0.637930 | 0.754791 | 0.974104 |
| Development | Strong | PrunedFullCalRefit | 0.691174 | 0.830012 | 1.054825 |
| Development | Strong | GlobalTSVD | 0.646188 | 0.770070 | 0.983732 |
| Development | Strong | ClusterModesNoQ | 0.656428 | 0.774324 | 0.936082 |
| Development | Strong | QAwareRidge | 0.637266 | 0.755312 | 0.978158 |
| Development | Strong | QAwareCluster | 0.657158 | 0.774712 | 0.930758 |

## 等模型权重汇总

| 方法 | 正常道误差 | 低能量区误差 | 深部误差 |
|---|---:|---:|---:|
| Observed | 0.999608 | 0.983096 | 1.021445 |
| FullCalRefit | 0.585307 | 0.832632 | 0.691830 |
| PrunedFullCalRefit | 0.620815 | 0.895080 | 0.740175 |
| GlobalTSVD | 0.590531 | 0.838899 | 0.701168 |
| ClusterModesNoQ | 0.601589 | 0.815415 | 0.709007 |
| QAwareRidge | 0.583853 | 0.834906 | 0.691779 |
| QAwareCluster | 0.605345 | 0.813551 | 0.715921 |

## 预登记工程门槛

{}

均值与胜负是描述统计；道、同模型噪声条件不视为独立地质模型。求解器通过不等于物理模型或保幅成立。低能量区评价使用真值，只在评分程序内读取。
11/11类统计谬误按本实验范围检查：分模型呈现以避免聚合掩盖退化；不从区域均值推断单反射；生成地质的选择局限、支持选择和同族验证均披露；不按结果筛选/控制案例；诊断基率不适用；不将开发到验证变化归因于重复测量改善；保留全部失败；不做多重探索的显著性宣称；规则与源码冻结；归因限于配对消融；受控计算不涉及观察性反向因果。

原始运行：C:\Users\qq175\Documents\ChatGPT\研究\q_identifiable_cluster_20261023\runs\development_20261023T064919Z；全部CSV、模式路径和源码哈希随结果保存。