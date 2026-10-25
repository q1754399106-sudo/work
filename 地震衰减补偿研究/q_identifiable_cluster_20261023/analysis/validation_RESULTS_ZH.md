# validation 结果

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-10-23
- Verification Status: ANALYZED
- Version Label: q_cluster_validation_v1

| 模型 | 条件 | 方法 | 正常道误差 | 深部误差 | 低能量区误差 |
|---|---|---|---:|---:|---:|
| FreshFoldLowQ | Moderate | Observed | 1.014591 | 1.016713 | 1.015395 |
| FreshFoldLowQ | Moderate | FullCalRefit | 0.636730 | 0.691915 | 1.042474 |
| FreshFoldLowQ | Moderate | PrunedFullCalRefit | 0.652925 | 0.707900 | 1.098728 |
| FreshFoldLowQ | Moderate | GlobalTSVD | 0.638716 | 0.693994 | 1.050758 |
| FreshFoldLowQ | Moderate | ClusterModesNoQ | 0.649326 | 0.705605 | 1.012727 |
| FreshFoldLowQ | Moderate | QAwareRidge | 0.625821 | 0.679519 | 1.053412 |
| FreshFoldLowQ | Moderate | QAwareCluster | 0.651486 | 0.706709 | 1.012489 |
| FreshFoldLowQ | Strong | Observed | 1.014522 | 1.015377 | 1.017375 |
| FreshFoldLowQ | Strong | FullCalRefit | 0.711230 | 0.752015 | 1.321670 |
| FreshFoldLowQ | Strong | PrunedFullCalRefit | 0.756301 | 0.801440 | 1.568346 |
| FreshFoldLowQ | Strong | GlobalTSVD | 0.715799 | 0.754806 | 1.342925 |
| FreshFoldLowQ | Strong | ClusterModesNoQ | 0.737166 | 0.789943 | 1.307656 |
| FreshFoldLowQ | Strong | QAwareRidge | 0.707015 | 0.747433 | 1.332216 |
| FreshFoldLowQ | Strong | QAwareCluster | 0.735183 | 0.786506 | 1.288811 |
| FreshDenseQReverse | Moderate | Observed | 0.856329 | 0.958194 | 0.860706 |
| FreshDenseQReverse | Moderate | FullCalRefit | 0.471586 | 0.593371 | 0.621641 |
| FreshDenseQReverse | Moderate | PrunedFullCalRefit | 0.475637 | 0.603918 | 0.646042 |
| FreshDenseQReverse | Moderate | GlobalTSVD | 0.472369 | 0.594860 | 0.624102 |
| FreshDenseQReverse | Moderate | ClusterModesNoQ | 0.479038 | 0.600414 | 0.642625 |
| FreshDenseQReverse | Moderate | QAwareRidge | 0.467286 | 0.587178 | 0.622935 |
| FreshDenseQReverse | Moderate | QAwareCluster | 0.482412 | 0.601723 | 0.644472 |
| FreshDenseQReverse | Strong | Observed | 0.864774 | 0.962975 | 0.876841 |
| FreshDenseQReverse | Strong | FullCalRefit | 0.568604 | 0.710895 | 0.824610 |
| FreshDenseQReverse | Strong | PrunedFullCalRefit | 0.590450 | 0.748656 | 0.911351 |
| FreshDenseQReverse | Strong | GlobalTSVD | 0.571768 | 0.717260 | 0.834551 |
| FreshDenseQReverse | Strong | ClusterModesNoQ | 0.581954 | 0.728013 | 0.831657 |
| FreshDenseQReverse | Strong | QAwareRidge | 0.569350 | 0.713881 | 0.834338 |
| FreshDenseQReverse | Strong | QAwareCluster | 0.582984 | 0.731423 | 0.831301 |
| FreshThinLowFrequency | Moderate | Observed | 0.970126 | 1.007360 | 0.896204 |
| FreshThinLowFrequency | Moderate | FullCalRefit | 0.429320 | 0.459252 | 0.647330 |
| FreshThinLowFrequency | Moderate | PrunedFullCalRefit | 0.429854 | 0.458471 | 0.673530 |
| FreshThinLowFrequency | Moderate | GlobalTSVD | 0.429865 | 0.460183 | 0.648101 |
| FreshThinLowFrequency | Moderate | ClusterModesNoQ | 0.437124 | 0.465755 | 0.650237 |
| FreshThinLowFrequency | Moderate | QAwareRidge | 0.427843 | 0.458042 | 0.643685 |
| FreshThinLowFrequency | Moderate | QAwareCluster | 0.442112 | 0.471206 | 0.637989 |
| FreshThinLowFrequency | Strong | Observed | 0.974024 | 1.008136 | 0.932561 |
| FreshThinLowFrequency | Strong | FullCalRefit | 0.536182 | 0.588098 | 0.984982 |
| FreshThinLowFrequency | Strong | PrunedFullCalRefit | 0.542063 | 0.605735 | 0.993511 |
| FreshThinLowFrequency | Strong | GlobalTSVD | 0.537631 | 0.590214 | 0.988111 |
| FreshThinLowFrequency | Strong | ClusterModesNoQ | 0.538886 | 0.591226 | 0.984044 |
| FreshThinLowFrequency | Strong | QAwareRidge | 0.535453 | 0.588102 | 0.981908 |
| FreshThinLowFrequency | Strong | QAwareCluster | 0.539852 | 0.591568 | 0.977958 |

## 等模型权重汇总

| 方法 | 正常道误差 | 低能量区误差 | 深部误差 |
|---|---:|---:|---:|
| Observed | 0.949061 | 0.933180 | 0.994793 |
| FullCalRefit | 0.558942 | 0.907118 | 0.632591 |
| PrunedFullCalRefit | 0.574538 | 0.981918 | 0.654353 |
| GlobalTSVD | 0.561024 | 0.914758 | 0.635220 |
| ClusterModesNoQ | 0.570582 | 0.904824 | 0.646826 |
| QAwareRidge | 0.555461 | 0.911416 | 0.629026 |
| QAwareCluster | 0.572338 | 0.898836 | 0.648189 |

## 预登记工程门槛

{
  "overall_improves_2pct": false,
  "each_model_within_2pct": false,
  "better_than_noQ_1pct": false,
  "better_than_Qridge_1pct": false,
  "low_energy_within_2pct": true,
  "numerical": true,
  "all_pass": false
}

均值与胜负是描述统计；道、同模型噪声条件不视为独立地质模型。求解器通过不等于物理模型或保幅成立。低能量区评价使用真值，只在评分程序内读取。
11/11类统计谬误按本实验范围检查：分模型呈现以避免聚合掩盖退化；不从区域均值推断单反射；生成地质的选择局限、支持选择和同族验证均披露；不按结果筛选/控制案例；诊断基率不适用；不将开发到验证变化归因于重复测量改善；保留全部失败；不做多重探索的显著性宣称；规则与源码冻结；归因限于配对消融；受控计算不涉及观察性反向因果。

原始运行：C:\Users\qq175\Documents\ChatGPT\研究\q_identifiable_cluster_20261023\runs\validation_20261023T065539Z；全部CSV、模式路径和源码哈希随结果保存。