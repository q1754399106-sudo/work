# 地震衰减补偿（Q 补偿）研究：目录总览

2026-10-19 至 2026-10-26 全部实验都在本文件夹。课题：把被地层吸收衰减的地震道集恢复成无吸收的弹性响应。

- 想看结论：[最终报告](final_report_20261025/REPORT_ZH.md) 和 [勘误](final_report_20261025/ERRATA_ZH.md)。
- 想直接用方法：见下面「怎么选方法」。本研究方法的一行调用入口在 [compensation_package_20261024](compensation_package_20261024/README_ZH.md)（不含 SC）；SC 和 QScan＋SC 的实现在 [novelty_evidence_20261025](novelty_evidence_20261025/) 和 [qscan_sc_20261026](qscan_sc_20261026/RESULTS_ZH.md)。
- 每个目录的标签同时存在 [目录标签.csv](目录标签.csv) 里，可以用表格软件筛选。
- 想看地下模型：[MATLAB 在线模型图册（手机可看）](model_visualization_20261026/README.md)，覆盖 29 个主要模型配置，提供 PNG / PDF / 可编辑 FIG。
- 想看补偿后的道集：[MATLAB 补偿道集在线图册（手机可看）](compensated_gathers_20261026/README.md)，覆盖同一批 29 个模型，并给出 Q 低估场景及自研方法对照。

各目录保留原英文名。约 280 个脚本按英文名互相引用，冻结方案也记录了脚本哈希，改名会让旧脚本无法直接重跑。中文名只写在本页和 CSV 里。

## 标签说明

| 标签 | 含义 |
|---|---|
| **推荐** | 当前推荐使用的方法或组成部分 |
| **不如文献方法** | 在本研究流程内通过了门槛或在用，但在中立数据上不如文献基线 SC |
| **在用·未过门槛** | 后续实验实际在用，但没有通过自己的预登记门槛 |
| **基线** | 用来对比的最初方法 |
| **已取代** | 当时有效，后来被更好的版本吸收或替代，代码仍是后续方法的基础 |
| **已淘汰** | 未通过门槛或效果不如现有方法，不再使用 |
| **机制研究** | 不是新方法，用来解释原因或确定使用规则 |
| **无法判定** | 实测资料上没能分出优劣 |
| **对照实验** | 与文献方法的比较，不提出新方法 |
| **进行中** | 方案已冻结，还在运行 |
| 数据 / 工具 / 图件 / 报告 | 辅助材料，不涉及方法优劣 |

误差都是越低越好。“平方误差”与“相对误差”两种口径不能直接比较，见最终报告摘要。

## 怎么选方法（摘自最终报告）

- **Q 可信**（来自测井、VSP 等）：用文献方法 SC（马雄等 2023 式的结构约束）。在两种方法都没调过参的二维数据 H1–H6 上，SC 12/12 胜过相干标准版，平均相对误差 0.351 对 0.427（`sc_neutral_20261026`）。只在本研究自己的方法里选，就用相干软收缩（标准版或 v1），不要叠加 QScan。
- **Q 可能有误**（只来自地面资料标定）：用 QScan＋SC。在 H1–H6 上它 6/6 胜过 QScan＋相干标准版，平均误差比 0.922（`qscan_sc_20261026`）。一维数据上的复核（`qscan_sc_round3_20261026`）出现分化：Kennett 的 F3/F4 上 QScan＋SC 同样 24/24 更好；但在 V4/V5（相干规则的开发生成器）上，QScan＋相干标准版 24/24 更好。所以这条规则在三种生成器中的两种上成立，其中包括唯一中立的二维数据。拆开三项生成器差别的 24 个新模型（`generator_factors_20261026`）也全部是 SC 一侧更好，用 V4/V5 原始生成代码和新种子生成的 9 个同家族场景（`vfamily_20261027`）上，标准版一侧在 Q 有误和 Q 正确时都稳定更好，所以例外是 V 家族可重复的特性：数据更像 V 家族时宜用 QScan＋标准版。不过 SC 当时的参数是在 Kennett F1/F2 上选的；给 SC 同样的开发预算、在 V1–V3 上重新选参后（`sc_vtuned_20261028`），W 类 Q 有误时差距收回 89%、已无稳定优劣，Q 正确时和 X 类标准版仍稳定更好但只领先约 7–11%，所以 V 家族上的优势一部分来自调参，一部分来自构造。Q 正确时 SC 与 QScan＋SC 没有稳定差别（SC 平均好 9.2%，赢 8/12）。
- 合成数据上，QScan＋相干标准版在 F3/F4 上比最初分步法误差低 55%（28/28 行更好）。这是描述性结果，不是预登记确认。实测资料上至今无法判定优劣。

---

## 一、起点与基线（10-19 → 10-22）

| 中文名 | 目录 | 标签 | 效果 | 优点 |
|---|---|---|---|---|
| 最初研究交接包 | [attenuation_research_handoff_20261019](attenuation_research_handoff_20261019/START_HERE.md) | 基线 | 含默认的 Frozen 分步法；四路对比中最好的组合只比基线低 1.55% | 起点数据、模型和评价工具都在这里 |
| 首轮候选执行与交付框架 | [attenuation_run_20261019](attenuation_run_20261019/README_EXECUTION_ZH.md) | 已取代 | 候选没有替换 Frozen | 建立了“冻结→测试→严格评价→交付”流程和求解器审计脚本，后续沿用 |
| 最初固定参数分步法（FixedSplit）原始归档 | [original_fixed_split_20261022](original_fixed_split_20261022/START_HERE.md) | 基线 | F3/F4 平均误差 1.051，V5 为 0.918 | 简单，参数固定，无需工具箱；复杂模型、中等噪声、Q 低估 30% 时曾优于 FullCalRefit |
| 最初方法在新数据上的补跑 | [fixed_baseline_round3_20261024](fixed_baseline_round3_20261024/RESULTS_ZH.md) | 基线 | 在 F3/F4、V5 上按相同观测重跑 FixedSplit | 让新方法与最初方法同场对比 |

## 二、P1 不确定性加权（10-21）

| 中文名 | 目录 | 标签 | 效果 | 优点 |
|---|---|---|---|---|
| P1 方案审查与历史回放 | [p1_review_20261021](p1_review_20261021/EFFECTS_ZH.md) | 已淘汰 | CombinedUncertainty 在历史 126 对观测上暂定通过（43/43），误差降 6.07%；单独的噪声项、Q 误差项、加强阻尼都未通过 | 第一次实证“Q 不确定性＋噪声不确定性”联合加权有用 |
| P1 MATLAB 修订版 | [p1_fixed_20261021](p1_fixed_20261021/RESULTS_ZH.md) | 已淘汰 | 相干噪声场景误差降 31%，其他场景 −0.7%～+3.9%；整体验收失败，有速度失配 | 修好了校准、求解器的若干缺陷，改成纯 MATLAB 实现 |
| 校准可靠性检查与保守回退 | [p1_guarded_20261021](p1_guarded_20261021/README_ZH.md) | 已淘汰 | 速度 +3% 时能止损，−3% 时过度拒绝 | 校准异常时明确返回“未补偿”，不冒充成功 |
| P1 完整重跑 | [p1_complete_20261021](p1_complete_20261021/README_ZH.md) | 已淘汰 | 9 模型、40,824 行指标；CombinedUncertainty 工程验收未通过（314/328） | 覆盖最全的回归验证，推翻了早先的暂定通过 |
| 校准收缩改进（CalShrink025Light） | [p1_improved_20261022](p1_improved_20261022/SELECTED_DEVELOPMENT_RESULTS_ZH.md) | 已取代 | 开发集上比 Frozen 误差低 9.70%，通过开发门槛，没有独立确认 | Q、主频的几何收缩让剖面域校准更稳，后续被沿用 |

## 三、联合正则与 FullCalRefit（10-22）

| 中文名 | 目录 | 标签 | 效果 | 优点 |
|---|---|---|---|---|
| 联合正则改进方向（文档） | [JOINT_REGULARIZATION_DIRECTIONS_20261022.md](JOINT_REGULARIZATION_DIRECTIONS_20261022.md) | 报告 | 建议先拆清各项职责、解除权重绑定，不要先加更多惩罚项 | 给出下一轮实验的方向 |
| 联合正则重构与交叉消融 | [joint_refactor_20261022](joint_refactor_20261022/PLAN_ZH.md)、[结果](joint_refactor_20261022_results/RESULTS_ZH.md) | 已淘汰 | DecoupledWeak 只比 Legacy 好 0.75%；物理 Q 惩罚反而差 12.10% | DecoupledWeak 解除了权重绑定，成为后面的中间对照 |
| 统一方法（源谱惩罚） | [universal_method_20261022](universal_method_20261022/RESULTS_ZH.md) | 已淘汰 | 三个全局候选都没过预设门槛 | 一套规则、不读模型标签的设计思路 |
| 最初分步法与新方法的道集对比 | [fixed_vs_new_gathers_20261022](fixed_vs_new_gathers_20261022/START_HERE.md) | 已取代 | DecoupledWeak 整体误差 0.854 → 0.681（降 20%），60/60 组更好；但弱反射振幅误差高 10%，AVO 梯度误差高 2% | 同一份观测、同一份校准的公平对比 |
| 数据驱动联合选参，得到 FullCalRefit | [adaptive_joint_cv_20261022](adaptive_joint_cv_20261022/START_HERE.md) | 已取代 | FullCalRefit（反射系数估计＋振幅重拟合）开发集 59/60 通过；剪枝版在稀疏模型上整体误差更低（0.152 对 0.180） | **后续所有方法的核心**：用反射系数估计取代分步滤波 |

## 四、复杂模型与 Q 误差问题（10-23）

| 中文名 | 目录 | 标签 | 效果 | 优点 |
|---|---|---|---|---|
| 更换地质模型对比 | [changed_models_comparison_20261023](changed_models_comparison_20261023/RESULTS_ZH.md) | 已取代 | 4 个新地质模型 80 组：剪枝版波形误差 0.175，FullCalRefit 0.202，FixedSplit 0.838 | 证明换了地质模型仍明显优于最初方法 |
| 统一优化（剪枝版）交付核对 | [delivery_check_20261023](delivery_check_20261023/RESULTS_ZH.md) | 已取代 | PrunedFullCalRefit：开发 100 组 0.152（FixedSplit 0.838），全新 60 组 0.153（FixedSplit 0.862） | 稀疏模型上误差最低，尖灭尾部干净 |
| 复杂类实测地层模型 | [complex_fieldlike_20261023](complex_fieldlike_20261023/START_HERE.md) | 已取代 | FullCalRefit 全道误差比分步法低 11%（中等噪声）、16%（强噪声）；据此选未剪枝版 | 第一次在 60 层复杂模型上验证 |
| 实验版本筛选（文档） | [EXPERIMENT_VERSION_SELECTION_20261023_ZH.md](EXPERIMENT_VERSION_SELECTION_20261023_ZH.md) | 报告 | 复杂地层选未剪枝 FullCalRefit，稀疏模型保留剪枝版作对照 | 把正反两方面证据都列出 |
| 去干扰模块与真实炮集 | [denoise_and_field_20261023](denoise_and_field_20261023/START_HERE.md) | 已淘汰 | 去干扰后 FullCalRefit 全道误差最低（0.547 / 0.644），但只在开发资料上；长孔径检测在 4 个真实炮集上都没触发 | 可作可选前处理 |
| 实测炮集试点 | [field_pilot_20261023](field_pilot_20261023/README.md) | 无法判定 | 局部去干扰减去 0.8–1.3% 能量；Q 可估计性门槛未过，没法比较补偿方法 | 建立了开发炮/留出炮的实测检验流程 |
| Q 误差稳定性 | [q_error_stability_20261023](q_error_stability_20261023/README.md) | 机制研究 | 重拟合法能容忍的 Q 误差区间几乎只有 0%：Q 一偏就失控 | 暴露了关键弱点，直接促成后来的 QScan |
| 99 道高清道集 | [q_dense_gathers_20261023](q_dense_gathers_20261023/README_ZH.md) | 图件 | 不涉及新方法 | 直观看不同 Q 误差下的道集 |
| Q 不确定性感知的反射簇降秩 | [q_identifiable_cluster_20261023](q_identifiable_cluster_20261023/RESULTS_ZH.md) | 已淘汰 | 正常道误差比 FullCalRefit 差 2.40% | 低能量区略好（0.91%） |
| 输出保幅信赖域 | [q_output_trust_20261023](q_output_trust_20261023/RESULTS_ZH.md) | 已淘汰 | 只好 0.18%，未过门槛 | 输出空间约束能防止降秩损失振幅 |
| 带 Q 参考值的实测资料检索 | [q_reference_data_20261023](q_reference_data_20261023/CATALOG_ZH.md) | 数据 | Cape EGS 井中 DAS、Sumatra 自然地震，均附作者 Q | 有外部 Q 参考，可做实测适配 |
| 有外部 Q 参考的自然地震试验 | [q_supported_trial_20261023](q_supported_trial_20261023/README.md) | 无法判定 | 两步岭回归的回放误差低于稀疏版；Cape 频谱检查不支持补偿成功 | 第一次用带 Q 参考的实测资料 |
| 组系数软收缩 | [soft_refit_20261023](soft_refit_20261023/) | 已取代 | 介于未剪枝和硬剪枝之间，保住弱反射 | 相干软收缩的基础 |

## 五、相干收缩与 QScan（10-23 → 10-24）

| 中文名 | 目录 | 标签 | 效果 | 优点 |
|---|---|---|---|---|
| 相干软收缩 v1 | [coherence_soft_20261023](coherence_soft_20261023/) | 不如文献方法 | 5/5 预登记门槛通过（30 个新场景），比软收缩误差低 8.2%；但 Q 可信时不如文献方法 SC | 本研究唯一完整通过门槛的降噪版本；用邻道相干度决定收缩强度，保弱反射 |
| 相干软收缩 v2（标准版 / 稳健版） | [coherence_soft_v2_20261024](coherence_soft_v2_20261024/) | 在用·未过门槛；不如文献方法 | 比 v1 再低 2.1%，没到 3% 门槛（4/5）；后续实验实际用标准版。中立数据 H1–H6 上 SC 12/12 胜过它 | 组合中换成 v1 只差 0.6–1.6%，效果基本等同 |
| 文献方法对比：BAE-Q、QScan、SBL-TVD、Wang 反 Q 滤波 | [bae_q_20261024](bae_q_20261024/) | QScan 选择步骤 **推荐**（接 SC 使用）；QScanRefit 整体已被 QScan＋SC 取代；其余已淘汰 | QScanRefit 在验证集和独立集都通过全部门槛，开发 20 场景比最初方法低 50%，Q ±30% 内误差基本平坦。它不做收缩，H1–H6 上其中的 FullCalRefit（0.569）不如相干标准版（0.427）和 SC（0.351），且从未与 QScan＋SC 直接比较。BAE-Q 在 Q 正确时误差翻倍；SBL-TVD 对 Q 误差更敏感；Wang 反 Q 在 Q 正确时精度差 | **Q 有误时不再失控**；这里还生成了独立 Kennett 细层数据 F1/F2 |
| 推荐版本统一入口 `q_compensate` | [compensation_package_20261024](compensation_package_20261024/README_ZH.md) | 工具（本研究方法的封装） | 与已归档实验逐值一致；入口里不含 SC，两种模式在 H1–H6 上分别不如 SC、QScan＋SC | 一行调用：相干标准版，或加 `'mode','qscan'`；其中的 QScan 步骤可以接到 SC 前面 |
| 文献调研 | [literature_review_20261024](literature_review_20261024/LITERATURE_REVIEW_ZH.md) | 报告 | 本项目每个问题文献里都有成熟对应方法 | 指出 Q 进字典、SBL 自动定参等路线 |
| 阶段进度汇总 | [progress_summary_20261024](progress_summary_20261024/SUMMARY_ZH.md) | 报告 | 截至 10-24 的方法对比图表 | 一页看全前五天 |
| 第三轮：子波相位、QScan＋相干组合、门限 QScan | [round3_20261024](round3_20261024/) | 组合不如文献方法；门限 QScan 已淘汰 | QScan＋相干标准版：F3/F4 平均误差 0.468（最初方法 1.051，降 55%，28/28 行更好），V5 降 41%。描述性结果，登记切换规则曾淘汰该组合。H1–H6 上 QScan＋SC 6/6 胜过它。子波相位只贡献约 2%；门限 QScan 门限过保守 | 本研究自己的方法里 Q 有误时最稳；各道集最坏平均也只有 0.479 |

## 六、QScan 偏差与修正尝试（10-24 → 10-25）

| 中文名 | 目录 | 标签 | 效果 | 优点 |
|---|---|---|---|---|
| QScan 偏差—方差分析 | [qscan_bias_variance_20261024](qscan_bias_variance_20261024/) | 机制研究 | QScan 确定性地选 1.24–1.28 倍 Q0，无噪声时也一样 | 证明偏差是确定的，不是噪声造成 |
| 更换 QScan 选择准则（残差、BIC、AIC 等） | [qscan_criterion_20261024](qscan_criterion_20261024/) | 已淘汰 | 都不比 QScan 好 | 排除了“换准则”这条路 |
| 不依赖稀疏先验的 4 种准则 | [qcriteria_20261024](qcriteria_20261024/RESULTS_ZH.md) | 已淘汰 | 都比 QScan 差 | 确认原 QScan 偏差最小、最稳 |
| 分深度 QScan 上限诊断（第四轮） | [round4_20261024](round4_20261024/RESULTS_ZH.md) | 已淘汰 | 用真值挑逐层因子也只好 5.7%；只用数据时反而差 11.1% | 用上限诊断快速否定一个方向 |
| 按深度分段 Q 校正 | [qscan_layered_20261024](qscan_layered_20261024/RESULTS_ZH.md) | 已淘汰 | 所有分段变体都不如共同因子（误差为其 1.01–1.31 倍） | 预登记执行，结论干净 |
| QScan 偏差的物理来源 | [physics_bias_20261024](physics_bias_20261024/RESULTS_ZH.md) | 机制研究 | 主要来源是密集反射系数与 lasso 稀疏先验失配；细层合并成块状层后偏差降约 60% | 找到偏差原因；多次波、透射损失被排除 |
| 理论预测失真的新指标 | [theory_predict_diagnostic_20261024](theory_predict_diagnostic_20261024/RESULTS_ZH.md) | 已淘汰 | 部分有效（相关 0.583 → 0.650），完整指标未过继续开发门槛 | 反射系数输出敏感性这条诊断线有保留价值 |
| 无真值诊断 B0 | [qscan_bias_diagnostic_20261024](qscan_bias_diagnostic_20261024/) | 已淘汰 | 未过 F2/中等噪声门槛 | 不用真值判断是否该叠加 QScan 的首次尝试 |
| B0M 未见数据确认 | [b0_world_truth_20261024](b0_world_truth_20261024/README.md) | 已淘汰 | 在 U1–U4 上确认 39/39，但只适用于 FullCalRefit 且以 Q 正确为前提；作组合开关不比“始终不叠加”好 | 诊断本身在其适用范围内可靠 |
| 组合的偏差—方差分解与使用规则 | [combo_bias_variance_20261024](combo_bias_variance_20261024/) | 机制研究 | QScan 和相干收缩是互相替代而非互补；Q 有误时组合平方误差好 63–79%，Q 可信时差 21–50%。常数去偏 κ 已由用户决定结束 | 定下了“何时叠加 QScan”的使用规则 |
| 实测 4：U1320A（IODP 308）VSP | [field_vsp_u1320a_20261024](field_vsp_u1320a_20261024/) | 无法判定 | VSP 的 Q 不可识别，两条过井测线结论相反 | 实测检验流程 |
| QScan 偏移在新模型 G1–G8 上 | [qscan_offset_20261025](qscan_offset_20261025/) | 机制研究 | 偏移 −4 到 +4 步，因模型而异；常数去偏不可行 | 两条主要发现之一 |
| 无真值特征预测 QScan 偏移 | [qscan_predictor_20261025](qscan_predictor_20261025/) | 已淘汰 | 13 个规则无一通过 P1 | 排除“预测偏移再修正” |
| “有把握才修正”规则 | [qscan_abstain_20261025](qscan_abstain_20261025/) | 已淘汰 | 新模型 K1–K8 上 A1 通过（平均 −13.4%），A2 未过（最坏 +14.1%）；20 dB 下不安全 | 平均有收益，但最坏情况不可控 |
| 组合层面的准则对照 | [critcombo_20261025](critcombo_20261025/) | 已淘汰 | 7 个准则无一通过 W1；A1、A3 平均持平但最坏差 39–62%（平方误差） | 选择准则路线在组合层面结束 |

## 七、复杂模型瓶颈、创新性与复验（10-25 → 10-26）

| 中文名 | 目录 | 标签 | 效果 | 优点 |
|---|---|---|---|---|
| 真 Q 上限诊断 | [trueq_20261025](trueq_20261025/) | 机制研究 | 换成每道真实衰减后误差降 17.5%：四层 Q 标定是主要瓶颈；子波只贡献 0–1% | 找到复杂模型剩余误差的来源 |
| 共享鬼波标定（SharedGhost） | [calq_20261025](calq_20261025/) | 已淘汰 | 开发集剖面误差降 6.9%，但参数门槛未过 | 指出标定偏差来自模板缺少震源鬼波 |
| SharedGhost 新模型确认 | [calq_confirm_20261025](calq_confirm_20261025/) | 已淘汰 | C1–C3 上 K2–K4 全部未过，最坏差 17.4% | 及时否定，没有采用 |
| 联合事件收缩（组间交叉项） | [joint_event_shrinkage_20261025](joint_event_shrinkage_20261025/) | 已淘汰 | 开发集只好 0.004%，诊断集 0.033%，远低于 2%/3% 门槛 | 排除了“加交叉项”这个创新方向 |
| 创新方向评审 | [innovation_direction_20261025](innovation_direction_20261025/DIRECTION_ZH.md) | 报告 | 文献核对：方法组成部分都有先例 | 明确了能主张的只是两条机制发现 |
| 与文献基线 SC 对照（已见数据） | [novelty_evidence_20261025](novelty_evidence_20261025/) | 对照实验 | 相干标准版未过优势门槛：标准版 / SC 平均误差比 1.055，8 个场景只赢 4 个；在独立 Kennett 数据上 0/4 | 提供了 SC 的独立实现，以及收缩增益的理论分析 |
| 二维黏声有限差分复验 | [fd2d_20261025](fd2d_20261025/RESULTS_ZH.md) | 机制研究（复验通过） | E2–E4 成立，两条机制发现在独立二维正演上复现；E1 偏移不稳定 | 证明发现不是一维 Kennett 生成器的特性 |
| 实测 5：Volve 15/9-F-15A | [field_vsp_volve_20261025](field_vsp_volve_20261025/) | 无法判定 | VSP 的 Q 未过可识别性门槛；井上指标分辨力与阈值相当 | 用了更深、更长的 VSP（191 级、0.79 s） |
| Octave 兼容层 v1 / v2 | [octave_compat_20261025](octave_compat_20261025/README.md)、[v2](octave_compat_v2_20261025/README.md) | 工具 | 在 Octave 上逐位复现 MATLAB 的 rand / randn | 没有 MATLAB 时也能复跑实验 |
| 最终报告与勘误 | [final_report_20261025](final_report_20261025/REPORT_ZH.md) | 报告 | 研究线定稿，两次更正 | 全部结论与出处 |
| 中立数据上相干标准版对 SC | [sc_neutral_20261026](sc_neutral_20261026/RESULTS_ZH.md) | 对照实验·已完成 | 中立数据上 SC 稳定优于相干标准版：12/12 个主场景胜出，平均相对误差 0.351 对 0.427，最坏场景比值 0.868；不带结构项的 SC0（0.402）也更好 | 预登记、两种方法都没在这批数据上调过参，回答了“结果是否依赖生成器” |
| QScan＋SC 对 QScan＋相干标准版 | [qscan_sc_20261026](qscan_sc_20261026/RESULTS_ZH.md) | 对照实验·已完成 | Q 有误时 QScan＋SC 稳定优于 QScan＋相干标准版（6/6，比值 0.922）；Q 正确时 SC 与 QScan＋SC 无稳定差别（SC 平均好 9.2%，8/12） | 预登记；给出 Q 有误时目前最好的组合 |
| 一维数据上 QScan＋SC 对 QScan＋相干标准版（7 档 Q 误差） | [qscan_sc_round3_20261026](qscan_sc_round3_20261026/RESULTS_ZH.md) | 对照实验·已完成 | 按生成器分化：Kennett 的 F3/F4 上 QScan＋SC 稳定更好（24/24，比值 0.804）；V4/V5 上 QScan＋相干标准版稳定更好（24/24，QScan＋SC 比值 1.206）；Q 正确时 SC 与 QScan＋SC 无稳定差别 | 预登记；给出 QScan＋SC 规则的适用范围；外部验证方案第一阶段的正交矩阵 |
| 文献选 Q 准则 Yan 2023 的偏移 | [yan_criterion_20261026](yan_criterion_20261026/RESULTS_ZH.md) | 对照实验·已完成 | 文献准则 Yan 2023 在 0–1 全部 11 个权重下选 Q 偏移都因模型而异（14 个模型，极差 4–9 步），且比 QScan 偏得多（平均 2.5–6.6 步对 2.1 步）；QScan 偏移因模型而异不是 QScan 特有的 | 把第二条机制发现推广到“残差＋稀疏度”这一类准则；两种准则偏移的模型排序高度一致（秩相关 0.83–0.89） |
| 文献选 Q 准则 Cheng 2025、Ma 2018 型的偏移 | [cheng_criterion_20261026](cheng_criterion_20261026/RESULTS_ZH.md) | 对照实验·已完成 | 按摘要实现的 Cheng 浅深一致性准则退化：13/14 个模型选补偿最强的网格端点，因为真值本身深层振幅只有浅层的 9–58%；判定不作解释。7 个准则变体（含 Ma 型 lp）都没有跨模型稳定的偏移，QScan 的选择最接近误差最优点（描述性） | 第二条机制发现的推广范围不变，仍只覆盖“残差＋稀疏度”类准则；真正检验 Cheng 2025 需要原文的振幅处理细节 |
| V 组翻转原因：生成器因子实验 | [generator_factors_20261026](generator_factors_20261026/RESULTS_ZH.md) | 机制研究·已完成 | 24 个新模型（稀疏/密集 × 有无构造 × Kennett/求解器算子）全部是 SC 一侧更好，包括仿 V4 的格子（QScan＋SC/QScan＋标准版约 0.84）；复现不通过，三项都不能解释 V4/V5 上的翻转 | 排除了三项可分离的物理解释；SC 一侧更好的证据扩大到 24 个新模型 |
| V 组翻转原因：场景家族复验 | [vfamily_20261027](vfamily_20261027/RESULTS_ZH.md) | 机制研究·已完成 | 用 V4/V5 原始生成代码和新种子生成 W1–W6、X1–X3；两类在 Q 有误（0.763、23/24；0.876、12/12）和 Q 正确（0.731、12/12；0.883、6/6）时都是标准版一侧稳定更好，翻转重现 | 确认 V4/V5 上的翻转是场景家族特性，不是两个具体场景的偶然结果 |
| V 组翻转原因：SC 在 V 家族上调参后复比 | [sc_vtuned_20261028](sc_vtuned_20261028/RESULTS_ZH.md) | 机制研究·已完成 | SC 在 V1–V3 上按同等预算重新选参（a 1/2、mu 1）后：W 类 Q 有误差距收回 89%（1.310→1.033，无稳定优劣），Q 正确仍是标准版 12/12（1.123）；X 类收回 41–46%，仍是标准版稳定更好（1.077、1.078）；结论“部分解释” | 把 vfamily 的翻转拆成“谁调过参”与“家族构造”两部分；剩余差距是构造优势的上限 |
| 实测 6：Volve 15/9-F-11 T2 | [field_vsp_volve_f11_20261026](field_vsp_volve_f11_20261026/RESULTS_ZH.md) | 无法判定 | 登记的 Q 不可用（剔除规则连锁误删，Q≈10）；事后 Q 为 79/135，但不补偿的输入在走廊叠加上最高，两条走廊叠加相差 0.064（阈值 0.02），F-11 与 F-15A 上 SC 对 QScan＋SC 的方向相反 | 同一工区第二口井；RSVSP 上段 Q≈79（69–90）是三次实测中最接近通过门槛的估计 |
| Stratton 3D 适用性核查 | [field_stratton_20261026](field_stratton_20261026/README.md) | 数据（不适用） | 21 口井都没有声波；VSP 只有一口井的上行波，没有下行波，给不出独立 Q | 记下了不适用的原因，避免重复下载 6.6 GB 炮集 |
| 实测 7：SC / QScan＋SC 的实测适配与应用 | [field_application_20261028](field_application_20261028/RESULTS_ZH.md) | 机制研究·实测应用 | 数据估计子波＋按信噪比限带；两口井与 inline 10141 上频带从约 5–29 Hz 拓宽到 10–62 Hz，井震标定在分辨力内不变，各频带信噪比不降（相干标准版降到 1/2–1/6）；登记门槛 A1 设计错误未过，判定“未证明可用” | 第一次给出实测补偿剖面和井上相关的自助区间（宽 0.2–0.45） |
| 实测 8：Volve 半合成检验（真实噪声＋测井反射系数） | [semisynthetic_volve_20261028](semisynthetic_volve_20261028/RESULTS_ZH.md) | 机制研究 | 补偿比不补偿明显更准（H1 成立）；QScan＋SC-a 与 30 dB 常规反 Q 大体持平；QScan 找不回 Q（H4 不成立），Q 误差是最大误差来源；井上相关只能判断方向 | 第一次在与 Volve 信噪比和频谱衰减一致、但答案已知的数据上检验 QScan＋SC-a 与常规反 Q 滤波 |
| 实测 9：Volve 多井标定 | [multiwell_tie_volve_20261028](multiwell_tie_volve_20261028/RESULTS_ZH.md) | 机制研究·实测应用 | 8 组独立井上补偿后井震相关全部提高（平均 +0.25，区间 0.18–0.32，H1 成立）；但固定在 2.5 s 的时不变滤波同样好，c = 0.77 与 SR 的 c 也分不出；只用测井的检查显示单井地层颜色可使 SR 偏 ±25–50% | 第一次在实测上明确得到“补偿比不补偿好”；井上评价从 2 口井扩到 10 口井、8 组 |
| 实测 10：频谱衰减估 Q（替代 QScan） | [specslope_q_volve_20261028](specslope_q_volve_20261028/RESULTS_ZH.md) | 推荐 | SR 在新数据上 92% 落在 ±25% 内、比 QScan 准（|log 误差| 0.09 对 0.24），补偿零时移相关高 0.06；实测 c≈1.0，VSP Q 基本正确，QScan 的 0.77 偏强（实测 9 显示地层颜色可致约 ±30% 偏差，此判断把握降低） | 用多窗口相干谱衰减独立估计 Q 系数，在新构造、新噪声、新 c 上检验，并用到实测 |
| 实测 11：改进 QScan＋SC（SR 相位＋保守振幅） | [srsc_volve_20261028](srsc_volve_20261028/RESULTS_ZH.md) | 已淘汰 | 新数据上比 QScan＋SC-a 零时移相关高 0.12（0.07–0.17），允许时移也高 0.02，颜色压力组不变差；接近真 Q 的 SC-a；只加 SR 相位校正就能得到大部分改进（实测 12：调参后的 SC-a 和稳定化反 Q 都明显更好，增益来自抵消过强稀疏正则） | 把 Q 的作用拆成振幅和相位：相位按 SR 的 Q 校正，振幅保守（c_SR×1.3）少放大噪声；保留原 SC-a，额外计算很小 |
| 实测 12：改进方案对经典保守基线 | [baselines_volve_20261029](baselines_volve_20261029/RESULTS_ZH.md) | 推荐 | 五项判定全不成立：调好的稳定化反 Q（70 dB，S 0.888）和稀疏项极小的 SC-a（0.864）都比实验 11 改进方案（0.729）高；解耦加在调好的方法上无增益（−0.004、−0.005）；正则/增益参数比 Q 准确度更重要 | 用调好参数的 SC-a、稳定化反 Q、相位反 Q＋时变 Wiener 作基线，检验“振幅、相位用不同 Q”是否算方法贡献 |
| 实测 13：实验 12 方法的实测多井检验 | [field_baselines_volve_20261029](field_baselines_volve_20261029/RESULTS_ZH.md) | 推荐 | 三项支持、一项无法判定、无相反：70 dB 反 Q（+0.033）、小 a 的 SC-a（+0.049）都高于实验 11 改进方案；用户 QScan＋SC-a 调小 a 后 +0.063，组平均 T2 0.523 居首；70 dB 与 30 dB 分不出；所有补偿比不补偿高 0.23–0.30（8/8） | 参数沿用实验 12（70 dB 反 Q、小 a 的 SC-a、Wiener），连同用户的 QScan＋SC-a（原样／调小 a）在 8 组井的 104 条 inline 上整剖面处理，检验半合成结论在实测井上的方向 |

## 八、交付压缩包与其他文件

按原样保存，没有重新打包；内容对应文件名里的阶段和日期，较大的文件用 Git LFS 保存。

| 文件 | 对应阶段 |
|---|---|
| `original_fixed_split_method_20261022.zip` | 最初分步法 |
| `p1_fixed_20261021.zip`、`p1_guarded_20261021.zip`、`p1_complete_20261021_results_full_20261021T145453Z.zip` | P1 各版本 |
| `new_joint_method_complete_20261022.zip`、`joint_refactor_20261022_results.zip`、`fixed_vs_new_gathers_comparison_20261022.zip` | 联合正则阶段 |
| `changed_models_comparison_20261023.zip`、`complex_fieldlike_comparison_20261023.zip`、`original_models_20261023.zip`、`avo_comparison_20261023.zip`、`optimized_seismic_method_20261023.zip`、`optimized_seismic_experiments_20261023.zip`、`experiment_results.pdf` | 复杂模型阶段 |
| `denoise_and_real_seismic_20261023.zip`、`field_pilot_20261023.zip`、`q_reference_data_20261023.zip`、`q_supported_trial_20261023.zip`、`q_trial.pdf` | 去干扰与实测 |
| `q_error_stability_20261023.zip`、`q_error_gathers_20261023.zip`、`q_stability.pdf`、`dense_gathers_hd_20261023.zip` | Q 误差稳定性与图件 |
| `q_identifiable_cluster_experiment_20261023.zip`、`q_output_trust_experiment_20261023.zip`、`soft_refit_results_20261023.zip` | 反射簇、信赖域、软收缩 |
| `prepare_changed_models.py`、`prepare_soft_trial.py` | 两轮实验的输入准备脚本 |
| `UPLOAD_MANIFEST_20261023.csv` | 10-23 研究快照的文件清单（路径相对于当时的仓库根目录，即本文件夹） |
