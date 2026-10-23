# 校准可靠性与保守回退实验

独立目录：`p1_guarded_20261021`。此前的原包、检查版本与修订版本均未修改。

**试验结论：不提升为推荐默认版。** 本目录 API 的实验默认用于复现该策略，不代表部署建议。+3% 条件能够止损，但 −3% 条件出现过度拒绝；完整对比见 `RESULTS_ZH.md`。

本次尝试三个变化：

1. 对 Q／主频搜索触边、优化未收敛、归一化稳健拟合损失过大进行可靠性检查。损失指标为 `2*minimum_loss/(样本数*校准道数)`，门槛 3 是探索性工程设置，不是统计检验。
2. 校准正常时，Offset 道集使用 `NoiseUncertaintyOnly`；Section 剖面保留 `Frozen`，避免为小幅整体改善而改变弱反射与断层恢复。
3. 校准异常时，返回原始观测 `D`，显式标记 `REJECTED_CALIBRATION_PASSTHROUGH` 和 `CompensationApplied=0`。这不是成功补偿，也没有恢复被衰减的信息。

默认 API `p1_compensate` 使用 `GuardedNoise`，并在拒绝时跳过传播算子构建和反演。显式请求旧分支仍允许运行，用于比较失败情况。道集数据与校准诊断决定分支，不读取真值或重建误差进行选择。

## 执行范围

SoftLowQ 一模型、校准种子 20261030；观测种子 20261101 与新增的 20261107；5 个主场景；另加 +3% 与 −3% 速度条件各一个 SmoothAVO 场景。速度改变时厚度固定，包含垂直走时与校准失配，所有相关算子重新构建。

14 组观测 × 5 方法 × 2 路线，共 140 行指标。方法为 Frozen、NoiseUncertaintyOnly、CombinedUncertainty、GuardedNoise、NoCompensation。新增种子和负速度条件提供额外检查，但仍是同一模型族下的探索性实验；规则设计已参考先前结果，不能声称盲测或独立外部验证。

```powershell
python scripts/run_matlab.py --action guard_trial --output runs/my_trial --matlab 'D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe'
```

Python 仅使用标准库启动 MATLAB、记录执行及源码哈希。执行前会运行已有回归检查和新增拒绝／回退检查。运行结果、实际配置与每种压力条件的校准诊断另存。

CSV 中记录 `CalibrationReliable`、`CalibrationFitDeviance`、`FallbackUsed`、`CompensationApplied`、`EffectiveMethod`、`PipelineStatus`。恒等返回的残差 0 对应 `I*S=D`，不表示衰减反演残差为 0。验收脚本新增补偿覆盖检查：有回退记录不能被判为全部观测成功补偿。

## 如何理解结果

需要同时看恢复误差、局部保幅和实际补偿覆盖率。若异常输入上误差下降只是因为停止错误反演，应描述为“避免失控”，不能描述为“解决速度误差反演”。Q 分支是否值得保留，以其相对仅噪声分支的增益判断；本试验不会根据真值在线选择最优分支。
