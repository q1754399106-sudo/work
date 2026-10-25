# P1 MATLAB 修订版（2026-10-21）

这是独立修订副本。原始 ZIP、`p1_review_20261021/` 和历史结果保持不变。

本版以 MATLAB 为唯一数值实现，不再依赖缺失的 `verify_numpy_solver`、SciPy 或旧机器的绝对路径。Python 仅负责启动 MATLAB、记录源码哈希和重算 CSV 验收。

## 修改内容

1. 修正 `uq_case` 返回值；用统一的 `p1_prepare_calibration` 在每次校准后重新构建传播算子和 AVO 投影。
2. 无结构分支在尾端对角正则项非零时使用分块 Cholesky 求解；残差超限不能再返回成功标志。
3. 纯 Q 分支关闭噪声尾端项和额外 AVO 噪声收缩项。
4. 将 Q／主频的 profile-loss 协方差传至求解权重：在 log(Q,f) 空间使用五点正权求积，平均原有经验模态惩罚。零协方差恢复点估计公式。协方差仍含经验方差下限，不能解读为校准过的置信区间。
5. 速度压力测试真实将公开层速度乘以配置中的 1.03，保持厚度不变，用错误速度重新校准并构建算子；真实合成数据保持原速度。不是修改主频冒充速度误差。
6. 增加 `TargetMatched`（Frozen 先验、候选残差目标）和 `CombinedFixedTarget`（候选权重、Frozen 残差目标）两项对照。
7. 全量模式在模型、校准种子、场景、观测种子上完整交叉；增加配置快照、校准与算子变化记录、逐场景指标、源码哈希及数值质量断言。
8. 统一 JSON 对象字段，避免 MATLAB 将不同字段的对象数组解码为 cell 后产生点索引错误。

## 方法解释与适用范围

本版是“校准不确定性平均的经验模态正则化”，不是完整的前向模型误差协方差反演。Q、主频不确定性现在确实影响模态惩罚，但没有据此构造完整的数据误差协方差，也没有宣称实现 block bootstrap。保留旧分支名称以便对照历史 CSV。

修改基于已观察的历史结果，因此配置明确标记 `frozen_before_test=false`；本次数据属于回归与探索性验证，不能重新标为盲测。整体误差、弱反射和 AVO 局部误差仍需分别报告，修复代码不会保证所有局部指标改善。

## 运行

需要 MATLAB R2024b（本机已验证）和 Python 3 的标准库，无需额外 MATLAB 工具箱。输出必须为新目录。

```powershell
python scripts/run_matlab.py --action smoke --output runs/my_smoke --matlab 'D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe'
```

`smoke`：先运行原始基础验证和新增回归检查，再运行 SoftLowQ、1 次校准、1 个观测种子、5 个主场景，加 1 个真实速度误差场景，9 种方法、2 条路线，共 108 行指标。该范围保持与此前五场景结果可比，并补查速度压力测试。

```powershell
python scripts/run_matlab.py --action validate --output runs/my_validation
python scripts/run_matlab.py --action full --output runs/my_full --timeout 21600
```

`full`：9 模型 × 3 校准 × 6 观测种子；Main 5 场景、LowerSNR 2 场景和 4 类压力集各 2 场景；9 方法 × 2 路线，共 43,740 行。耗时显著高于 smoke。本次修订交付不把未运行的全量模式称为已验证结果。

验收独立于进程退出状态：

```powershell
python scripts/evaluate_candidate.py --metrics runs/my_smoke/matlab/metrics.csv --candidate CombinedUncertainty --output runs/my_smoke/acceptance
```

108 行指标对应 6 组观测的 9 方法 × 2 路线，不是 108 次独立实验。检查日志和 `checks.json` 中的残差与求解标志，不能只看进程成功退出。

## 检查与审计

- `validation/original_units/checks.json`：原始公开输入一致性、同模板校准和求解验证。
- `validation/regression_checks.json`：非零尾端正则的无结构求解、纯 Q 消融隔离、不确定性权重敏感性、零协方差退化、重新校准后算子变化。
- `matlab/velocity_stress_*.json`：实际速度值、算子变化及压力条件下校准结果。
- `matlab/executed_spec.json`：本次实际使用的配置，不用全量配置替代小范围执行记录。
- `execution.json`：命令、时长、退出状态及实际执行源码哈希。
