# 交给下一位AI的执行提示

请接手本包任务书。先验证原版，再按P1实现误差模型；不要重新询问已有参数，也不要启动MATLAB新实例。

```text
python scripts/verify_package.py
python scripts/summarize_existing.py
python scripts/run_matlab.py --action check --output 新检查目录 --engine-lib 官方Engine的lib目录
python scripts/run_matlab.py --action validate --output 新验证目录 --engine-lib 官方Engine的lib目录
python scripts/run_matlab.py --action regression --output 新回归目录 --engine-lib 官方Engine的lib目录
python scripts/run_matlab.py --action replay --output 新复跑目录 --engine-lib 官方Engine的lib目录
```

此机器Python：C:/Users/qq175/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe
此机器Engine lib：C:/Users/qq175/Documents/Codex/2026-10-17/new-chat/matlab_engine_runtime/lib
MATLAB安装：D:/Program Files (x86)/MATLAB2024B。只需连接共享ResearchMATLAB；这些路径用于定位，不写死在核心求解器里。

`core/frozen_compensate.m`：public输入+独立校准资料的生产原型，variant仅允许Baseline或Frozen，同时输出direct/split。无实测输入时不伪造实测结论。
`core/run_research.m`：validate沿用小规模数学/输入隔离检查；regression复算一个保存的完整案例并检查路径和参数一致性；replay复跑固定历史输入的四方法。

候选实现放在新的候选目录，不修改evidence。明确复用/替换核心函数，不让MATLAB路径遮蔽未知版本。
新实验CSV保持历史字段，增加CalibrationSeed；Baseline/Frozen与Candidate均保留Direct/Split，同一观测键唯一。

```text
python scripts/evaluate_candidate.py --metrics 新实验metrics.csv --candidate Candidate --baseline Frozen --route Split --output 新评价目录
```

代码字段说明：Method是候选名称，Route是Direct或Split，Model/Suite/Scenario/Seed/CalibrationSeed为配对键；仅后验诊断字段HTruthRatio允许读真值，选参禁止。
该命令自动生成EFFECTS_ZH.md、effects.html、acceptance.json与paired_differences.csv。若缺相同先验Direct、缺局部场景或配对不全，报告不能给出通过。
输出报告必须另外交代RUN_METADATA.json中的新地质独立性、校准种子重复层级、源码哈希、数据来源和真值隔离审计。数值门槛通过不等于研究完成。

先读 `reports/current_effects.html` 了解起点，不把历史Target08当成当前默认。不要为了过门槛改变测试集、丢弃失败或事后换指标。
