# 第二轮：Q 误差感知的输出保幅信赖域

本目录承接上一轮反射簇降秩失败诊断，在模式删除时增加累计输出空间约束。先查看 `RESULTS_ZH.md` 和 `index.html`；它们在全部运行完成后生成。方法定义与预登记门槛见 `PLAN_ZH.md`。

原始研究目录和上一轮冻结代码未改动。`core/qtrust_refit.m` 为本轮新增估计器；`core/qtrust_controls.m` 为数值控制；`core/run_output_trust_case.m` 为统一实验入口。参考波形由同一观测的完整支持Q加权重估得到，保护预算是近似噪声/模型误差尺度，不等于真实波形保真保证。

## 复现已有输入

需要 MATLAB R2024b、Python 3.12、NumPy、Matplotlib、Pillow。本机MATLAB路径记录在 `run_round2.py`；更换机器需调整可执行文件定位，并保留原源码快照。已有输入不必再生成。

```text
python run_round2.py development
python analyze_round2.py development
python run_round2.py validation
python analyze_round2.py validation
python run_round2.py sensitivity
python summarize_sensitivity.py
python build_report.py
python audit_completion.py
```

每次执行创建新运行目录。`SELECTED.json` 已存在时，重新评分开发集只允许复现相同选择；不能覆盖冻结选择。`prepare_sensitivity.py` 仅首次构建扰动输入时运行，交付包已有输入时跳过。`prepare_data.py validation` 具有防覆盖保护。

若MATLAB安装位置不同，可直接在MATLAB中进入本目录，使用相同冻结核心复算单个已登记案例，而无需修改冻结的Python调度器：

```matlab
addpath('core');
run_output_trust_case('replay_V2DenseLowQ', 'inputs/V2DenseLowQ.json', 'validation', 0.5);
```

输出目录必须是尚不存在的新路径。此处0.5是本轮实际选定的kappa；完整自动化运行记录仍以随包SOURCE_MANIFEST和execution.json为准。

`FROZEN_CANDIDATE.json` 记录验证前的源码、选择、输入和真值哈希；`runs/*/SOURCE_MANIFEST.json` 对应实际运行源码和输入快照；`COMPLETION_AUDIT.json` 区分计算完成与方法晋级。

`analysis/` 保留全部开发候选、逐模型/逐道评分、保护预算使用量、模式删除和Q扰动结果；`runs/` 保留精确源码、原始输出、日志、PID和退出码。`prepared.mat` 是不含真值的计算缓存，方便后续诊断；不是独立实验数据。

新地质配置仍来自同族准1D正演，公开粗四层几何由真属性空间平均构造，不能称为实测测井资料或完整二维粘弹性验证。若利用本轮验证结果继续开发，需将这些模型转为开发/回归资料。
