# Q 不确定性感知的反射簇自适应反演

这是按“密集反射簇不整组删除，而按 Q 误差下的可辨识组合选择自由度”方向实现的首轮研究实验。它是待检验的创新候选，不是已确认首创的方法或生产默认版。

先看 `RESULTS_ZH.md`（运行和评分结束后生成），图件入口为 `index.html`。

## 方法文件

- `core/qcluster_refit.m`：低秩模型误差白化、簇内组合模式、单调复杂度目标和联合振幅重估。
- `core/qcluster_calibration_cov.m`：独立校准资料的局部 logQ 敏感性，剖除事件幅度/相位并考虑主频 nuisance 参数。
- `core/run_qcluster_case.m`：原 FullCalRefit / 原整组剪枝与四个配对候选的统一批处理。
- `core/qcluster_controls.m`：重复字典压缩、Q误差为零的退化检查；已知支持仅用于数值控制，不进入性能评价。
- `core/matlab_src/`：118个原 MATLAB 文件中的相关数值依赖，原文件清单及哈希见 `FROZEN_BASELINE_HASHES.json`。

## 复现已冻结输入

需要 MATLAB R2024b，以及 Python 3.12、NumPy、Matplotlib、Pillow。MATLAB 路径在 `run_experiment.py` 中为本机安装位置；换机器时仅调整可执行文件定位，并保留原冻结源码版本。

```text
python run_experiment.py development
python analyze.py development
python run_experiment.py validation
python analyze.py validation
python run_experiment.py sensitivity
python summarize_sensitivity.py
python build_report.py
python audit_completion.py
```

各次执行都建立新的 `runs/` 目录。不要覆盖已有输入、真值或冻结清单；`prepare_data.py validation` 和 `prepare_sensitivity.py` 具有防覆盖保护，交付包已有输入时不需要再次运行它们。

## 证据与适用范围

`SPEC.json` 和 `PLAN_ZH.md` 登记模型、空间抽样、对照和晋级门槛；`FROZEN_CANDIDATE.json` 在新模型验证前冻结源码、输入和真值；运行内 `SOURCE_MANIFEST.json` 对应实际执行版本。

新模型来自与既有复杂模型相同的准1D正演族。校准需要假定已知的粗层位几何（本实验来自真属性空间平均）；不存在实测或完整二维波动方程验证。未确定的创新性与实际效果必须分别判断。

图件使用同一案例内共同振幅尺度，波形不逐道归一化；完整CSV、未改善结果、运行日志均保留。Python运行本身使用本机依赖路径，不包含外部数据上传。
