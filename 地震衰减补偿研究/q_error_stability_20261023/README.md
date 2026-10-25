# 合成资料 Q 误差稳定性

先看 `q_stability_results.pdf` 或 `index.html`，完整数值和指标定义见 `RESULTS_ZH.md`。

本轮停止真实资料实验，使用原反射模型算法。已知真 Q 主模型测试 −50% 至 +50% 的均匀误差，以及两种随深度变化的误差；复杂模型只扰动原已标定的四层 Q，不能将其横轴解释为空间真 Q 误差。

## 文件

- `PLAN.json`：计算前注册的设计、误差网格、评分与容差规则。
- `inputs/primary.json`：真 Q、反射率、固定噪声实现、11 个均匀误差与 2 个深度误差模式。
- `inputs/complex.json`：旧复杂模型固定的 11 道及原纯净参考、低能量掩膜、冻结标定。
- `source/`：65 个未修改的核心 MATLAB 文件；`FROZEN_CORE_HASHES.json` 为哈希。
- `matlab_primary/`、`matlab_complex/`：所有输出、参考、噪声尺度及诊断。
- `analysis/metrics.csv`：全部指标；主模型含弱反射振幅和噪声响应，复杂模型含固定低能量区评分。
- `analysis/depth_pattern_metrics.csv`：非均匀 Q 误差结果。
- `analysis/summary.json`：完成检查、指标与连续通过区间。

## 用冻结输入复算

先复制整个目录，在副本中将两个已有输出目录重命名存档；入口拒绝覆盖旧输出。

在 MATLAB R2024b 中进入此目录：

```matlab
run_q_stability('primary')
run_q_stability('complex')
```

随后用具有 NumPy、Matplotlib 的 Python 运行：

```powershell
python analyze.py
```

如需重建 PDF / 图集 / 压缩包，运行 `build_report.py`，需 ReportLab、pypdf、pypdfium2、Pillow 与 Windows 微软雅黑字体。

本包已含全部冻结输入，不需要再次执行 `prepare.py`。准备脚本引用原工作区此前复杂试验，用于说明来源。

## 解读

真 Q 基准并不保证在含噪资料上最低误差：正则化、分步算子近似和稀疏先验也影响结果。稳定区间基于预先选定的相对误差容差和有限测试网格，不是置信区间或通用安全范围。绝对误差与相对变化应同时阅读。

首次主模型入口的结构体初始化错误已修正，仅修改试验入口，科学参数未改变，失败日志保留。
