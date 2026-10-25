# 有外部 Q 参考的真实资料试验

先打开 `q_supported_trial_results.pdf` 或 `index.html`；完整方法、定义和限制见 `RESULTS_ZH.md`。

结论：本次自然地震适配中，两步岭回归的受控额外衰减回放误差低于稀疏 / 剪枝版本。此回放使用匹配正反算子，目标是原实测波形，不是真实无衰减信号。Cape 的简化条件频谱检查未支持补偿成功。未替换原反射补偿入口。

## 复算已有试验

压缩包已经包含 `inputs/trial.json`、7 个核心 MATLAB 快照、原始 SAC、作者 Q/t* 表、Cape Q 及频谱样本。

1. 在解压目录的副本中，将已有 `matlab` 输出文件夹重命名为 `matlab_archived`；原始数据和输入不变。入口拒绝覆盖已有输出。
2. 在 MATLAB 中进入该目录，执行 `run_trial`。本轮使用 R2024b，`-nojvm -singleCompThread`。
3. 安装 Python、NumPy、Matplotlib，运行 `cape_spectral_trial.py`，再运行 `analyze_trial.py`。
4. 需要重建图集 / PDF 时，运行 `build_report.py`；需 ReportLab、pypdf、pypdfium2、Pillow 和 Windows 微软雅黑字体。

`prepare_trial.py` 是首次准备输入的来源脚本，默认引用同一工作区的 `q_reference_data_20261023` 和 `denoise_and_field_20261023`。本包独立复算不需要重新运行准备脚本，直接使用已冻结的 `inputs/trial.json`。

## 输入与方法

- 作者 α=0、Z 分量的 Q 和 t*；没有在试验波形上重新估计 Q。
- SAC 的 P 到时 / Q 与作者 t* 已核对。
- 单路径平稳吸收，偏移阶段为单位阵；未知源波形字典 W=I。属于实测适配，不是原版分层 Ricker 反射模型。
- 全路径 t* 因子 0.75 / 1 / 1.25；额外 t*=0.03 / 0.06 / 0.12 s；新增噪声 0 / 1σ，固定种子。
- 4 条记录来自 2 个地震，重复条件和同事件台站不构成独立地震样本。

## 数据与许可

Sumatra：https://zenodo.org/records/6932666

Cape：https://gdr.openei.org/submissions/1823

原数据 CC BY 4.0。保留原作者及来源。本包新增了窗口截取、带通、模型回放与补偿试验，不将作者估计 Q 宣称为精确物理真值。
