# 定方向前的理论推导（2026-11-21）

标签：机制研究（理论推导，没有新的偏移或补偿运行）。

- `DERIVATION_ZH.md`：正文。期望算子的风险分解（定理 1）、相关类指标下的最优算子与实测分辨下限（定理 2）、谱比校准定 σ 的价值上限、随机极性估噪的偏差，以及对分块翻转第一步的事先预测（定理 3）。
- `scripts/t_bayes.py`：定理 1 的蒙特卡罗核对；用 `uq_rtm_20261117/analysis/stage2_C*.json` 的缓存分数估校准上限。
- `scripts/t_floor.py`：相关指标下 Q 误差的分辨下限，以及频散参考频率错位造成的时移。
- `scripts/t_blockflip.py`：分块翻转估噪的几何模型预测（约 6 分钟）。
- `analysis/*.json`：上述脚本的输出。

运行：`D:\ChatGPT\研究\.venv\Scripts\python.exe scripts\<name>.py`。
