# 联合正则研究重构

独立于 p1_improved_20261022 的新研究工程；旧确认实验未修改。

本轮 Legacy 指相同 CalShrink025Light 配置下的原绑定求解，不能与最早的 Combined 或 Frozen 混称；Old/New 是校准估计的交叉。这里没有重新运行 Frozen 的全网格比较。

- `core_p1/matlab_src/p1_solve_branch.m` 新增第四输出 ctx 与显式 prepare-only 模式；普通调用沿用原求解。
- `core_p1/matlab_src/joint_select.m` 独立 alpha/beta、H 整体尺度归一化、共同物理坐标正则、完整选参残差记录。
- `core_p1/matlab_src/joint_physical_ensemble.m` 校准 Q 附近的八等权分层情景、均值算子及 PSD 的 RQ。
- `core_p1/joint_compensate.m` 真实输入研究入口，强制显式给出 mode 和 beta_factor；返回估计、配置、先导、选参及求解诊断。只接受白名单 public 参数，不能传入真实 Q。

MATLAB 调用示例（D/校准资料由调用者提供）：

```matlab
addpath('C:/Users/qq175/Documents/ChatGPT/研究/joint_refactor_20261022/core_p1');
options=struct('mode','Physical','beta_factor',0.3);
out=joint_compensate(D,offsets,public,calibrationD,calibrationOffsets,'Offset',options);
S=out.gather;
```

这只是研究选项，实验失败时也保留输出和风险，不作为推荐默认。

本轮已完成640行配对试验。新校准下较小固定结构系数版只比同校准 Legacy 改善0.75%，物理Q惩罚版平均误差增加12.10%；均未通过预定工程筛选。原有 `combined_uncertainty_v2` 入口保留，新增 `joint_compensate` 用于研究比较。

复现实验（Python 可执行文件使用本机可用的真实解释器）：

```powershell
& 'C:/Users/qq175/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' joint_refactor_20261022/scripts/run_local.py --entry run_joint_trial --label paired_development --timeout 3600 --run-pointer joint_refactor_20261022/LATEST_RUN.json
& 'C:/Users/qq175/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' joint_refactor_20261022/scripts/analyze_joint.py
```

每次新建时间戳运行，保存实际运行源码、源码 SHA256 清单、命令、PID、退出码、日志和 CSV。设计见 PLAN_ZH.md，结果见 RESULTS_ZH.md。`test_joint_algebra` 使用独立小矩阵稠密解核对四方法、两路线的法方程。
