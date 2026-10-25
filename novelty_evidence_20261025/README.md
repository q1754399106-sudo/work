# 创新性证据核查（2026-10-25）

主报告：`RESULTS_ZH.md`；数学条件与证明：`THEORY_ZH.md`；先行研究和基线适配：`METHOD_COMPARISON_ZH.md`。

本轮是已见模型上的描述性比较，不是首次未见模型确认。SC是2023年论文式17结构的加权独立实现，采用FISTA，不是作者原版IRLS代码。候选方法和基线均使用相同Q、源子波及观测；基线只在F1/F2选择参数，再固定评估F3/F4/V4/V5。

## 重现

需要Python 3.12、NumPy和SciPy。当前实际版本见`analysis/environment.json`。捆绑Python若没有SciPy，可在本目录运行：

```powershell
python -m pip install --target cache/deps --no-deps scipy
python scripts/restore_prepared.py
python scripts/run_stage.py theory
python scripts/check_operators.py
python scripts/audit_archives.py
python scripts/run_stage.py dev
python scripts/run_stage.py test
python scripts/make_report.py
python scripts/regional_diagnostic.py
python scripts/finalize_report.py
```

`inputs/prepared/`包含12个实际使用的归档MATLAB重拟合状态及SHA256，用于独立重建收缩结果，避免依赖未上传缓存。其它输入从固定Git提交`198f04b`取出并校验LFS对象；新检出环境应先取得该提交及其LFS对象，例如`git lfs fetch origin 198f04b76d54c5ca77cd7b59df3a5653d6001a3e`。本轮没有重新运行MATLAB稀疏支持发现阶段。

`cache/`仅放可重建缓存和任务依赖，不上传。日志记录1800秒硬超时、PID、命令、源码哈希及退出码；只对本任务进程实施超时。首次评估因空支持边界实现问题被主动停止，修复后重跑；该例外和原日志均保留。

## 审计轨迹

- `PLAN_ZH.md`：计算前初始方案。
- `DEV_ADDENDUM_ZH.md`：开发网格扩展、评估措辞门槛和实现修复。
- `inputs/initial_sources/`、`logs/dev_initial*`：初始开发源码和结果记录。
- `inputs/test_attempt1_sources/`、`logs/test_attempt1*`：空支持修复前的评估尝试。
- `analysis/selection.json`：评估前固定的基线参数。
- `analysis/comparison.csv`、`summary.json`：全部评估场景及主张门槛。
- `analysis/operator_checks.json`：伴随、有限差分梯度、物理坐标与空支持检查。
- `analysis/archive_audit.*`：对原680实现研究的聚合表进行算术核对，不冒充新的680次重跑。

本轮没有改动或中断`fd2d_20261025`，没有接管其后续实验。
