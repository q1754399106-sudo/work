# 执行记录与复现

工作目录限定 D:\ChatGPT\研究\joint_event_shrinkage_20261025。未改动上游源码、Git配置、其他实验或C盘旧副本；未运行支持发现；未改 Q。提交、推送和清理由主智能体最终处理。

1. 先读取主目录 AGENTS.md 和指定上游代码；利用 git ls-tree、git grep 定位 round3_20261024/run_combo.m 的内嵌 refit_state。完整来源SHA与快照见 inputs/source_manifest.json。
2. 在运行主候选前写 PLAN_ZH.md，冻结候选、开发分割、eta网格和全部门槛。由真实SVD ridge推导全N。
3. 首次数值检查通过；检查输出发现 bootstrap 模块导入时重复写入同内容的源码快照，遂增加 main guard，消除导入写入。该问题不改变任何计算或候选；原检查日志保留。随后增加主审要求的正交负对照，并完成第二次数值检查。之后只增加耗时分项记录，再进入开发；代数和求解器实现未变。
4. 开发阶段只计算四个场景和 eta=[0,.25,.5,1]，选定 eta=.25；开发改善0.003940%，已知未达2%门槛，仍按预登记完成冻结诊断。selection.json落盘后不修改scripts，不增加候选或参数。
5. SC输入从既有已收敛缓存按字节复制到本目录 inputs/sc，逐场景与 ff95a69 提交的开发/评估表核验到1e-12；不是重跑SC。每份输入含原求解元数据，最终输入清单记录全部SHA256，后续复现不依赖原SC缓存。

## 命令（PowerShell）

```powershell
Set-Location -LiteralPath 'D:\ChatGPT\研究\joint_event_shrinkage_20261025'
& 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts\bootstrap.py
& 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts\run_stage.py checks
& 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts\run_stage.py dev
& 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts\run_stage.py frozen
```

上列完整顺序仅用于从零复现。已交付目录含selection.json，为防止覆盖冻结选择，dev阶段会主动拒绝重写。**重跑已冻结结果只需最后一条frozen命令；不要再运行bootstrap覆盖已冻结的common/reference。** 从零开发应在另建复现目录复制本任务后，仅将该新目录的analysis/selection.json移入备份，再按上列顺序执行。bootstrap使用相对目录定位，因此复现目录仍应放在同一仓库下。

全新clone需恢复本任务全部inputs/sc与inputs/prepared的LFS实体，并取得198f04b提交中round3_20261024/dev下F1/F2/F3/F4/V4/V5的观测与truth，以及round3_20261024/results下这12个Q+0场景所需LFS实体。读取器验证每个LFS对象的SHA256与长度，不将指针文件当数据。保留ff95a69与198f04b可读取的Git历史及相应本地LFS对象；若对象缺失，先用git lfs fetch获取指定提交所需对象。inputs/source_manifest.json与analysis/input_manifest.json列出实际使用的来源/输入。

Python为3.12.14，NumPy 2.3.5，SciPy 1.18.1；线程限制为1。NumPy随当前运行时，SciPy从既有 novelty_evidence_20261025/cache/deps 只读加载，该ignored依赖目录未作为本项目成果上传。脚本不安装/修改外部依赖，不声称无依赖自包含。全新环境可使用 `python -m pip install numpy==2.3.5 scipy==1.18.1` 安装到对应Python环境；代码添加的上游deps路径若不存在，正常从Python环境找包。也可用pip的--target安装到该精确deps路径，注意避免混用不一致的已有版本。实际版本见analysis/environment.json，数据/源码内容由哈希清单限定。

冻结结果完成后事后汇总使用 `python analysis/posthoc_aggregation.py`；它只读取已保存检查/波形，未包含于预先冻结的候选scripts集合，独立源码与输入哈希见analysis/posthoc_provenance.json。

每次stage产生 logs/<stage>_<UTC>.log 与 .json，记录精确argv、工作目录、开始/结束时间、退出码和源代码SHA256。selection.json同时冻结PLAN、开发输入和scripts哈希；frozen会拒绝任何哈希变化。缓存原观测/真值从198f04b的Git/LFS读入本目录cache；12份重拟合状态从ff95a69归档复制本目录inputs/prepared，不写上游cache。
