# 未知Q下的补偿参数重调

先读CONCLUSIONS_ZH.md了解结果与代价，完整表格见RESULTS_ZH.md，执行前方案在PLAN_ZH.md。真实Q仅生成合成资料和事后评分；补偿从独立校准观测估计有效Q，只有带误差的近似主频参考。

## MATLAB

```matlab
run_unknown_q('validate')
run_unknown_q('calibrate')
run_unknown_q('tune')
run_unknown_q('search2')
run_unknown_q('test')
```

默认新时间戳输出；tune/search2/test默认读取交付固定路径results/matlab_calibrate、matlab_tune及matlab_search2。重新执行整条流程时请显式给出新的根目录并同步调整这些读取路径，或复制本工程到新的研究目录，避免覆盖交付结果。默认时间戳不表示下游会自动查找最新目录。

完整复跑可用replay_unknown_q.py自动建立不含旧结果的新工程，并顺序执行所有阶段，不必手改读取目录：

```text
python python_src/replay_unknown_q.py --output-root 新工程绝对路径 --engine-python 官方Engine兼容Python路径 --analysis-python 已安装requirements依赖的Python路径
```

该入口仍只连接现有ResearchMATLAB。--analysis-python可省略，此时只执行MATLAB流程。复跑会使用相同冻结种子，属于复现，不算新的独立验证。

Engine脚本只连接已共享的ResearchMATLAB，不启动/关闭实例，不清空工作区或关闭图窗。需要兼容的官方Engine Python，可设置MATLAB_ENGINE_PYTHON：

```text
python python_src/connect_unknown_q.py --action calibrate --output results/matlab_calibrate
python python_src/connect_unknown_q.py --action tune --output results/matlab_tune
python python_src/connect_unknown_q.py --action search2 --output results/matlab_search2
python python_src/connect_unknown_q.py --action test --output results/matlab_test
```

## 真实输入的入口

```matlab
% public只含已知运动学与采样：dt,nt,t,offset,vp,h,fref,nfft,
% lambda_x,lambda_fixed,lambda_grid,cov_floor,nominal_frequency。
% 不得含真实Q、真实fpeak、rho或反射振幅。
out = compensate_unknown_q(D, offsets, public, calibrationD, calibrationOffsets, 'Offset', preset);
S = out.split;
```

calibrationD须是同一区域/吸收背景下可用的校准炮检距道集，并包含本原型所假设的三个运动学事件；本轮合成实验用独立校准道集，不能说成任意资料无需先验的全盲估计。若无独立校准道集，需要另行设计同道集交叉验证或外部约束，不能偷用真实Q替代。

preset在本轮由两次独立噪声观测的对称预测风险选择；实际资料没有重复观测时，本原型的调参步骤不能原样套用。默认Frozen只是未调优参照，不表示稳定最佳。

第一轮300条全选Frozen后，在未运行测试前追加5个独立改变单项参数的候选，再运行300条并按同一规则合并选择；原选择保留。test读取results/matlab_search2的最终选择。追加原因、范围和时点见ADDITIONAL_SEARCH_ZH.md。初始PLAN_ZH.md按原样保留，没有倒改成“十候选都预先指定”。

## 边界

本轮估计单个有效Q，非逐层Q恢复。假定速度/界面几何已知、子波可由共同Ricker形状参数及自由正交分量近似；没有解决任意子波形状、未知速度、实测Q估计或真实NMO/偏移。数据仍为Shuey二项与列时变吸收算子的合成模型。

独立NumPy只核对三套校准数据的Q/主频粗网格目标值；不声称完整补偿网格跨语言复现。补偿求解器沿用已核验版本，Frozen显式参数回归检查单独保留。

图表中英且有图例；全部种子、选参输入和结果保留。任何退化都计入报告。
