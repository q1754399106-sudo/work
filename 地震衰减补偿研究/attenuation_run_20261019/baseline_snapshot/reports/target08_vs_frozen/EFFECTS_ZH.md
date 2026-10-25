# Tuned 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 126 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 -6.92%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|SoftLowQ|0.892977|1.176436|-31.74|2.368708|7/12|
|LowerSNR|ThinHighFrequency|0.796037|0.735203|7.64|0.769416|0/12|
|LowerSNR|VelocityReversal|0.675363|0.622337|7.85|0.647312|0/12|
|Main|SoftLowQ|0.763758|1.004657|-31.54|2.855816|12/30|
|Main|ThinHighFrequency|0.735211|0.678583|7.70|0.716286|0/30|
|Main|VelocityReversal|0.559791|0.518034|7.46|0.560500|0/30|

未通过项目：
- LowerSNR/SoftLowQ/mean_noninferiority：{"Suite": "LowerSNR", "Model": "SoftLowQ", "n": 12, "baseline": 0.8929767611181485, "candidate": 1.1764356291161597, "improvement_percent": -31.743140509398682, "baseline_median": 0.8758417134806935, "candidate_median": 1.179179503043715, "baseline_p90": 1.059815075090963, "candidate_p90": 1.622387718266536, "baseline_p95": 1.073449544443812, "candidate_p95": 1.970349066283285, "baseline_worst": 1.08724302222858, "candidate_worst": 2.36870778737757, "baseline_failures": 4, "candidate_failures": 7}
- LowerSNR/SoftLowQ/failure_count：{"baseline": 4, "candidate": 7, "n": 12}
- Main/SoftLowQ/mean_noninferiority：{"Suite": "Main", "Model": "SoftLowQ", "n": 30, "baseline": 0.763758422316647, "candidate": 1.0046567626000358, "improvement_percent": -31.541169726507402, "baseline_median": 0.668457795299802, "candidate_median": 0.7444674885367305, "baseline_p90": 1.134532447400935, "candidate_p90": 1.538861375618679, "baseline_p95": 1.1795736563589974, "candidate_p95": 1.796388458524499, "baseline_worst": 1.24779254442409, "candidate_worst": 2.85581590953063, "baseline_failures": 6, "candidate_failures": 12}
- Main/SoftLowQ/failure_count：{"baseline": 6, "candidate": 12, "n": 30}
- Main/SoftLowQ/WeakEventError：{"baseline": 0.5034851739117435, "candidate": 0.5719460295241808, "allowance": 0.025174258695587176, "n": 6}
- Main/SoftLowQ/FaultNeighborhoodError：{"baseline": 0.6585558833873612, "candidate": 0.7742468207493395, "allowance": 0.03292779416936806, "n": 6}
- Main/SoftLowQ/PinchoutTailError：{"baseline": 0.0008695287878698672, "candidate": 0.0053829900425701295, "allowance": 0.001, "n": 6}
- Main/ThinHighFrequency/WeakEventError：{"baseline": 1.178976534487705, "candidate": 1.355316364278945, "allowance": 0.058948826724385256, "n": 6}
- Main/ThinHighFrequency/PinchoutTailError：{"baseline": 0.006006143020115632, "candidate": 0.012675339615280232, "allowance": 0.001, "n": 6}
- Main/VelocityReversal/WeakEventError：{"baseline": 0.4959318440919577, "candidate": 0.5627478150983601, "allowance": 0.024796592204597886, "n": 6}
- Main/VelocityReversal/PinchoutTailError：{"baseline": 0.023182958910829634, "candidate": 0.03626076813577985, "allowance": 0.0023182958910829637, "n": 6}
- main_equal_model_improvement：{"baseline": 0.6862534005355426, "candidate": 0.7337581548860752, "relative_improvement": -0.06922334273820807}

触边/目标失配诊断：{'scale endpoint': 3, 'residual target mismatch': 1}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
