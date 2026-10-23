# StrongerDamping 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 126 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 -2.90%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|SoftLowQ|0.950922|0.913912|3.89|1.139515|4/12|
|LowerSNR|ThinHighFrequency|0.789360|0.816200|-3.40|0.825497|0/12|
|LowerSNR|VelocityReversal|0.668010|0.697298|-4.38|0.709497|0/12|
|Main|SoftLowQ|0.773074|0.779261|-0.80|1.185782|6/30|
|Main|ThinHighFrequency|0.729228|0.757219|-3.84|0.778205|0/30|
|Main|VelocityReversal|0.555080|0.580489|-4.58|0.617773|0/30|

未通过项目：
- LowerSNR/SoftLowQ/EventAmplitudeError：{"baseline": 0.6542885219054589, "candidate": 0.716846434577576, "allowance": 0.03271442609527295, "n": 12}
- LowerSNR/SoftLowQ/AVOInterceptError：{"baseline": 0.633898381311064, "candidate": 0.6959856805770381, "allowance": 0.031694919065553204, "n": 12}
- LowerSNR/SoftLowQ/AVOGradientError：{"baseline": 0.7147761837975675, "candidate": 0.7706087984920512, "allowance": 0.03573880918987837, "n": 12}
- LowerSNR/ThinHighFrequency/mean_noninferiority：{"Suite": "LowerSNR", "Model": "ThinHighFrequency", "n": 12, "baseline": 0.7893604330538123, "candidate": 0.8161997830666345, "improvement_percent": -3.400138756510551, "baseline_median": 0.7897390251155469, "candidate_median": 0.816285870353409, "baseline_p90": 0.8010175943595527, "candidate_p90": 0.8250500328850957, "baseline_p95": 0.8012588619886777, "candidate_p95": 0.8252523030396148, "baseline_worst": 0.801399317315715, "candidate_worst": 0.825497031169318, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/VelocityReversal/mean_noninferiority：{"Suite": "LowerSNR", "Model": "VelocityReversal", "n": 12, "baseline": 0.6680101337703328, "candidate": 0.6972982437132886, "improvement_percent": -4.384381083809319, "baseline_median": 0.6676537018741515, "candidate_median": 0.697109313265769, "baseline_p90": 0.6800993614090992, "candidate_p90": 0.7085951858414025, "baseline_p95": 0.6806823238278298, "candidate_p95": 0.709027827539214, "baseline_worst": 0.681388951402094, "candidate_worst": 0.709496684586262, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/VelocityReversal/EventAmplitudeError：{"baseline": 0.5674835016032596, "candidate": 0.6154198159235614, "allowance": 0.02837417508016298, "n": 12}
- LowerSNR/VelocityReversal/AVOInterceptError：{"baseline": 0.5574330685312989, "candidate": 0.5939901195301037, "allowance": 0.027871653426564946, "n": 12}
- Main/SoftLowQ/EventAmplitudeError：{"baseline": 0.515668396337273, "candidate": 0.5797351032780286, "allowance": 0.02578341981686365, "n": 30}
- Main/SoftLowQ/AVOInterceptError：{"baseline": 0.5484923196893657, "candidate": 0.5916744003819857, "allowance": 0.027424615984468288, "n": 12}
- Main/SoftLowQ/AVOGradientError：{"baseline": 0.6331500335347208, "candidate": 0.675928219283517, "allowance": 0.03165750167673604, "n": 12}
- Main/ThinHighFrequency/mean_noninferiority：{"Suite": "Main", "Model": "ThinHighFrequency", "n": 30, "baseline": 0.7292277060186926, "candidate": 0.7572188032381204, "improvement_percent": -3.8384577256737273, "baseline_median": 0.729663196758499, "candidate_median": 0.756498767374693, "baseline_p90": 0.7497770984826575, "candidate_p90": 0.7758731784124848, "baseline_p95": 0.751430280341858, "candidate_p95": 0.7775976906477714, "baseline_worst": 0.751644974949284, "candidate_worst": 0.778205339298211, "baseline_failures": 0, "candidate_failures": 0}
- Main/ThinHighFrequency/FaultNeighborhoodError：{"baseline": 0.7235684502794169, "candidate": 0.761404059530193, "allowance": 0.03617842251397085, "n": 6}
- Main/VelocityReversal/mean_noninferiority：{"Suite": "Main", "Model": "VelocityReversal", "n": 30, "baseline": 0.5550795714794274, "candidate": 0.5804892090668835, "improvement_percent": -4.577656770854133, "baseline_median": 0.546325639483688, "candidate_median": 0.5807065747497665, "baseline_p90": 0.5902649159019921, "candidate_p90": 0.6155225922750431, "baseline_p95": 0.5935314511286296, "candidate_p95": 0.6160012226764574, "baseline_worst": 0.595894416691639, "candidate_worst": 0.617772776950928, "baseline_failures": 0, "candidate_failures": 0}
- Main/VelocityReversal/EventAmplitudeError：{"baseline": 0.4169279969608545, "candidate": 0.4585611740566997, "allowance": 0.020846399848042724, "n": 30}
- Main/VelocityReversal/FaultNeighborhoodError：{"baseline": 0.5431739818503365, "candidate": 0.5816089607235718, "allowance": 0.027158699092516825, "n": 6}
- Main/VelocityReversal/AVOInterceptError：{"baseline": 0.46417704930116116, "candidate": 0.4933831410884042, "allowance": 0.02320885246505806, "n": 12}
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.7056564092745623, "relative_improvement": -0.028963053171668474}
- candidate_solver_quality："Direct and requested route included"

触边/目标失配诊断：{'residual target mismatch': 40}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
