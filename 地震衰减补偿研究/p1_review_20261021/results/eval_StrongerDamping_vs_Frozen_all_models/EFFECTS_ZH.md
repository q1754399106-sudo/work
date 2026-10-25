# StrongerDamping 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 246 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 -2.90%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|DevMain|Dev_DeepContrast|0.555567|0.594036|-6.92|0.620268|0/15|
|DevMain|Dev_ModerateQ|0.623227|0.658733|-5.70|0.692330|0/15|
|LowerSNR|SoftLowQ|0.950922|0.913912|3.89|1.139515|4/12|
|LowerSNR|ThinHighFrequency|0.789360|0.816200|-3.40|0.825497|0/12|
|LowerSNR|VelocityReversal|0.668010|0.697298|-4.38|0.709497|0/12|
|Main|SoftLowQ|0.773074|0.779261|-0.80|1.185782|6/30|
|Main|ThinHighFrequency|0.729228|0.757219|-3.84|0.778205|0/30|
|Main|VelocityReversal|0.555080|0.580489|-4.58|0.617773|0/30|
|MultiCalMain|SoftLowQ|0.895024|0.847101|5.35|1.317269|2/10|
|MultiCalMain|ThinHighFrequency|0.736741|0.762912|-3.55|0.790135|0/10|
|MultiCalMain|VelocityReversal|0.550469|0.576603|-4.75|0.612230|0/10|
|NewGeologyMain|Holdout_DeepLowQ|0.608064|0.620857|-2.10|0.772537|0/15|
|NewGeologyMain|Holdout_HighContrastReversal|0.585202|0.613607|-4.85|0.642425|0/15|
|NewGeologyMain|Holdout_StrongAtten|0.884106|0.868773|1.73|1.633246|3/15|
|NewGeologyMain|Holdout_ThinInterbed|0.624575|0.663052|-6.16|0.685097|0/15|

未通过项目：
- DevMain/Dev_DeepContrast/mean_noninferiority：{"Suite": "DevMain", "Model": "Dev_DeepContrast", "n": 15, "baseline": 0.5555668623913042, "candidate": 0.5940360612921255, "improvement_percent": -6.9243148763840745, "baseline_median": 0.549017068621116, "candidate_median": 0.589293196179889, "baseline_p90": 0.5797920780690694, "candidate_p90": 0.6186333908820032, "baseline_p95": 0.5819320112784258, "candidate_p95": 0.6200576550132723, "baseline_worst": 0.585230233916147, "candidate_worst": 0.620268439787053, "baseline_failures": 0, "candidate_failures": 0}
- DevMain/Dev_DeepContrast/EventAmplitudeError：{"baseline": 0.420858700635473, "candidate": 0.4761989195434105, "allowance": 0.021042935031773652, "n": 15}
- DevMain/Dev_DeepContrast/WeakEventError：{"baseline": 0.3603009485759657, "candidate": 0.400560991441756, "allowance": 0.018015047428798284, "n": 3}
- DevMain/Dev_DeepContrast/FaultNeighborhoodError：{"baseline": 0.554433488107266, "candidate": 0.5985105206457847, "allowance": 0.0277216744053633, "n": 3}
- DevMain/Dev_DeepContrast/AVOInterceptError：{"baseline": 0.43007012294736047, "candidate": 0.471440564032795, "allowance": 0.021503506147368026, "n": 6}
- DevMain/Dev_DeepContrast/AVOGradientError：{"baseline": 0.46940032957964267, "candidate": 0.5092649730050784, "allowance": 0.023470016478982135, "n": 6}
- DevMain/Dev_ModerateQ/mean_noninferiority：{"Suite": "DevMain", "Model": "Dev_ModerateQ", "n": 15, "baseline": 0.6232270560216668, "candidate": 0.6587327378657055, "improvement_percent": -5.697070032659868, "baseline_median": 0.626045224838796, "candidate_median": 0.65594697746645, "baseline_p90": 0.6659203374297996, "candidate_p90": 0.689202618652147, "baseline_p95": 0.6688007105805127, "candidate_p95": 0.6908331108404827, "baseline_worst": 0.671160680427789, "candidate_worst": 0.692329601352306, "baseline_failures": 0, "candidate_failures": 0}
- DevMain/Dev_ModerateQ/EventAmplitudeError：{"baseline": 0.5399860686998048, "candidate": 0.5844107727822216, "allowance": 0.02699930343499024, "n": 15}
- DevMain/Dev_ModerateQ/WeakEventError：{"baseline": 0.535253903402684, "candidate": 0.5634706766149646, "allowance": 0.0267626951701342, "n": 3}
- DevMain/Dev_ModerateQ/FaultNeighborhoodError：{"baseline": 0.613503929088807, "candidate": 0.656723765459964, "allowance": 0.03067519645444035, "n": 3}
- DevMain/Dev_ModerateQ/AVOInterceptError：{"baseline": 0.558857441253484, "candidate": 0.5908342426741955, "allowance": 0.027942872062674203, "n": 6}
- DevMain/Dev_ModerateQ/AVOGradientError：{"baseline": 0.6153958272613365, "candidate": 0.6497841042698692, "allowance": 0.030769791363066824, "n": 6}
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
- MultiCalMain/SoftLowQ/EventAmplitudeError：{"baseline": 0.4728685813820803, "candidate": 0.5431914067022839, "allowance": 0.023643429069104018, "n": 10}
- MultiCalMain/SoftLowQ/AVOInterceptError：{"baseline": 0.517327248578216, "candidate": 0.5665031369245105, "allowance": 0.0258663624289108, "n": 4}
- MultiCalMain/SoftLowQ/AVOGradientError：{"baseline": 0.605060998461212, "candidate": 0.6538734776707348, "allowance": 0.0302530499230606, "n": 4}
- MultiCalMain/ThinHighFrequency/mean_noninferiority：{"Suite": "MultiCalMain", "Model": "ThinHighFrequency", "n": 10, "baseline": 0.7367413506371769, "candidate": 0.7629117562910369, "improvement_percent": -3.5521836301473098, "baseline_median": 0.735395474015435, "candidate_median": 0.7612650155230685, "baseline_p90": 0.7517426952601529, "candidate_p90": 0.7780217529270276, "baseline_p95": 0.7588454247493779, "candidate_p95": 0.7840784116008873, "baseline_worst": 0.765948154238603, "candidate_worst": 0.790135070274747, "baseline_failures": 0, "candidate_failures": 0}
- MultiCalMain/VelocityReversal/mean_noninferiority：{"Suite": "MultiCalMain", "Model": "VelocityReversal", "n": 10, "baseline": 0.5504688036637976, "candidate": 0.5766033803651397, "improvement_percent": -4.747694424715121, "baseline_median": 0.5437912255586885, "candidate_median": 0.5708541834545935, "baseline_p90": 0.5902780097626055, "candidate_p90": 0.611707041009962, "baseline_p95": 0.5906141951787637, "candidate_p95": 0.6119687347133075, "baseline_worst": 0.590950380594922, "candidate_worst": 0.612230428416653, "baseline_failures": 0, "candidate_failures": 0}
- MultiCalMain/VelocityReversal/EventAmplitudeError：{"baseline": 0.4136944185357116, "candidate": 0.4560638950840862, "allowance": 0.020684720926785582, "n": 10}
- MultiCalMain/VelocityReversal/FaultNeighborhoodError：{"baseline": 0.5349012080567885, "candidate": 0.5793559646124715, "allowance": 0.02674506040283943, "n": 2}
- MultiCalMain/VelocityReversal/AVOInterceptError：{"baseline": 0.46920885002901924, "candidate": 0.49724985040473024, "allowance": 0.023460442501450962, "n": 4}
- NewGeologyMain/Holdout_DeepLowQ/mean_noninferiority：{"Suite": "NewGeologyMain", "Model": "Holdout_DeepLowQ", "n": 15, "baseline": 0.6080635027111908, "candidate": 0.6208565673803732, "improvement_percent": -2.1039027358395357, "baseline_median": 0.626496011685541, "candidate_median": 0.650213805408801, "baseline_p90": 0.7455588832135221, "candidate_p90": 0.7665549374296278, "baseline_p95": 0.751236124241705, "candidate_p95": 0.7717644384008345, "baseline_worst": 0.753943775607822, "candidate_worst": 0.772537282778505, "baseline_failures": 0, "candidate_failures": 0}
- NewGeologyMain/Holdout_DeepLowQ/WeakEventError：{"baseline": 0.1434549630476593, "candidate": 0.212122589454489, "allowance": 0.007172748152382966, "n": 3}
- NewGeologyMain/Holdout_HighContrastReversal/mean_noninferiority：{"Suite": "NewGeologyMain", "Model": "Holdout_HighContrastReversal", "n": 15, "baseline": 0.5852019654660128, "candidate": 0.613606697232773, "improvement_percent": -4.853833965533716, "baseline_median": 0.579513989873878, "candidate_median": 0.611908296089393, "baseline_p90": 0.6132063879191676, "candidate_p90": 0.6370523758775574, "baseline_p95": 0.6139269611592918, "candidate_p95": 0.6394872741543003, "baseline_worst": 0.614165228685408, "candidate_worst": 0.642425293954422, "baseline_failures": 0, "candidate_failures": 0}
- NewGeologyMain/Holdout_HighContrastReversal/EventAmplitudeError：{"baseline": 0.44913930507433075, "candidate": 0.4942082042719209, "allowance": 0.02245696525371654, "n": 15}
- NewGeologyMain/Holdout_HighContrastReversal/FaultNeighborhoodError：{"baseline": 0.586825223293971, "candidate": 0.6257709655052237, "allowance": 0.029341261164698554, "n": 3}
- NewGeologyMain/Holdout_HighContrastReversal/AVOInterceptError：{"baseline": 0.4979074242475095, "candidate": 0.5302017095853169, "allowance": 0.024895371212375478, "n": 6}
- NewGeologyMain/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 0.5060567977737731, "candidate": 0.5814986546251819, "allowance": 0.025302839888688658, "n": 15}
- NewGeologyMain/Holdout_StrongAtten/AVOInterceptError：{"baseline": 0.5064186329684983, "candidate": 0.5625365477031097, "allowance": 0.025320931648424916, "n": 6}
- NewGeologyMain/Holdout_StrongAtten/AVOGradientError：{"baseline": 0.6080300364778738, "candidate": 0.6621022517401178, "allowance": 0.03040150182389369, "n": 6}
- NewGeologyMain/Holdout_ThinInterbed/mean_noninferiority：{"Suite": "NewGeologyMain", "Model": "Holdout_ThinInterbed", "n": 15, "baseline": 0.624575450612793, "candidate": 0.6630523795305954, "improvement_percent": -6.160493320711113, "baseline_median": 0.627387951646851, "candidate_median": 0.668057105662356, "baseline_p90": 0.6405421959378956, "candidate_p90": 0.6791878191618566, "baseline_p95": 0.6428412071868024, "candidate_p95": 0.6812377637007815, "baseline_worst": 0.64643368126862, "candidate_worst": 0.685096972640728, "baseline_failures": 0, "candidate_failures": 0}
- NewGeologyMain/Holdout_ThinInterbed/EventAmplitudeError：{"baseline": 0.5351884293839958, "candidate": 0.584571658732514, "allowance": 0.026759421469199792, "n": 15}
- NewGeologyMain/Holdout_ThinInterbed/FaultNeighborhoodError：{"baseline": 0.6249763950455043, "candidate": 0.6695927697668796, "allowance": 0.031248819752275215, "n": 3}
- NewGeologyMain/Holdout_ThinInterbed/AVOInterceptError：{"baseline": 0.5300507337882416, "candidate": 0.569202509604954, "allowance": 0.026502536689412083, "n": 6}
- NewGeologyMain/Holdout_ThinInterbed/AVOGradientError：{"baseline": 0.5544148575163168, "candidate": 0.5948737233733545, "allowance": 0.027720742875815843, "n": 6}
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.7056564092745623, "relative_improvement": -0.028963053171668474}
- candidate_solver_quality："Direct and requested route included"

触边/目标失配诊断：{'residual target mismatch': 97}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
