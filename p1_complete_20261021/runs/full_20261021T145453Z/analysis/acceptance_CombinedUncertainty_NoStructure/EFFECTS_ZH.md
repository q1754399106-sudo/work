# CombinedUncertainty_NoStructure 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 2916 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 0.69%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.677322|0.657051|2.99|0.709191|0/90|
|LowerSNR|Dev_ModerateQ|0.724330|0.706859|2.41|0.765511|0/90|
|LowerSNR|Holdout_DeepLowQ|0.667303|0.661155|0.92|0.804807|0/90|
|LowerSNR|Holdout_HighContrastReversal|0.694350|0.680817|1.95|0.708335|0/90|
|LowerSNR|Holdout_StrongAtten|0.853844|0.973048|-13.96|2.635167|18/90|
|LowerSNR|Holdout_ThinInterbed|0.737769|0.711510|3.56|0.738114|0/90|
|LowerSNR|SoftLowQ|0.821895|0.909637|-10.68|1.974243|18/90|
|LowerSNR|ThinHighFrequency|0.810706|0.792545|2.24|0.828180|0/90|
|LowerSNR|VelocityReversal|0.661767|0.647593|2.14|0.683309|0/90|
|Main|Dev_DeepContrast|0.565929|0.544154|3.85|0.592282|0/90|
|Main|Dev_ModerateQ|0.635546|0.615638|3.13|0.660922|0/90|
|Main|Holdout_DeepLowQ|0.610412|0.589792|3.38|0.763113|0/90|
|Main|Holdout_HighContrastReversal|0.594396|0.578494|2.68|0.597824|0/90|
|Main|Holdout_StrongAtten|0.803759|0.862350|-7.29|1.758914|18/90|
|Main|Holdout_ThinInterbed|0.634161|0.608611|4.03|0.637151|0/90|
|Main|SoftLowQ|0.763472|0.802251|-5.08|1.459785|18/90|
|Main|ThinHighFrequency|0.740526|0.722007|2.50|0.760952|0/90|
|Main|VelocityReversal|0.560572|0.544962|2.78|0.574407|0/90|
|StressColoredNoise|Dev_DeepContrast|0.599550|0.568651|5.15|0.594643|0/36|
|StressColoredNoise|Dev_ModerateQ|0.683567|0.657424|3.82|0.682238|0/36|
|StressColoredNoise|Holdout_DeepLowQ|0.626961|0.597333|4.73|0.614230|0/36|
|StressColoredNoise|Holdout_HighContrastReversal|0.630106|0.605299|3.94|0.620020|0/36|
|StressColoredNoise|Holdout_StrongAtten|0.965173|1.123500|-16.40|1.674285|18/36|
|StressColoredNoise|Holdout_ThinInterbed|0.660187|0.629857|4.59|0.650538|0/36|
|StressColoredNoise|SoftLowQ|0.904496|0.984197|-8.81|1.339086|18/36|
|StressColoredNoise|ThinHighFrequency|0.753189|0.732806|2.71|0.750737|0/36|
|StressColoredNoise|VelocityReversal|0.597856|0.573047|4.15|0.584745|0/36|
|StressNonRickerPhase|Dev_DeepContrast|0.697046|0.674116|3.29|0.692239|0/36|
|StressNonRickerPhase|Dev_ModerateQ|0.777500|0.752934|3.16|0.781088|0/36|
|StressNonRickerPhase|Holdout_DeepLowQ|0.738732|0.721092|2.39|0.736479|0/36|
|StressNonRickerPhase|Holdout_HighContrastReversal|0.722117|0.711216|1.51|0.718442|0/36|
|StressNonRickerPhase|Holdout_StrongAtten|1.076898|1.226409|-13.88|1.831767|18/36|
|StressNonRickerPhase|Holdout_ThinInterbed|0.738735|0.712272|3.58|0.732239|0/36|
|StressNonRickerPhase|SoftLowQ|1.008869|1.114693|-10.49|1.563206|18/36|
|StressNonRickerPhase|ThinHighFrequency|0.809304|0.787345|2.71|0.803708|0/36|
|StressNonRickerPhase|VelocityReversal|0.701777|0.692528|1.32|0.700527|0/36|
|StressResidualNMO|Dev_DeepContrast|0.627995|0.607726|3.23|0.632768|0/36|
|StressResidualNMO|Dev_ModerateQ|0.690217|0.659150|4.50|0.700346|0/36|
|StressResidualNMO|Holdout_DeepLowQ|0.650529|0.626432|3.70|0.642154|0/36|
|StressResidualNMO|Holdout_HighContrastReversal|0.722035|0.719068|0.41|0.742620|0/36|
|StressResidualNMO|Holdout_StrongAtten|1.001713|1.150335|-14.84|1.747644|18/36|
|StressResidualNMO|Holdout_ThinInterbed|0.656751|0.626243|4.65|0.654580|0/36|
|StressResidualNMO|SoftLowQ|0.934197|1.040391|-11.37|1.470031|18/36|
|StressResidualNMO|ThinHighFrequency|0.744920|0.718712|3.52|0.739859|0/36|
|StressResidualNMO|VelocityReversal|0.694868|0.690434|0.64|0.708584|0/36|
|StressVelocityError|Dev_DeepContrast|0.920253|0.920543|-0.03|0.921561|0/36|
|StressVelocityError|Dev_ModerateQ|12.236479|13.091979|-6.99|19.514555|36/36|
|StressVelocityError|Holdout_DeepLowQ|0.847196|0.844750|0.29|0.847728|0/36|
|StressVelocityError|Holdout_HighContrastReversal|0.850872|0.851147|-0.03|0.867436|0/36|
|StressVelocityError|Holdout_StrongAtten|14.457904|16.332217|-12.96|33.438722|36/36|
|StressVelocityError|Holdout_ThinInterbed|0.710914|0.686242|3.47|0.703297|0/36|
|StressVelocityError|SoftLowQ|14.120586|16.051815|-13.68|31.614998|36/36|
|StressVelocityError|ThinHighFrequency|0.806192|0.790387|1.96|0.801612|0/36|
|StressVelocityError|VelocityReversal|0.850682|0.851571|-0.10|0.853409|0/36|

未通过项目：
- LowerSNR/Holdout_StrongAtten/mean_noninferiority：{"Suite": "LowerSNR", "Model": "Holdout_StrongAtten", "n": 90, "baseline": 0.8538435327156635, "candidate": 0.9730483920402647, "improvement_percent": -13.960972327736476, "baseline_median": 0.7903174650752609, "candidate_median": 0.7901774414177365, "baseline_p90": 1.0992761746676751, "candidate_p90": 1.639206264389255, "baseline_p95": 1.191877494681728, "candidate_p95": 1.7193217378643744, "baseline_worst": 1.25211561641426, "candidate_worst": 2.63516661407356, "baseline_failures": 18, "candidate_failures": 18}
- LowerSNR/SoftLowQ/mean_noninferiority：{"Suite": "LowerSNR", "Model": "SoftLowQ", "n": 90, "baseline": 0.8218947235076453, "candidate": 0.9096371749505904, "improvement_percent": -10.675631432269306, "baseline_median": 0.762669089980514, "candidate_median": 0.7608423879859165, "baseline_p90": 1.0111943640167251, "candidate_p90": 1.455721289355464, "baseline_p95": 1.056320706978729, "candidate_p95": 1.7918552909362613, "baseline_worst": 1.09467958397136, "candidate_worst": 1.97424269582199, "baseline_failures": 12, "candidate_failures": 18}
- LowerSNR/SoftLowQ/failure_count：{"baseline": 12, "candidate": 18, "n": 90}
- Main/Holdout_StrongAtten/mean_noninferiority：{"Suite": "Main", "Model": "Holdout_StrongAtten", "n": 90, "baseline": 0.8037586111625233, "candidate": 0.8623503380889875, "improvement_percent": -7.289716851894079, "baseline_median": 0.687899660693441, "candidate_median": 0.687117793118228, "baseline_p90": 1.2892501099165043, "candidate_p90": 1.6250016515055854, "baseline_p95": 1.3879030639368484, "candidate_p95": 1.694266890210384, "baseline_worst": 1.45637753462654, "candidate_worst": 1.75891395212566, "baseline_failures": 18, "candidate_failures": 18}
- Main/SoftLowQ/mean_noninferiority：{"Suite": "Main", "Model": "SoftLowQ", "n": 90, "baseline": 0.7634720014977895, "candidate": 0.8022512394260035, "improvement_percent": -5.079326792879946, "baseline_median": 0.668367737070901, "candidate_median": 0.668736156031714, "baseline_p90": 1.137851319931971, "candidate_p90": 1.37423998505106, "baseline_p95": 1.202098995205776, "candidate_p95": 1.4188412004945434, "baseline_worst": 1.26524829534274, "candidate_worst": 1.45978495254148, "baseline_failures": 18, "candidate_failures": 18}
- StressColoredNoise/Holdout_StrongAtten/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "Holdout_StrongAtten", "n": 36, "baseline": 0.9651726744680841, "candidate": 1.1234995725850732, "improvement_percent": -16.403997160844263, "baseline_median": 0.912454857072414, "candidate_median": 1.075948560096345, "baseline_p90": 1.2667105774230651, "candidate_p90": 1.597705973480065, "baseline_p95": 1.3491484252610024, "candidate_p95": 1.62622563472886, "baseline_worst": 1.37392313280136, "candidate_worst": 1.67428512881587, "baseline_failures": 18, "candidate_failures": 18}
- StressColoredNoise/SoftLowQ/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "SoftLowQ", "n": 36, "baseline": 0.9044957958406755, "candidate": 0.9841966501575939, "improvement_percent": -8.811633473966696, "baseline_median": 0.877754520303871, "candidate_median": 0.9363028398084, "baseline_p90": 1.150025279201295, "candidate_p90": 1.307210011065425, "baseline_p95": 1.1617222756397925, "candidate_p95": 1.3220855383010874, "baseline_worst": 1.18851627840058, "candidate_worst": 1.33908617858329, "baseline_failures": 18, "candidate_failures": 18}
- StressNonRickerPhase/Holdout_StrongAtten/mean_noninferiority：{"Suite": "StressNonRickerPhase", "Model": "Holdout_StrongAtten", "n": 36, "baseline": 1.0768983805240055, "candidate": 1.2264092979428463, "improvement_percent": -13.88347499845719, "baseline_median": 1.022924865346023, "candidate_median": 1.19556032036815, "baseline_p90": 1.48195089647486, "candidate_p90": 1.7459394412594, "baseline_p95": 1.4916990782616226, "candidate_p95": 1.77829112249784, "baseline_worst": 1.52369627637407, "candidate_worst": 1.8317672117848, "baseline_failures": 18, "candidate_failures": 18}
- StressNonRickerPhase/SoftLowQ/mean_noninferiority：{"Suite": "StressNonRickerPhase", "Model": "SoftLowQ", "n": 36, "baseline": 1.008868960114915, "candidate": 1.1146929031924653, "improvement_percent": -10.489364551913315, "baseline_median": 0.983329240390935, "candidate_median": 1.091902836128452, "baseline_p90": 1.2883764862927798, "candidate_p90": 1.5060973624241951, "baseline_p95": 1.3064663257533426, "candidate_p95": 1.5261832753344875, "baseline_worst": 1.33289968920696, "candidate_worst": 1.56320576275987, "baseline_failures": 18, "candidate_failures": 18}
- StressResidualNMO/Holdout_StrongAtten/mean_noninferiority：{"Suite": "StressResidualNMO", "Model": "Holdout_StrongAtten", "n": 36, "baseline": 1.0017132732092255, "candidate": 1.1503349072406106, "improvement_percent": -14.836744007118963, "baseline_median": 0.955561922817037, "candidate_median": 1.1259166434155357, "baseline_p90": 1.428304907310015, "candidate_p90": 1.68944734692309, "baseline_p95": 1.4526408557919275, "candidate_p95": 1.7190297170430848, "baseline_worst": 1.46656310487474, "candidate_worst": 1.74764370476141, "baseline_failures": 18, "candidate_failures": 18}
- StressResidualNMO/SoftLowQ/mean_noninferiority：{"Suite": "StressResidualNMO", "Model": "SoftLowQ", "n": 36, "baseline": 0.9341967879837275, "candidate": 1.040391070830283, "improvement_percent": -11.367442514521397, "baseline_median": 0.910158440540896, "candidate_median": 1.023462111752997, "baseline_p90": 1.23187413130327, "candidate_p90": 1.4377451193273498, "baseline_p95": 1.2492454491544476, "candidate_p95": 1.45265791239448, "baseline_worst": 1.28032603150826, "candidate_worst": 1.47003051397246, "baseline_failures": 18, "candidate_failures": 18}
- StressVelocityError/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Dev_ModerateQ", "n": 36, "baseline": 12.236478804956636, "candidate": 13.091978821052383, "improvement_percent": -6.991390495027128, "baseline_median": 13.1079701613949, "candidate_median": 13.707374293503701, "baseline_p90": 15.78237672772405, "candidate_p90": 16.9679065563926, "baseline_p95": 16.010052777633923, "candidate_p95": 17.842464738639425, "baseline_worst": 17.3511560774854, "candidate_worst": 19.5145554240714, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Dev_ModerateQ/EventAmplitudeError：{"baseline": 3.237604398825965, "candidate": 3.4917939348385185, "allowance": 0.16188021994129825, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOInterceptError：{"baseline": 2.7282603219423742, "candidate": 3.004702449485712, "allowance": 0.13641301609711873, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOGradientError：{"baseline": 2.3894831285401654, "candidate": 2.833009637763759, "allowance": 0.11947415642700827, "n": 36}
- StressVelocityError/Holdout_StrongAtten/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Holdout_StrongAtten", "n": 36, "baseline": 14.457904392050219, "candidate": 16.332216688055002, "improvement_percent": -12.96392786381535, "baseline_median": 12.676617824936663, "candidate_median": 14.168928688330116, "baseline_p90": 23.9993398175421, "candidate_p90": 27.853920444773948, "baseline_p95": 26.6926891176284, "candidate_p95": 33.429974970216946, "baseline_worst": 26.6926892998688, "candidate_worst": 33.4387215617004, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 3.7553818613177787, "candidate": 4.280212984961468, "allowance": 0.18776909306588896, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOInterceptError：{"baseline": 3.037779074454735, "candidate": 3.4770189797825846, "allowance": 0.15188895372273675, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOGradientError：{"baseline": 2.281710106016304, "candidate": 2.9108314816226324, "allowance": 0.11408550530081521, "n": 36}
- StressVelocityError/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "SoftLowQ", "n": 36, "baseline": 14.120586199193193, "candidate": 16.05181492117317, "improvement_percent": -13.67668944289524, "baseline_median": 13.015832476676675, "candidate_median": 14.86124150253865, "baseline_p90": 21.57603703907565, "candidate_p90": 24.944940308791253, "baseline_p95": 25.36939763507525, "candidate_p95": 31.5479593114605, "baseline_worst": 25.4147687872566, "candidate_worst": 31.6149979080187, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/SoftLowQ/EventAmplitudeError：{"baseline": 4.558573031099418, "candidate": 5.287389418406782, "allowance": 0.2279286515549709, "n": 36}
- StressVelocityError/SoftLowQ/AVOInterceptError：{"baseline": 3.354832433157156, "candidate": 3.956568372704084, "allowance": 0.1677416216578578, "n": 36}
- StressVelocityError/SoftLowQ/AVOGradientError：{"baseline": 2.1350922144952227, "candidate": 2.927470849359556, "allowance": 0.10675461072476114, "n": 36}
- main_equal_model_improvement：{"baseline": 0.6565302490183482, "candidate": 0.6520289186542897, "relative_improvement": 0.0068562421469368084}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
