# QModelErrorOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 2916 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 2.82%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.677322|0.642804|5.10|0.698597|0/90|
|LowerSNR|Dev_ModerateQ|0.724330|0.694537|4.11|0.743311|0/90|
|LowerSNR|Holdout_DeepLowQ|0.667303|0.653898|2.01|0.798299|0/90|
|LowerSNR|Holdout_HighContrastReversal|0.694350|0.671157|3.34|0.701291|0/90|
|LowerSNR|Holdout_StrongAtten|0.853844|0.814771|4.58|1.184408|11/90|
|LowerSNR|Holdout_ThinInterbed|0.737769|0.704699|4.48|0.734378|0/90|
|LowerSNR|SoftLowQ|0.821895|0.793099|3.50|1.100618|9/90|
|LowerSNR|ThinHighFrequency|0.810706|0.788498|2.74|0.824591|0/90|
|LowerSNR|VelocityReversal|0.661767|0.639967|3.29|0.675821|0/90|
|Main|Dev_DeepContrast|0.565929|0.539604|4.65|0.582891|0/90|
|Main|Dev_ModerateQ|0.635546|0.614131|3.37|0.671568|0/90|
|Main|Holdout_DeepLowQ|0.610412|0.601907|1.39|0.756212|0/90|
|Main|Holdout_HighContrastReversal|0.594396|0.574793|3.30|0.600718|0/90|
|Main|Holdout_StrongAtten|0.803759|0.782267|2.67|1.467910|18/90|
|Main|Holdout_ThinInterbed|0.634161|0.608599|4.03|0.644330|0/90|
|Main|SoftLowQ|0.763472|0.756427|0.92|1.319847|18/90|
|Main|ThinHighFrequency|0.740526|0.721648|2.55|0.766407|0/90|
|Main|VelocityReversal|0.560572|0.542682|3.19|0.578561|0/90|
|StressColoredNoise|Dev_DeepContrast|0.599550|0.569807|4.96|0.602276|0/36|
|StressColoredNoise|Dev_ModerateQ|0.683567|0.661617|3.21|0.695039|0/36|
|StressColoredNoise|Holdout_DeepLowQ|0.626961|0.614183|2.04|0.658433|0/36|
|StressColoredNoise|Holdout_HighContrastReversal|0.630106|0.606760|3.71|0.628995|0/36|
|StressColoredNoise|Holdout_StrongAtten|0.965173|0.905150|6.22|1.240554|18/36|
|StressColoredNoise|Holdout_ThinInterbed|0.660187|0.631000|4.42|0.658677|0/36|
|StressColoredNoise|SoftLowQ|0.904496|0.881935|2.49|1.152887|15/36|
|StressColoredNoise|ThinHighFrequency|0.753189|0.734107|2.53|0.754559|0/36|
|StressColoredNoise|VelocityReversal|0.597856|0.575530|3.73|0.596488|0/36|
|StressNonRickerPhase|Dev_DeepContrast|0.697046|0.672791|3.48|0.700013|0/36|
|StressNonRickerPhase|Dev_ModerateQ|0.777500|0.758624|2.43|0.796258|0/36|
|StressNonRickerPhase|Holdout_DeepLowQ|0.738732|0.737551|0.16|0.778539|0/36|
|StressNonRickerPhase|Holdout_HighContrastReversal|0.722117|0.705449|2.31|0.724014|0/36|
|StressNonRickerPhase|Holdout_StrongAtten|1.076898|1.047607|2.72|1.503761|18/36|
|StressNonRickerPhase|Holdout_ThinInterbed|0.738735|0.713139|3.46|0.741649|0/36|
|StressNonRickerPhase|SoftLowQ|1.008869|1.012587|-0.37|1.390221|18/36|
|StressNonRickerPhase|ThinHighFrequency|0.809304|0.788004|2.63|0.810442|0/36|
|StressNonRickerPhase|VelocityReversal|0.701777|0.687254|2.07|0.706909|0/36|
|StressResidualNMO|Dev_DeepContrast|0.627995|0.596226|5.06|0.629400|0/36|
|StressResidualNMO|Dev_ModerateQ|0.690217|0.663054|3.94|0.711347|0/36|
|StressResidualNMO|Holdout_DeepLowQ|0.650529|0.642506|1.23|0.691908|0/36|
|StressResidualNMO|Holdout_HighContrastReversal|0.722035|0.709530|1.73|0.732132|0/36|
|StressResidualNMO|Holdout_StrongAtten|1.001713|0.966114|3.55|1.473384|18/36|
|StressResidualNMO|Holdout_ThinInterbed|0.656751|0.618672|5.80|0.658766|0/36|
|StressResidualNMO|SoftLowQ|0.934197|0.932506|0.18|1.335472|18/36|
|StressResidualNMO|ThinHighFrequency|0.744920|0.715614|3.93|0.744445|0/36|
|StressResidualNMO|VelocityReversal|0.694868|0.683569|1.63|0.708556|0/36|
|StressVelocityError|Dev_DeepContrast|0.920253|0.919038|0.13|0.920006|0/36|
|StressVelocityError|Dev_ModerateQ|12.236479|13.610815|-11.23|20.535539|36/36|
|StressVelocityError|Holdout_DeepLowQ|0.847196|0.844017|0.38|0.848155|0/36|
|StressVelocityError|Holdout_HighContrastReversal|0.850872|0.848211|0.31|0.864774|0/36|
|StressVelocityError|Holdout_StrongAtten|14.457904|17.349282|-20.00|35.638661|36/36|
|StressVelocityError|Holdout_ThinInterbed|0.710914|0.683252|3.89|0.707203|0/36|
|StressVelocityError|SoftLowQ|14.120586|16.951143|-20.05|33.769727|36/36|
|StressVelocityError|ThinHighFrequency|0.806192|0.788913|2.14|0.803091|0/36|
|StressVelocityError|VelocityReversal|0.850682|0.849092|0.19|0.851798|0/36|

未通过项目：
- LowerSNR/Holdout_ThinInterbed/WeakEventError：{"baseline": 1.0292816908792286, "candidate": 1.1859410643091761, "allowance": 0.05146408454396143, "n": 18}
- Main/ThinHighFrequency/WeakEventError：{"baseline": 1.1419495008520673, "candidate": 1.2029134431845157, "allowance": 0.05709747504260337, "n": 18}
- StressVelocityError/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Dev_ModerateQ", "n": 36, "baseline": 12.236478804956636, "candidate": 13.610815139550665, "improvement_percent": -11.23146908927204, "baseline_median": 13.1079701613949, "candidate_median": 14.55674116655235, "baseline_p90": 15.78237672772405, "candidate_p90": 17.447085658771798, "baseline_p95": 16.010052777633923, "candidate_p95": 18.435741693050723, "baseline_worst": 17.3511560774854, "candidate_worst": 20.5355390430607, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Dev_ModerateQ/EventAmplitudeError：{"baseline": 3.237604398825965, "candidate": 3.5279353295057443, "allowance": 0.16188021994129825, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOInterceptError：{"baseline": 2.7282603219423742, "candidate": 3.019908298974721, "allowance": 0.13641301609711873, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOGradientError：{"baseline": 2.3894831285401654, "candidate": 2.741501644479322, "allowance": 0.11947415642700827, "n": 36}
- StressVelocityError/Holdout_StrongAtten/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Holdout_StrongAtten", "n": 36, "baseline": 14.457904392050219, "candidate": 17.3492818040553, "improvement_percent": -19.99859269781119, "baseline_median": 12.676617824936663, "candidate_median": 15.065821952004265, "baseline_p90": 23.9993398175421, "candidate_p90": 29.0846208120011, "baseline_p95": 26.6926891176284, "candidate_p95": 35.630902770260825, "baseline_worst": 26.6926892998688, "candidate_worst": 35.638661424219, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 3.7553818613177787, "candidate": 4.388511250547733, "allowance": 0.18776909306588896, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOInterceptError：{"baseline": 3.037779074454735, "candidate": 3.5490201725994655, "allowance": 0.15188895372273675, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOGradientError：{"baseline": 2.281710106016304, "candidate": 2.831644716715502, "allowance": 0.11408550530081521, "n": 36}
- StressVelocityError/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "SoftLowQ", "n": 36, "baseline": 14.120586199193193, "candidate": 16.951143307998215, "improvement_percent": -20.045606243788594, "baseline_median": 13.015832476676675, "candidate_median": 15.22478624595885, "baseline_p90": 21.57603703907565, "candidate_p90": 26.0822502457303, "baseline_p95": 25.36939763507525, "candidate_p95": 33.704547072309396, "baseline_worst": 25.4147687872566, "candidate_worst": 33.7697272547145, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/SoftLowQ/EventAmplitudeError：{"baseline": 4.558573031099418, "candidate": 5.383188687466484, "allowance": 0.2279286515549709, "n": 36}
- StressVelocityError/SoftLowQ/AVOInterceptError：{"baseline": 3.354832433157156, "candidate": 3.998836570234744, "allowance": 0.1677416216578578, "n": 36}
- StressVelocityError/SoftLowQ/AVOGradientError：{"baseline": 2.1350922144952227, "candidate": 2.7470157075430586, "allowance": 0.10675461072476114, "n": 36}
- main_equal_model_improvement：{"baseline": 0.6565302490183482, "candidate": 0.6380065377555891, "relative_improvement": 0.02821455872057066}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
