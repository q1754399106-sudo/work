# CombinedUncertainty 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 2916 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 5.65%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.677322|0.643950|4.93|0.698597|0/90|
|LowerSNR|Dev_ModerateQ|0.724330|0.695845|3.93|0.747548|0/90|
|LowerSNR|Holdout_DeepLowQ|0.667303|0.656589|1.61|0.798299|0/90|
|LowerSNR|Holdout_HighContrastReversal|0.694350|0.671859|3.24|0.701291|0/90|
|LowerSNR|Holdout_StrongAtten|0.853844|0.796797|6.68|1.200742|3/90|
|LowerSNR|Holdout_ThinInterbed|0.737769|0.705339|4.40|0.734378|0/90|
|LowerSNR|SoftLowQ|0.821895|0.782871|4.75|1.194673|6/90|
|LowerSNR|ThinHighFrequency|0.810706|0.788900|2.69|0.824591|0/90|
|LowerSNR|VelocityReversal|0.661767|0.640761|3.17|0.675821|0/90|
|Main|Dev_DeepContrast|0.565929|0.534362|5.58|0.582891|0/90|
|Main|Dev_ModerateQ|0.635546|0.608149|4.31|0.657916|0/90|
|Main|Holdout_DeepLowQ|0.610412|0.589918|3.36|0.756212|0/90|
|Main|Holdout_HighContrastReversal|0.594396|0.570955|3.94|0.583834|0/90|
|Main|Holdout_StrongAtten|0.803759|0.715341|11.00|0.994405|0/90|
|Main|Holdout_ThinInterbed|0.634161|0.604026|4.75|0.644330|0/90|
|Main|SoftLowQ|0.763472|0.694792|9.00|0.921299|0/90|
|Main|ThinHighFrequency|0.740526|0.719379|2.86|0.766407|0/90|
|Main|VelocityReversal|0.560572|0.538178|3.99|0.560340|0/90|
|StressColoredNoise|Dev_DeepContrast|0.599550|0.559487|6.68|0.584753|0/36|
|StressColoredNoise|Dev_ModerateQ|0.683567|0.651003|4.76|0.674608|0/36|
|StressColoredNoise|Holdout_DeepLowQ|0.626961|0.589057|6.05|0.601754|0/36|
|StressColoredNoise|Holdout_HighContrastReversal|0.630106|0.597500|5.17|0.610559|0/36|
|StressColoredNoise|Holdout_StrongAtten|0.965173|0.803733|16.73|0.966826|0/36|
|StressColoredNoise|Holdout_ThinInterbed|0.660187|0.621786|5.82|0.643575|0/36|
|StressColoredNoise|SoftLowQ|0.904496|0.778773|13.90|0.898170|0/36|
|StressColoredNoise|ThinHighFrequency|0.753189|0.729589|3.13|0.747983|0/36|
|StressColoredNoise|VelocityReversal|0.597856|0.564712|5.54|0.574424|0/36|
|StressNonRickerPhase|Dev_DeepContrast|0.697046|0.662605|4.94|0.682263|0/36|
|StressNonRickerPhase|Dev_ModerateQ|0.777500|0.746159|4.03|0.774252|0/36|
|StressNonRickerPhase|Holdout_DeepLowQ|0.738732|0.712852|3.50|0.724977|0/36|
|StressNonRickerPhase|Holdout_HighContrastReversal|0.722117|0.697585|3.40|0.706310|0/36|
|StressNonRickerPhase|Holdout_StrongAtten|1.076898|0.885756|17.75|1.074135|15/36|
|StressNonRickerPhase|Holdout_ThinInterbed|0.738735|0.703193|4.81|0.724025|0/36|
|StressNonRickerPhase|SoftLowQ|1.008869|0.864950|14.27|1.004164|1/36|
|StressNonRickerPhase|ThinHighFrequency|0.809304|0.782802|3.27|0.799370|0/36|
|StressNonRickerPhase|VelocityReversal|0.701777|0.678320|3.34|0.685926|0/36|
|StressResidualNMO|Dev_DeepContrast|0.627995|0.583757|7.04|0.606435|0/36|
|StressResidualNMO|Dev_ModerateQ|0.690217|0.648733|6.01|0.689550|0/36|
|StressResidualNMO|Holdout_DeepLowQ|0.650529|0.612039|5.92|0.627434|0/36|
|StressResidualNMO|Holdout_HighContrastReversal|0.722035|0.696779|3.50|0.714878|0/36|
|StressResidualNMO|Holdout_StrongAtten|1.001713|0.799634|20.17|1.000561|1/36|
|StressResidualNMO|Holdout_ThinInterbed|0.656751|0.609318|7.22|0.638103|0/36|
|StressResidualNMO|SoftLowQ|0.934197|0.781581|16.34|0.948246|0/36|
|StressResidualNMO|ThinHighFrequency|0.744920|0.710315|4.65|0.734230|0/36|
|StressResidualNMO|VelocityReversal|0.694868|0.666304|4.11|0.681555|0/36|
|StressVelocityError|Dev_DeepContrast|0.920253|0.918739|0.16|0.920197|0/36|
|StressVelocityError|Dev_ModerateQ|12.236479|12.803399|-4.63|18.905602|36/36|
|StressVelocityError|Holdout_DeepLowQ|0.847196|0.842679|0.53|0.846145|0/36|
|StressVelocityError|Holdout_HighContrastReversal|0.850872|0.847850|0.36|0.864209|0/36|
|StressVelocityError|Holdout_StrongAtten|14.457904|15.913836|-10.07|32.709105|36/36|
|StressVelocityError|Holdout_ThinInterbed|0.710914|0.679074|4.48|0.698193|0/36|
|StressVelocityError|SoftLowQ|14.120586|15.643921|-10.79|30.942101|36/36|
|StressVelocityError|ThinHighFrequency|0.806192|0.787089|2.37|0.800029|0/36|
|StressVelocityError|VelocityReversal|0.850682|0.848609|0.24|0.850883|0/36|

未通过项目：
- LowerSNR/Holdout_ThinInterbed/WeakEventError：{"baseline": 1.0292816908792286, "candidate": 1.1859410643091761, "allowance": 0.05146408454396143, "n": 18}
- Main/ThinHighFrequency/WeakEventError：{"baseline": 1.1419495008520673, "candidate": 1.2029134431845157, "allowance": 0.05709747504260337, "n": 18}
- StressVelocityError/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Dev_ModerateQ", "n": 36, "baseline": 12.236478804956636, "candidate": 12.803399145693646, "improvement_percent": -4.633034958613802, "baseline_median": 13.1079701613949, "candidate_median": 13.457355750030949, "baseline_p90": 15.78237672772405, "candidate_p90": 16.6376961284395, "baseline_p95": 16.010052777633923, "candidate_p95": 17.4024241709052, "baseline_worst": 17.3511560774854, "candidate_worst": 18.9056018216322, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Dev_ModerateQ/EventAmplitudeError：{"baseline": 3.237604398825965, "candidate": 3.4651123957071315, "allowance": 0.16188021994129825, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOInterceptError：{"baseline": 2.7282603219423742, "candidate": 2.982515134787779, "allowance": 0.13641301609711873, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOGradientError：{"baseline": 2.3894831285401654, "candidate": 2.8017961793107125, "allowance": 0.11947415642700827, "n": 36}
- StressVelocityError/Holdout_StrongAtten/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Holdout_StrongAtten", "n": 36, "baseline": 14.457904392050219, "candidate": 15.913835640285011, "improvement_percent": -10.07014024131565, "baseline_median": 12.676617824936663, "candidate_median": 13.48624467105379, "baseline_p90": 23.9993398175421, "candidate_p90": 27.48626766753935, "baseline_p95": 26.6926891176284, "candidate_p95": 32.70219546471733, "baseline_worst": 26.6926892998688, "candidate_worst": 32.7091049732802, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 3.7553818613177787, "candidate": 4.243865303375985, "allowance": 0.18776909306588896, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOInterceptError：{"baseline": 3.037779074454735, "candidate": 3.4748344437541947, "allowance": 0.15188895372273675, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOGradientError：{"baseline": 2.281710106016304, "candidate": 2.817732306279983, "allowance": 0.11408550530081521, "n": 36}
- StressVelocityError/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "SoftLowQ", "n": 36, "baseline": 14.120586199193193, "candidate": 15.643920657594458, "improvement_percent": -10.788039794610672, "baseline_median": 13.015832476676675, "candidate_median": 14.163460797908051, "baseline_p90": 21.57603703907565, "candidate_p90": 24.47655443576945, "baseline_p95": 25.36939763507525, "candidate_p95": 30.8736716429768, "baseline_worst": 25.4147687872566, "candidate_worst": 30.9421007362676, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/SoftLowQ/EventAmplitudeError：{"baseline": 4.558573031099418, "candidate": 5.2378895151726095, "allowance": 0.2279286515549709, "n": 36}
- StressVelocityError/SoftLowQ/AVOInterceptError：{"baseline": 3.354832433157156, "candidate": 3.9260610167413685, "allowance": 0.1677416216578578, "n": 36}
- StressVelocityError/SoftLowQ/AVOGradientError：{"baseline": 2.1350922144952227, "candidate": 2.7646271897218186, "allowance": 0.10675461072476114, "n": 36}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
