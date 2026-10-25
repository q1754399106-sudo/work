# NoiseUncertaintyOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 2916 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 4.01%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.677322|0.669066|1.22|0.711319|0/90|
|LowerSNR|Dev_ModerateQ|0.724330|0.717918|0.89|0.771192|0/90|
|LowerSNR|Holdout_DeepLowQ|0.667303|0.663289|0.60|0.804169|0/90|
|LowerSNR|Holdout_HighContrastReversal|0.694350|0.688564|0.83|0.712371|0/90|
|LowerSNR|Holdout_StrongAtten|0.853844|0.801973|6.07|1.107782|3/90|
|LowerSNR|Holdout_ThinInterbed|0.737769|0.729040|1.18|0.757206|0/90|
|LowerSNR|SoftLowQ|0.821895|0.782014|4.85|1.057517|3/90|
|LowerSNR|ThinHighFrequency|0.810706|0.804980|0.71|0.832895|0/90|
|LowerSNR|VelocityReversal|0.661767|0.656578|0.78|0.686691|0/90|
|Main|Dev_DeepContrast|0.565929|0.555478|1.85|0.596077|0/90|
|Main|Dev_ModerateQ|0.635546|0.626481|1.43|0.673956|0/90|
|Main|Holdout_DeepLowQ|0.610412|0.595316|2.47|0.760936|0/90|
|Main|Holdout_HighContrastReversal|0.594396|0.586165|1.38|0.606302|0/90|
|Main|Holdout_StrongAtten|0.803759|0.711030|11.54|0.935081|0/90|
|Main|Holdout_ThinInterbed|0.634161|0.622836|1.79|0.652782|0/90|
|Main|SoftLowQ|0.763472|0.688521|9.82|0.849389|0/90|
|Main|ThinHighFrequency|0.740526|0.733409|0.96|0.771869|0/90|
|Main|VelocityReversal|0.560572|0.552422|1.45|0.575800|0/90|
|StressColoredNoise|Dev_DeepContrast|0.599550|0.581998|2.93|0.605096|0/36|
|StressColoredNoise|Dev_ModerateQ|0.683567|0.668615|2.19|0.691057|0/36|
|StressColoredNoise|Holdout_DeepLowQ|0.626961|0.595864|4.96|0.608379|0/36|
|StressColoredNoise|Holdout_HighContrastReversal|0.630106|0.615733|2.28|0.629014|0/36|
|StressColoredNoise|Holdout_StrongAtten|0.965173|0.791536|17.99|0.932711|0/36|
|StressColoredNoise|Holdout_ThinInterbed|0.660187|0.643115|2.59|0.664813|0/36|
|StressColoredNoise|SoftLowQ|0.904496|0.760405|15.93|0.847731|0/36|
|StressColoredNoise|ThinHighFrequency|0.753189|0.743798|1.25|0.760641|0/36|
|StressColoredNoise|VelocityReversal|0.597856|0.582898|2.50|0.593254|0/36|
|StressNonRickerPhase|Dev_DeepContrast|0.697046|0.680469|2.38|0.698947|0/36|
|StressNonRickerPhase|Dev_ModerateQ|0.777500|0.761081|2.11|0.787412|0/36|
|StressNonRickerPhase|Holdout_DeepLowQ|0.738732|0.711440|3.69|0.721295|0/36|
|StressNonRickerPhase|Holdout_HighContrastReversal|0.722117|0.709622|1.73|0.717118|0/36|
|StressNonRickerPhase|Holdout_StrongAtten|1.076898|0.863516|19.81|1.011060|2/36|
|StressNonRickerPhase|Holdout_ThinInterbed|0.738735|0.721793|2.29|0.741774|0/36|
|StressNonRickerPhase|SoftLowQ|1.008869|0.842790|16.46|0.943207|0/36|
|StressNonRickerPhase|ThinHighFrequency|0.809304|0.798937|1.28|0.813353|0/36|
|StressNonRickerPhase|VelocityReversal|0.701777|0.689051|1.81|0.695065|0/36|
|StressResidualNMO|Dev_DeepContrast|0.627995|0.606418|3.44|0.628949|0/36|
|StressResidualNMO|Dev_ModerateQ|0.690217|0.670982|2.79|0.707144|0/36|
|StressResidualNMO|Holdout_DeepLowQ|0.650529|0.615416|5.40|0.629541|0/36|
|StressResidualNMO|Holdout_HighContrastReversal|0.722035|0.704031|2.49|0.718811|0/36|
|StressResidualNMO|Holdout_StrongAtten|1.001713|0.780708|22.06|0.941472|0/36|
|StressResidualNMO|Holdout_ThinInterbed|0.656751|0.635305|3.27|0.662620|0/36|
|StressResidualNMO|SoftLowQ|0.934197|0.759128|18.74|0.877534|0/36|
|StressResidualNMO|ThinHighFrequency|0.744920|0.731839|1.76|0.753104|0/36|
|StressResidualNMO|VelocityReversal|0.694868|0.673800|3.03|0.687292|0/36|
|StressVelocityError|Dev_DeepContrast|0.920253|0.920800|-0.06|0.922430|0/36|
|StressVelocityError|Dev_ModerateQ|12.236479|12.364954|-1.05|17.749838|36/36|
|StressVelocityError|Holdout_DeepLowQ|0.847196|0.845681|0.18|0.849239|0/36|
|StressVelocityError|Holdout_HighContrastReversal|0.850872|0.853046|-0.26|0.868476|0/36|
|StressVelocityError|Holdout_StrongAtten|14.457904|14.669856|-1.47|27.943938|36/36|
|StressVelocityError|Holdout_ThinInterbed|0.710914|0.700782|1.43|0.718786|0/36|
|StressVelocityError|SoftLowQ|14.120586|14.535037|-2.94|27.721657|36/36|
|StressVelocityError|ThinHighFrequency|0.806192|0.800525|0.70|0.811484|0/36|
|StressVelocityError|VelocityReversal|0.850682|0.852918|-0.26|0.855334|0/36|

未通过项目：
- LowerSNR/Holdout_ThinInterbed/WeakEventError：{"baseline": 1.0292816908792286, "candidate": 1.0814510795024395, "allowance": 0.05146408454396143, "n": 18}
- StressVelocityError/Dev_ModerateQ/AVOGradientError：{"baseline": 2.3894831285401654, "candidate": 2.635544111777002, "allowance": 0.11947415642700827, "n": 36}
- StressVelocityError/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 3.7553818613177787, "candidate": 4.087547364593786, "allowance": 0.18776909306588896, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOInterceptError：{"baseline": 3.037779074454735, "candidate": 3.3111294512025595, "allowance": 0.15188895372273675, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOGradientError：{"baseline": 2.281710106016304, "candidate": 2.457991448130243, "allowance": 0.11408550530081521, "n": 36}
- StressVelocityError/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "SoftLowQ", "n": 36, "baseline": 14.120586199193193, "candidate": 14.53503660886173, "improvement_percent": -2.935079350262515, "baseline_median": 13.015832476676675, "candidate_median": 13.032037108829691, "baseline_p90": 21.57603703907565, "candidate_p90": 21.910767312457452, "baseline_p95": 25.36939763507525, "candidate_p95": 27.660355924446776, "baseline_worst": 25.4147687872566, "candidate_worst": 27.7216565392708, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/SoftLowQ/EventAmplitudeError：{"baseline": 4.558573031099418, "candidate": 5.076247326632747, "allowance": 0.2279286515549709, "n": 36}
- StressVelocityError/SoftLowQ/AVOInterceptError：{"baseline": 3.354832433157156, "candidate": 3.73012329739852, "allowance": 0.1677416216578578, "n": 36}
- StressVelocityError/SoftLowQ/AVOGradientError：{"baseline": 2.1350922144952227, "candidate": 2.3035990855472424, "allowance": 0.10675461072476114, "n": 36}
- main_equal_model_improvement：{"baseline": 0.6565302490183482, "candidate": 0.6301841314720333, "relative_improvement": 0.04012932775863398}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
