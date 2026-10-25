# StrongerDamping 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 2916 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 -4.04%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.677322|0.725312|-7.09|0.762217|0/90|
|LowerSNR|Dev_ModerateQ|0.724330|0.769198|-6.19|0.807695|0/90|
|LowerSNR|Holdout_DeepLowQ|0.667303|0.690200|-3.43|0.826597|0/90|
|LowerSNR|Holdout_HighContrastReversal|0.694350|0.728432|-4.91|0.754450|0/90|
|LowerSNR|Holdout_StrongAtten|0.853844|0.859396|-0.65|1.134369|9/90|
|LowerSNR|Holdout_ThinInterbed|0.737769|0.781897|-5.98|0.807380|0/90|
|LowerSNR|SoftLowQ|0.821895|0.836134|-1.73|1.025390|3/90|
|LowerSNR|ThinHighFrequency|0.810706|0.841026|-3.74|0.865629|0/90|
|LowerSNR|VelocityReversal|0.661767|0.694347|-4.92|0.733154|0/90|
|Main|Dev_DeepContrast|0.565929|0.606887|-7.24|0.638785|0/90|
|Main|Dev_ModerateQ|0.635546|0.672949|-5.89|0.726930|0/90|
|Main|Holdout_DeepLowQ|0.610412|0.624375|-2.29|0.779420|0/90|
|Main|Holdout_HighContrastReversal|0.594396|0.624023|-4.98|0.658149|0/90|
|Main|Holdout_StrongAtten|0.803759|0.811327|-0.94|1.322755|18/90|
|Main|Holdout_ThinInterbed|0.634161|0.673555|-6.21|0.702705|0/90|
|Main|SoftLowQ|0.763472|0.777525|-1.84|1.171772|18/90|
|Main|ThinHighFrequency|0.740526|0.769004|-3.85|0.800239|0/90|
|Main|VelocityReversal|0.560572|0.587753|-4.85|0.630825|0/90|
|StressColoredNoise|Dev_DeepContrast|0.599550|0.627858|-4.72|0.654622|0/36|
|StressColoredNoise|Dev_ModerateQ|0.683567|0.707237|-3.46|0.732898|0/36|
|StressColoredNoise|Holdout_DeepLowQ|0.626961|0.636595|-1.54|0.677708|0/36|
|StressColoredNoise|Holdout_HighContrastReversal|0.630106|0.651950|-3.47|0.672860|0/36|
|StressColoredNoise|Holdout_StrongAtten|0.965173|0.940621|2.54|1.291914|18/36|
|StressColoredNoise|Holdout_ThinInterbed|0.660187|0.687312|-4.11|0.712665|0/36|
|StressColoredNoise|SoftLowQ|0.904496|0.892289|1.35|1.132793|15/36|
|StressColoredNoise|ThinHighFrequency|0.753189|0.770917|-2.35|0.787578|0/36|
|StressColoredNoise|VelocityReversal|0.597856|0.620044|-3.71|0.638280|0/36|
|StressNonRickerPhase|Dev_DeepContrast|0.697046|0.718009|-3.01|0.742641|0/36|
|StressNonRickerPhase|Dev_ModerateQ|0.777500|0.795570|-2.32|0.826267|0/36|
|StressNonRickerPhase|Holdout_DeepLowQ|0.738732|0.739220|-0.07|0.775536|0/36|
|StressNonRickerPhase|Holdout_HighContrastReversal|0.722117|0.736570|-2.00|0.750803|0/36|
|StressNonRickerPhase|Holdout_StrongAtten|1.076898|1.023143|4.99|1.380471|18/36|
|StressNonRickerPhase|Holdout_ThinInterbed|0.738735|0.761785|-3.12|0.787132|0/36|
|StressNonRickerPhase|SoftLowQ|1.008869|0.975093|3.35|1.243652|18/36|
|StressNonRickerPhase|ThinHighFrequency|0.809304|0.826523|-2.13|0.843542|0/36|
|StressNonRickerPhase|VelocityReversal|0.701777|0.715255|-1.92|0.731039|0/36|
|StressResidualNMO|Dev_DeepContrast|0.627995|0.654234|-4.18|0.683437|0/36|
|StressResidualNMO|Dev_ModerateQ|0.690217|0.716450|-3.80|0.755814|0/36|
|StressResidualNMO|Holdout_DeepLowQ|0.650529|0.656818|-0.97|0.703805|0/36|
|StressResidualNMO|Holdout_HighContrastReversal|0.722035|0.731836|-1.36|0.749641|0/36|
|StressResidualNMO|Holdout_StrongAtten|1.001713|0.954528|4.71|1.328207|18/36|
|StressResidualNMO|Holdout_ThinInterbed|0.656751|0.688128|-4.78|0.720420|0/36|
|StressResidualNMO|SoftLowQ|0.934197|0.906267|2.99|1.191936|18/36|
|StressResidualNMO|ThinHighFrequency|0.744920|0.767730|-3.06|0.789392|0/36|
|StressResidualNMO|VelocityReversal|0.694868|0.704166|-1.34|0.726792|0/36|
|StressVelocityError|Dev_DeepContrast|0.920253|0.922609|-0.26|0.923824|0/36|
|StressVelocityError|Dev_ModerateQ|12.236479|9.197156|24.84|12.568385|36/36|
|StressVelocityError|Holdout_DeepLowQ|0.847196|0.850806|-0.43|0.855080|0/36|
|StressVelocityError|Holdout_HighContrastReversal|0.850872|0.856717|-0.69|0.871724|0/36|
|StressVelocityError|Holdout_StrongAtten|14.457904|9.236201|36.12|14.956092|36/36|
|StressVelocityError|Holdout_ThinInterbed|0.710914|0.735467|-3.45|0.753767|0/36|
|StressVelocityError|SoftLowQ|14.120586|9.217707|34.72|14.496720|36/36|
|StressVelocityError|ThinHighFrequency|0.806192|0.820794|-1.81|0.832167|0/36|
|StressVelocityError|VelocityReversal|0.850682|0.855578|-0.58|0.857854|0/36|

未通过项目：
- LowerSNR/Dev_DeepContrast/mean_noninferiority：{"Suite": "LowerSNR", "Model": "Dev_DeepContrast", "n": 90, "baseline": 0.6773218340343108, "candidate": 0.7253116770550557, "improvement_percent": -7.085234907444882, "baseline_median": 0.6760093211513235, "candidate_median": 0.7238152546162204, "baseline_p90": 0.7026608957993608, "candidate_p90": 0.7506874689861031, "baseline_p95": 0.7084506599265681, "candidate_p95": 0.7579415039463897, "baseline_worst": 0.712593089456478, "candidate_worst": 0.76221669536163, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/Dev_DeepContrast/EventAmplitudeError：{"baseline": 0.5830717980545005, "candidate": 0.6523102108655262, "allowance": 0.029153589902725027, "n": 90}
- LowerSNR/Dev_DeepContrast/WeakEventError：{"baseline": 0.505421067966303, "candidate": 0.5750908329465325, "allowance": 0.025271053398315153, "n": 18}
- LowerSNR/Dev_DeepContrast/FaultNeighborhoodError：{"baseline": 0.6687898542252189, "candidate": 0.7234642527548737, "allowance": 0.03343949271126095, "n": 18}
- LowerSNR/Dev_DeepContrast/AVOInterceptError：{"baseline": 0.5735436651050138, "candidate": 0.6319327303806477, "allowance": 0.02867718325525069, "n": 36}
- LowerSNR/Dev_DeepContrast/AVOGradientError：{"baseline": 0.6001290712998716, "candidate": 0.6561311453026544, "allowance": 0.030006453564993582, "n": 36}
- LowerSNR/Dev_ModerateQ/mean_noninferiority：{"Suite": "LowerSNR", "Model": "Dev_ModerateQ", "n": 90, "baseline": 0.7243299149267702, "candidate": 0.7691975271250115, "improvement_percent": -6.1943613364053185, "baseline_median": 0.7215823618384665, "candidate_median": 0.7696169901754355, "baseline_p90": 0.7521609096766868, "candidate_p90": 0.7863805019621015, "baseline_p95": 0.7604186663309265, "candidate_p95": 0.7945035549112721, "baseline_worst": 0.778485768393536, "candidate_worst": 0.807695324203897, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/Dev_ModerateQ/EventAmplitudeError：{"baseline": 0.6595179124913435, "candidate": 0.7187186304717041, "allowance": 0.032975895624567174, "n": 90}
- LowerSNR/Dev_ModerateQ/WeakEventError：{"baseline": 0.5822509659059155, "candidate": 0.6315996238565342, "allowance": 0.029112548295295776, "n": 18}
- LowerSNR/Dev_ModerateQ/FaultNeighborhoodError：{"baseline": 0.7076547710582758, "candidate": 0.7607840866377351, "allowance": 0.03538273855291379, "n": 18}
- LowerSNR/Dev_ModerateQ/AVOInterceptError：{"baseline": 0.666995618134488, "candidate": 0.7164806905295601, "allowance": 0.033349780906724406, "n": 36}
- LowerSNR/Dev_ModerateQ/AVOGradientError：{"baseline": 0.7273412793984401, "candidate": 0.7746652649424473, "allowance": 0.036367063969922, "n": 36}
- LowerSNR/Holdout_DeepLowQ/mean_noninferiority：{"Suite": "LowerSNR", "Model": "Holdout_DeepLowQ", "n": 90, "baseline": 0.6673031034536241, "candidate": 0.6901995267150473, "improvement_percent": -3.4311878879212276, "baseline_median": 0.6912110095207755, "candidate_median": 0.7067779634204595, "baseline_p90": 0.8012854485024485, "candidate_p90": 0.8223000338530514, "baseline_p95": 0.8035810982744676, "candidate_p95": 0.8238920254761575, "baseline_worst": 0.806260734350437, "candidate_worst": 0.826597219001289, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/Holdout_DeepLowQ/EventAmplitudeError：{"baseline": 0.5685606642017438, "candidate": 0.6101880709565126, "allowance": 0.02842803321008719, "n": 90}
- LowerSNR/Holdout_DeepLowQ/WeakEventError：{"baseline": 0.33233418276523063, "candidate": 0.43485019284625337, "allowance": 0.01661670913826153, "n": 18}
- LowerSNR/Holdout_DeepLowQ/AVOInterceptError：{"baseline": 0.5112130145272216, "candidate": 0.5413438275640444, "allowance": 0.02556065072636108, "n": 36}
- LowerSNR/Holdout_DeepLowQ/AVOGradientError：{"baseline": 0.5817415278844313, "candidate": 0.6175179495912121, "allowance": 0.02908707639422157, "n": 36}
- LowerSNR/Holdout_HighContrastReversal/mean_noninferiority：{"Suite": "LowerSNR", "Model": "Holdout_HighContrastReversal", "n": 90, "baseline": 0.6943498499312114, "candidate": 0.7284323594884613, "improvement_percent": -4.908549999776968, "baseline_median": 0.70179579953033, "candidate_median": 0.7377222504275415, "baseline_p90": 0.7154669116895616, "candidate_p90": 0.7513705976220115, "baseline_p95": 0.7172120208876537, "candidate_p95": 0.7526086931190912, "baseline_worst": 0.719898801268482, "candidate_worst": 0.754449921146527, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/Holdout_HighContrastReversal/EventAmplitudeError：{"baseline": 0.5947325180564691, "candidate": 0.650938880486641, "allowance": 0.029736625902823458, "n": 90}
- LowerSNR/Holdout_HighContrastReversal/FaultNeighborhoodError：{"baseline": 0.7017306852055745, "candidate": 0.7453385537942689, "allowance": 0.035086534260278725, "n": 18}
- LowerSNR/Holdout_HighContrastReversal/AVOInterceptError：{"baseline": 0.6053989181571502, "candidate": 0.6419001223286895, "allowance": 0.03026994590785751, "n": 36}
- LowerSNR/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 0.664174832811131, "candidate": 0.7409215732155001, "allowance": 0.03320874164055655, "n": 90}
- LowerSNR/Holdout_StrongAtten/WeakEventError：{"baseline": 0.6207810888154383, "candidate": 0.6675821624674622, "allowance": 0.031039054440771915, "n": 18}
- LowerSNR/Holdout_StrongAtten/FaultNeighborhoodError：{"baseline": 0.7659726889137061, "candidate": 0.8200944637351892, "allowance": 0.038298634445685306, "n": 18}
- LowerSNR/Holdout_StrongAtten/AVOInterceptError：{"baseline": 0.6489023637511987, "candidate": 0.7097183029897305, "allowance": 0.032445118187559935, "n": 36}
- LowerSNR/Holdout_StrongAtten/AVOGradientError：{"baseline": 0.7587911522402302, "candidate": 0.8107448440430659, "allowance": 0.03793955761201151, "n": 36}
- LowerSNR/Holdout_ThinInterbed/mean_noninferiority：{"Suite": "LowerSNR", "Model": "Holdout_ThinInterbed", "n": 90, "baseline": 0.7377687256006437, "candidate": 0.7818967236583437, "improvement_percent": -5.981277943406171, "baseline_median": 0.7390499825204024, "candidate_median": 0.78309973079607, "baseline_p90": 0.7516859557533527, "candidate_p90": 0.7994978925311161, "baseline_p95": 0.7526029990569079, "candidate_p95": 0.8036955150172649, "baseline_worst": 0.767110660354087, "candidate_worst": 0.807379836931177, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/Holdout_ThinInterbed/EventAmplitudeError：{"baseline": 0.6700815330140849, "candidate": 0.7281917358804382, "allowance": 0.03350407665070425, "n": 90}
- LowerSNR/Holdout_ThinInterbed/FaultNeighborhoodError：{"baseline": 0.7399552037385844, "candidate": 0.7975002574248602, "allowance": 0.03699776018692922, "n": 18}
- LowerSNR/Holdout_ThinInterbed/AVOInterceptError：{"baseline": 0.6566925243124815, "candidate": 0.7038294643541969, "allowance": 0.03283462621562407, "n": 36}
- LowerSNR/Holdout_ThinInterbed/AVOGradientError：{"baseline": 0.6870368894847895, "candidate": 0.7317718063445398, "allowance": 0.034351844474239474, "n": 36}
- LowerSNR/SoftLowQ/EventAmplitudeError：{"baseline": 0.6531756909326478, "candidate": 0.729458496547411, "allowance": 0.032658784546632395, "n": 90}
- LowerSNR/SoftLowQ/WeakEventError：{"baseline": 0.6428976449629661, "candidate": 0.6804523795878891, "allowance": 0.03214488224814831, "n": 18}
- LowerSNR/SoftLowQ/FaultNeighborhoodError：{"baseline": 0.7477894044780047, "candidate": 0.802052049053194, "allowance": 0.03738947022390023, "n": 18}
- LowerSNR/SoftLowQ/AVOInterceptError：{"baseline": 0.66097937002109, "candidate": 0.7196883762988314, "allowance": 0.0330489685010545, "n": 36}
- LowerSNR/SoftLowQ/AVOGradientError：{"baseline": 0.7491640969725849, "candidate": 0.8003225806206666, "allowance": 0.037458204848629245, "n": 36}
- LowerSNR/ThinHighFrequency/mean_noninferiority：{"Suite": "LowerSNR", "Model": "ThinHighFrequency", "n": 90, "baseline": 0.8107064572090406, "candidate": 0.841026269359856, "improvement_percent": -3.7399248373073624, "baseline_median": 0.813738128141657, "candidate_median": 0.8459217666021415, "baseline_p90": 0.8247776521654163, "candidate_p90": 0.8576815256552599, "baseline_p95": 0.8287747643022241, "candidate_p95": 0.8588794998144965, "baseline_worst": 0.836755799502568, "candidate_worst": 0.865628773137444, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/ThinHighFrequency/WeakEventError：{"baseline": 1.23523008516326, "candidate": 1.3860639594122417, "allowance": 0.06176150425816301, "n": 18}
- LowerSNR/VelocityReversal/mean_noninferiority：{"Suite": "LowerSNR", "Model": "VelocityReversal", "n": 90, "baseline": 0.661766900146344, "candidate": 0.6943468126372235, "improvement_percent": -4.923170452265713, "baseline_median": 0.6733797678797959, "candidate_median": 0.711870802667921, "baseline_p90": 0.6875102796038308, "candidate_p90": 0.7243451544379711, "baseline_p95": 0.6891892302083844, "candidate_p95": 0.7286550479999307, "baseline_worst": 0.694331992301704, "candidate_worst": 0.733154256764099, "baseline_failures": 0, "candidate_failures": 0}
- LowerSNR/VelocityReversal/EventAmplitudeError：{"baseline": 0.5586107521830113, "candidate": 0.6135646327618676, "allowance": 0.027930537609150565, "n": 90}
- LowerSNR/VelocityReversal/FaultNeighborhoodError：{"baseline": 0.6835934665422284, "candidate": 0.7271672250290074, "allowance": 0.03417967332711142, "n": 18}
- LowerSNR/VelocityReversal/AVOInterceptError：{"baseline": 0.5657088517258687, "candidate": 0.6029935591174831, "allowance": 0.02828544258629344, "n": 36}
- Main/Dev_DeepContrast/mean_noninferiority：{"Suite": "Main", "Model": "Dev_DeepContrast", "n": 90, "baseline": 0.5659290023635604, "candidate": 0.6068867198376546, "improvement_percent": -7.237253666632637, "baseline_median": 0.5660321734871454, "candidate_median": 0.6068306655046065, "baseline_p90": 0.5901045277525944, "candidate_p90": 0.6301496386966784, "baseline_p95": 0.5950770981912546, "candidate_p95": 0.6325673645507331, "baseline_worst": 0.605208248559707, "candidate_worst": 0.638785100259991, "baseline_failures": 0, "candidate_failures": 0}
- Main/Dev_DeepContrast/EventAmplitudeError：{"baseline": 0.44080392825913395, "candidate": 0.4979989469705819, "allowance": 0.022040196412956698, "n": 90}
- Main/Dev_DeepContrast/WeakEventError：{"baseline": 0.37567243112387844, "candidate": 0.4198181278375331, "allowance": 0.018783621556193925, "n": 18}
- Main/Dev_DeepContrast/FaultNeighborhoodError：{"baseline": 0.557682611845715, "candidate": 0.6021444524345253, "allowance": 0.02788413059228575, "n": 18}
- Main/Dev_DeepContrast/AVOInterceptError：{"baseline": 0.44630932950394914, "candidate": 0.4941061571220282, "allowance": 0.022315466475197457, "n": 36}
- Main/Dev_DeepContrast/AVOGradientError：{"baseline": 0.48221019127747694, "candidate": 0.5288713310622615, "allowance": 0.02411050956387385, "n": 36}
- Main/Dev_ModerateQ/mean_noninferiority：{"Suite": "Main", "Model": "Dev_ModerateQ", "n": 90, "baseline": 0.6355455024245654, "candidate": 0.6729486873864688, "improvement_percent": -5.885209606426711, "baseline_median": 0.632725437759406, "candidate_median": 0.671545205299608, "baseline_p90": 0.672558431227955, "candidate_p90": 0.7029954286675283, "baseline_p95": 0.6825379291195556, "candidate_p95": 0.709816326595398, "baseline_worst": 0.701959625326598, "candidate_worst": 0.726929794477053, "baseline_failures": 0, "candidate_failures": 0}
- Main/Dev_ModerateQ/EventAmplitudeError：{"baseline": 0.55506200417623, "candidate": 0.6016369372514204, "allowance": 0.0277531002088115, "n": 90}
- Main/Dev_ModerateQ/WeakEventError：{"baseline": 0.5217222669785073, "candidate": 0.5481395418724904, "allowance": 0.02608611334892537, "n": 18}
- Main/Dev_ModerateQ/FaultNeighborhoodError：{"baseline": 0.6141606918324858, "candidate": 0.6586499051795254, "allowance": 0.03070803459162429, "n": 18}
- Main/Dev_ModerateQ/AVOInterceptError：{"baseline": 0.5713120699614594, "candidate": 0.6095205774763688, "allowance": 0.028565603498072974, "n": 36}
- Main/Dev_ModerateQ/AVOGradientError：{"baseline": 0.6294720841621477, "candidate": 0.6698710770681049, "allowance": 0.03147360420810739, "n": 36}
- Main/Holdout_DeepLowQ/mean_noninferiority：{"Suite": "Main", "Model": "Holdout_DeepLowQ", "n": 90, "baseline": 0.6104121926560702, "candidate": 0.624374547943413, "improvement_percent": -2.2873650715574234, "baseline_median": 0.639423153037722, "candidate_median": 0.65755924916011, "baseline_p90": 0.7561506405465741, "candidate_p90": 0.7740937213807484, "baseline_p95": 0.76045112933405, "candidate_p95": 0.7770895306927644, "baseline_worst": 0.762579950948439, "candidate_worst": 0.779419631032414, "baseline_failures": 0, "candidate_failures": 0}
- Main/Holdout_DeepLowQ/WeakEventError：{"baseline": 0.16639957411731934, "candidate": 0.23044258774417734, "allowance": 0.008319978705865967, "n": 18}
- Main/Holdout_HighContrastReversal/mean_noninferiority：{"Suite": "Main", "Model": "Holdout_HighContrastReversal", "n": 90, "baseline": 0.5943960820764177, "candidate": 0.6240229845639781, "improvement_percent": -4.984370419142747, "baseline_median": 0.5918295697688305, "candidate_median": 0.6259097233305545, "baseline_p90": 0.6259072996816196, "candidate_p90": 0.6523990678154112, "baseline_p95": 0.6289010884996135, "candidate_p95": 0.654280276946713, "baseline_worst": 0.631492822001082, "candidate_worst": 0.658149263489411, "baseline_failures": 0, "candidate_failures": 0}
- Main/Holdout_HighContrastReversal/EventAmplitudeError：{"baseline": 0.4578174675134457, "candidate": 0.5056730428078096, "allowance": 0.022890873375672287, "n": 90}
- Main/Holdout_HighContrastReversal/FaultNeighborhoodError：{"baseline": 0.5830236499322279, "candidate": 0.6219125517704531, "allowance": 0.029151182496611396, "n": 18}
- Main/Holdout_HighContrastReversal/AVOInterceptError：{"baseline": 0.5082133636135573, "candidate": 0.5423492776159282, "allowance": 0.025410668180677866, "n": 36}
- Main/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 0.5389363319939608, "candidate": 0.611176762520759, "allowance": 0.026946816599698043, "n": 90}
- Main/Holdout_StrongAtten/WeakEventError：{"baseline": 0.45967325423196304, "candidate": 0.5080755651201063, "allowance": 0.022983662711598154, "n": 18}
- Main/Holdout_StrongAtten/FaultNeighborhoodError：{"baseline": 0.6769442916958105, "candidate": 0.7258016524454319, "allowance": 0.03384721458479053, "n": 18}
- Main/Holdout_StrongAtten/AVOInterceptError：{"baseline": 0.539769153831161, "candidate": 0.5942659497747741, "allowance": 0.02698845769155805, "n": 36}
- Main/Holdout_StrongAtten/AVOGradientError：{"baseline": 0.6580356372659678, "candidate": 0.710500541596928, "allowance": 0.03290178186329839, "n": 36}
- Main/Holdout_ThinInterbed/mean_noninferiority：{"Suite": "Main", "Model": "Holdout_ThinInterbed", "n": 90, "baseline": 0.6341612263689489, "candidate": 0.6735547731209793, "improvement_percent": -6.211913487298526, "baseline_median": 0.6359556021593185, "candidate_median": 0.6751535003603335, "baseline_p90": 0.6511575608003802, "candidate_p90": 0.689065018080407, "baseline_p95": 0.6568226811319722, "candidate_p95": 0.6908730241416838, "baseline_worst": 0.671607337019677, "candidate_worst": 0.702704638168964, "baseline_failures": 0, "candidate_failures": 0}
- Main/Holdout_ThinInterbed/EventAmplitudeError：{"baseline": 0.5474911078812102, "candidate": 0.5974933090715773, "allowance": 0.027374555394060515, "n": 90}
- Main/Holdout_ThinInterbed/FaultNeighborhoodError：{"baseline": 0.6199522865508265, "candidate": 0.6661511968414751, "allowance": 0.030997614327541324, "n": 18}
- Main/Holdout_ThinInterbed/AVOInterceptError：{"baseline": 0.545563575122299, "candidate": 0.5867864944416339, "allowance": 0.027278178756114952, "n": 36}
- Main/Holdout_ThinInterbed/AVOGradientError：{"baseline": 0.5748636068855963, "candidate": 0.6169045997862972, "allowance": 0.02874318034427982, "n": 36}
- Main/SoftLowQ/EventAmplitudeError：{"baseline": 0.5314970537855667, "candidate": 0.599669008265029, "allowance": 0.026574852689278335, "n": 90}
- Main/SoftLowQ/WeakEventError：{"baseline": 0.503443284430463, "candidate": 0.536533428037339, "allowance": 0.025172164221523154, "n": 18}
- Main/SoftLowQ/FaultNeighborhoodError：{"baseline": 0.6581669517333751, "candidate": 0.7042435817290771, "allowance": 0.03290834758666876, "n": 18}
- Main/SoftLowQ/AVOInterceptError：{"baseline": 0.5597252482007052, "candidate": 0.6095402423842199, "allowance": 0.027986262410035264, "n": 36}
- Main/SoftLowQ/AVOGradientError：{"baseline": 0.6535189089265551, "candidate": 0.7024239620912903, "allowance": 0.03267594544632776, "n": 36}
- Main/ThinHighFrequency/mean_noninferiority：{"Suite": "Main", "Model": "ThinHighFrequency", "n": 90, "baseline": 0.7405256035292337, "candidate": 0.7690041715211023, "improvement_percent": -3.845723612545471, "baseline_median": 0.7411572435643266, "candidate_median": 0.7709714683187325, "baseline_p90": 0.760220402549043, "candidate_p90": 0.7857961651552136, "baseline_p95": 0.7653512275718178, "candidate_p95": 0.7919318127504743, "baseline_worst": 0.775344957616933, "candidate_worst": 0.800239394135207, "baseline_failures": 0, "candidate_failures": 0}
- Main/ThinHighFrequency/FaultNeighborhoodError：{"baseline": 0.7352406525147257, "candidate": 0.7726229002411402, "allowance": 0.03676203262573629, "n": 18}
- Main/VelocityReversal/mean_noninferiority：{"Suite": "Main", "Model": "VelocityReversal", "n": 90, "baseline": 0.5605720190860246, "candidate": 0.5877530454806525, "improvement_percent": -4.848801843328676, "baseline_median": 0.5550234503195111, "candidate_median": 0.5899081305220445, "baseline_p90": 0.5925880557991208, "candidate_p90": 0.6201116966035246, "baseline_p95": 0.5972796904416281, "candidate_p95": 0.6235831236139988, "baseline_worst": 0.607793053235718, "candidate_worst": 0.630824562233564, "baseline_failures": 0, "candidate_failures": 0}
- Main/VelocityReversal/EventAmplitudeError：{"baseline": 0.425322338008518, "candidate": 0.46951592654098967, "allowance": 0.021266116900425902, "n": 90}
- Main/VelocityReversal/FaultNeighborhoodError：{"baseline": 0.5551339998944124, "candidate": 0.594045338965497, "allowance": 0.02775669999472062, "n": 18}
- Main/VelocityReversal/AVOInterceptError：{"baseline": 0.46821783880638385, "candidate": 0.5010437464496141, "allowance": 0.023410891940319194, "n": 36}
- Main/VelocityReversal/AVOGradientError：{"baseline": 0.69691627015373, "candidate": 0.7320226537647444, "allowance": 0.0348458135076865, "n": 36}
- StressColoredNoise/Dev_DeepContrast/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "Dev_DeepContrast", "n": 36, "baseline": 0.5995501818009573, "candidate": 0.6278584898792892, "improvement_percent": -4.7215911090708, "baseline_median": 0.6032937163598815, "candidate_median": 0.6302986452019685, "baseline_p90": 0.6191072416977146, "candidate_p90": 0.645204701988235, "baseline_p95": 0.622156945765993, "candidate_p95": 0.6485437870725643, "baseline_worst": 0.631280710736535, "candidate_worst": 0.654621587389403, "baseline_failures": 0, "candidate_failures": 0}
- StressColoredNoise/Dev_DeepContrast/EventAmplitudeError：{"baseline": 0.499648187482267, "candidate": 0.5371581328471781, "allowance": 0.024982409374113352, "n": 36}
- StressColoredNoise/Dev_DeepContrast/AVOInterceptError：{"baseline": 0.4810318582290118, "candidate": 0.5170389224997219, "allowance": 0.024051592911450592, "n": 36}
- StressColoredNoise/Dev_DeepContrast/AVOGradientError：{"baseline": 0.5058003384988475, "candidate": 0.5411078548616775, "allowance": 0.025290016924942378, "n": 36}
- StressColoredNoise/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "Dev_ModerateQ", "n": 36, "baseline": 0.6835671380829454, "candidate": 0.7072373221990814, "improvement_percent": -3.4627445933867795, "baseline_median": 0.6853133335681595, "candidate_median": 0.708812621088462, "baseline_p90": 0.704399386334412, "candidate_p90": 0.7252481512408075, "baseline_p95": 0.7096988325518285, "candidate_p95": 0.7302285891673248, "baseline_worst": 0.713983269096809, "candidate_worst": 0.732898245025426, "baseline_failures": 0, "candidate_failures": 0}
- StressColoredNoise/Dev_ModerateQ/EventAmplitudeError：{"baseline": 0.6097706865905674, "candidate": 0.6404054089295675, "allowance": 0.030488534329528375, "n": 36}
- StressColoredNoise/Holdout_HighContrastReversal/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "Holdout_HighContrastReversal", "n": 36, "baseline": 0.6301063795865991, "candidate": 0.6519502359638779, "improvement_percent": -3.466693416373623, "baseline_median": 0.6326926147007605, "candidate_median": 0.6552578936689946, "baseline_p90": 0.6441443045378995, "candidate_p90": 0.6650093280948896, "baseline_p95": 0.648656932643033, "candidate_p95": 0.669155087906362, "baseline_worst": 0.651716770337137, "candidate_worst": 0.67285977940395, "baseline_failures": 0, "candidate_failures": 0}
- StressColoredNoise/Holdout_HighContrastReversal/EventAmplitudeError：{"baseline": 0.5136002944747018, "candidate": 0.5472807434612899, "allowance": 0.02568001472373509, "n": 36}
- StressColoredNoise/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 0.613444603901771, "candidate": 0.6496620087200812, "allowance": 0.03067223019508855, "n": 36}
- StressColoredNoise/Holdout_StrongAtten/AVOInterceptError：{"baseline": 0.5895256232729431, "candidate": 0.625163682578458, "allowance": 0.02947628116364716, "n": 36}
- StressColoredNoise/Holdout_ThinInterbed/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "Holdout_ThinInterbed", "n": 36, "baseline": 0.6601867952343154, "candidate": 0.6873123928786554, "improvement_percent": -4.10877615852836, "baseline_median": 0.6621475024485545, "candidate_median": 0.6892962614732909, "baseline_p90": 0.6803130522774015, "candidate_p90": 0.7062273025725485, "baseline_p95": 0.6848622569519648, "candidate_p95": 0.70971543120562, "baseline_worst": 0.686327757285409, "candidate_worst": 0.712665305657129, "baseline_failures": 0, "candidate_failures": 0}
- StressColoredNoise/Holdout_ThinInterbed/EventAmplitudeError：{"baseline": 0.5927127155972729, "candidate": 0.6251989620771492, "allowance": 0.029635635779863646, "n": 36}
- StressColoredNoise/Holdout_ThinInterbed/AVOInterceptError：{"baseline": 0.5766255504609928, "candidate": 0.6077714775357131, "allowance": 0.028831277523049644, "n": 36}
- StressColoredNoise/Holdout_ThinInterbed/AVOGradientError：{"baseline": 0.6077315791429707, "candidate": 0.6389390882126998, "allowance": 0.03038657895714854, "n": 36}
- StressColoredNoise/SoftLowQ/EventAmplitudeError：{"baseline": 0.6203895374023074, "candidate": 0.6552039818924891, "allowance": 0.031019476870115373, "n": 36}
- StressColoredNoise/SoftLowQ/AVOInterceptError：{"baseline": 0.6035799780697939, "candidate": 0.6373076145575696, "allowance": 0.0301789989034897, "n": 36}
- StressColoredNoise/ThinHighFrequency/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "ThinHighFrequency", "n": 36, "baseline": 0.7531894214660506, "candidate": 0.7709171671676781, "improvement_percent": -2.3536902134288162, "baseline_median": 0.7533893255874344, "candidate_median": 0.771410096653449, "baseline_p90": 0.765676525426827, "candidate_p90": 0.782580778096519, "baseline_p95": 0.769286429090343, "candidate_p95": 0.7853240453567453, "baseline_worst": 0.771910430837063, "candidate_worst": 0.787578187843314, "baseline_failures": 0, "candidate_failures": 0}
- StressColoredNoise/VelocityReversal/mean_noninferiority：{"Suite": "StressColoredNoise", "Model": "VelocityReversal", "n": 36, "baseline": 0.5978563604746288, "candidate": 0.6200442217260159, "improvement_percent": -3.711236129322537, "baseline_median": 0.595330088669566, "candidate_median": 0.6181512173333714, "baseline_p90": 0.614897746363197, "candidate_p90": 0.635295735892486, "baseline_p95": 0.6158077489526275, "candidate_p95": 0.636379287076505, "baseline_worst": 0.617993662182543, "candidate_worst": 0.638280350937947, "baseline_failures": 0, "candidate_failures": 0}
- StressColoredNoise/VelocityReversal/EventAmplitudeError：{"baseline": 0.48104980243129164, "candidate": 0.5143402186794412, "allowance": 0.024052490121564584, "n": 36}
- StressColoredNoise/VelocityReversal/AVOInterceptError：{"baseline": 0.4885963120673107, "candidate": 0.5136184495554318, "allowance": 0.024429815603365535, "n": 36}
- StressNonRickerPhase/Dev_DeepContrast/mean_noninferiority：{"Suite": "StressNonRickerPhase", "Model": "Dev_DeepContrast", "n": 36, "baseline": 0.6970461785413431, "candidate": 0.7180087964214914, "improvement_percent": -3.0073499468880494, "baseline_median": 0.698521334253823, "candidate_median": 0.7200088784244665, "baseline_p90": 0.7159453808279226, "candidate_p90": 0.7349601550572324, "baseline_p95": 0.7191823403586628, "candidate_p95": 0.7382231959526888, "baseline_worst": 0.723846783050844, "candidate_worst": 0.742640879379643, "baseline_failures": 0, "candidate_failures": 0}
- StressNonRickerPhase/Dev_DeepContrast/EventAmplitudeError：{"baseline": 0.525246186941242, "candidate": 0.5695222589573243, "allowance": 0.026262309347062102, "n": 36}
- StressNonRickerPhase/Dev_DeepContrast/AVOInterceptError：{"baseline": 0.5095454932711367, "candidate": 0.552314058814737, "allowance": 0.02547727466355684, "n": 36}
- StressNonRickerPhase/Dev_DeepContrast/AVOGradientError：{"baseline": 0.5361556680477148, "candidate": 0.5788080550691752, "allowance": 0.026807783402385738, "n": 36}
- StressNonRickerPhase/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressNonRickerPhase", "Model": "Dev_ModerateQ", "n": 36, "baseline": 0.7774995119434641, "candidate": 0.7955700688807702, "improvement_percent": -2.324188846387343, "baseline_median": 0.77987519659109, "candidate_median": 0.798489944109203, "baseline_p90": 0.801075263674935, "candidate_p90": 0.81545608967874, "baseline_p95": 0.8023597529158675, "candidate_p95": 0.8171905344112163, "baseline_worst": 0.812413381728924, "candidate_worst": 0.826267324827494, "baseline_failures": 0, "candidate_failures": 0}
- StressNonRickerPhase/Dev_ModerateQ/EventAmplitudeError：{"baseline": 0.6584103483091689, "candidate": 0.6924182530000613, "allowance": 0.03292051741545845, "n": 36}
- StressNonRickerPhase/Holdout_HighContrastReversal/mean_noninferiority：{"Suite": "StressNonRickerPhase", "Model": "Holdout_HighContrastReversal", "n": 36, "baseline": 0.722117047796975, "candidate": 0.7365695126975219, "improvement_percent": -2.0014019811107175, "baseline_median": 0.720399863583691, "candidate_median": 0.7357020531631295, "baseline_p90": 0.735156585372071, "candidate_p90": 0.7484621435759085, "baseline_p95": 0.7360888644832271, "candidate_p95": 0.7488375640772338, "baseline_worst": 0.737390856051195, "candidate_worst": 0.750802967613043, "baseline_failures": 0, "candidate_failures": 0}
- StressNonRickerPhase/Holdout_HighContrastReversal/EventAmplitudeError：{"baseline": 0.5148252455753705, "candidate": 0.5581461161816753, "allowance": 0.025741262278768523, "n": 36}
- StressNonRickerPhase/Holdout_HighContrastReversal/AVOInterceptError：{"baseline": 0.5195806129130156, "candidate": 0.555772452321909, "allowance": 0.025979030645650783, "n": 36}
- StressNonRickerPhase/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 0.6195875304439878, "candidate": 0.6690867840222839, "allowance": 0.030979376522199394, "n": 36}
- StressNonRickerPhase/Holdout_StrongAtten/AVOInterceptError：{"baseline": 0.6065856867951095, "candidate": 0.6539418227951519, "allowance": 0.030329284339755475, "n": 36}
- StressNonRickerPhase/Holdout_StrongAtten/AVOGradientError：{"baseline": 0.7113499355613742, "candidate": 0.7576550778278113, "allowance": 0.03556749677806871, "n": 36}
- StressNonRickerPhase/Holdout_ThinInterbed/mean_noninferiority：{"Suite": "StressNonRickerPhase", "Model": "Holdout_ThinInterbed", "n": 36, "baseline": 0.7387352067237536, "candidate": 0.7617850145166639, "improvement_percent": -3.12017182653781, "baseline_median": 0.7416091537432945, "candidate_median": 0.7632708208498165, "baseline_p90": 0.7559164122214339, "candidate_p90": 0.7767742539984079, "baseline_p95": 0.75984352732183, "candidate_p95": 0.780371249864933, "baseline_worst": 0.765867545723415, "candidate_worst": 0.787131809845918, "baseline_failures": 0, "candidate_failures": 0}
- StressNonRickerPhase/Holdout_ThinInterbed/EventAmplitudeError：{"baseline": 0.6290841352441819, "candidate": 0.6661244251663321, "allowance": 0.0314542067622091, "n": 36}
- StressNonRickerPhase/Holdout_ThinInterbed/AVOInterceptError：{"baseline": 0.6187161954133497, "candidate": 0.6535842289239038, "allowance": 0.030935809770667483, "n": 36}
- StressNonRickerPhase/Holdout_ThinInterbed/AVOGradientError：{"baseline": 0.6436716665165896, "candidate": 0.6792877931892497, "allowance": 0.03218358332582948, "n": 36}
- StressNonRickerPhase/SoftLowQ/EventAmplitudeError：{"baseline": 0.6468803991412277, "candidate": 0.6911087572460192, "allowance": 0.03234401995706139, "n": 36}
- StressNonRickerPhase/SoftLowQ/AVOInterceptError：{"baseline": 0.6443418563125929, "candidate": 0.6854496219014323, "allowance": 0.03221709281562964, "n": 36}
- StressNonRickerPhase/SoftLowQ/AVOGradientError：{"baseline": 0.7257188351653601, "candidate": 0.76663275714588, "allowance": 0.036285941758268006, "n": 36}
- StressNonRickerPhase/ThinHighFrequency/mean_noninferiority：{"Suite": "StressNonRickerPhase", "Model": "ThinHighFrequency", "n": 36, "baseline": 0.8093039625392088, "candidate": 0.8265228445809991, "improvement_percent": -2.1276161786933256, "baseline_median": 0.8098065066620934, "candidate_median": 0.8271297351775735, "baseline_p90": 0.8225407536792115, "candidate_p90": 0.8378910236115089, "baseline_p95": 0.8257324309162533, "candidate_p95": 0.8409837405807302, "baseline_worst": 0.828298245765374, "candidate_worst": 0.843541808223879, "baseline_failures": 0, "candidate_failures": 0}
- StressNonRickerPhase/VelocityReversal/EventAmplitudeError：{"baseline": 0.48222824966455685, "candidate": 0.5253899808093847, "allowance": 0.024111412483227845, "n": 36}
- StressNonRickerPhase/VelocityReversal/AVOInterceptError：{"baseline": 0.4778531042538219, "candidate": 0.5138958431868914, "allowance": 0.023892655212691096, "n": 36}
- StressNonRickerPhase/VelocityReversal/AVOGradientError：{"baseline": 0.6813330167089052, "candidate": 0.7175647608790879, "allowance": 0.034066650835445265, "n": 36}
- StressResidualNMO/Dev_DeepContrast/mean_noninferiority：{"Suite": "StressResidualNMO", "Model": "Dev_DeepContrast", "n": 36, "baseline": 0.6279954773716842, "candidate": 0.6542336489966765, "improvement_percent": -4.178082895565671, "baseline_median": 0.629638575924385, "candidate_median": 0.655563782093272, "baseline_p90": 0.651402979533825, "candidate_p90": 0.674300788449072, "baseline_p95": 0.6534199642681975, "candidate_p95": 0.6772975128963028, "baseline_worst": 0.659985861154803, "candidate_worst": 0.683437057199226, "baseline_failures": 0, "candidate_failures": 0}
- StressResidualNMO/Dev_DeepContrast/EventAmplitudeError：{"baseline": 0.4927963138985242, "candidate": 0.5371828982160187, "allowance": 0.02463981569492621, "n": 36}
- StressResidualNMO/Dev_DeepContrast/AVOInterceptError：{"baseline": 0.4388264977867649, "candidate": 0.4877741066381382, "allowance": 0.021941324889338247, "n": 36}
- StressResidualNMO/Dev_DeepContrast/AVOGradientError：{"baseline": 0.3665095712197684, "candidate": 0.4321173257251207, "allowance": 0.01832547856098842, "n": 36}
- StressResidualNMO/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressResidualNMO", "Model": "Dev_ModerateQ", "n": 36, "baseline": 0.6902167890378823, "candidate": 0.7164504603458821, "improvement_percent": -3.8007871910168856, "baseline_median": 0.6899995079230785, "candidate_median": 0.7166348205061766, "baseline_p90": 0.718774159427331, "candidate_p90": 0.740955501658109, "baseline_p95": 0.7228488490486082, "candidate_p95": 0.744964551038701, "baseline_worst": 0.73522144306917, "candidate_worst": 0.755813505828756, "baseline_failures": 0, "candidate_failures": 0}
- StressResidualNMO/Dev_ModerateQ/EventAmplitudeError：{"baseline": 0.5981187136427658, "candidate": 0.6362898568818303, "allowance": 0.029905935682138292, "n": 36}
- StressResidualNMO/Dev_ModerateQ/AVOInterceptError：{"baseline": 0.5745308065764096, "candidate": 0.6122839911736041, "allowance": 0.028726540328820483, "n": 36}
- StressResidualNMO/Dev_ModerateQ/AVOGradientError：{"baseline": 0.5760832767855549, "candidate": 0.6230240839122586, "allowance": 0.02880416383927775, "n": 36}
- StressResidualNMO/Holdout_DeepLowQ/AVOGradientError：{"baseline": 0.4315197777518419, "candidate": 0.45924276389203433, "allowance": 0.021575988887592098, "n": 36}
- StressResidualNMO/Holdout_HighContrastReversal/EventAmplitudeError：{"baseline": 0.5484467584635176, "candidate": 0.5835124035579572, "allowance": 0.027422337923175884, "n": 36}
- StressResidualNMO/Holdout_HighContrastReversal/AVOInterceptError：{"baseline": 0.4695075026061066, "candidate": 0.5097616805305341, "allowance": 0.02347537513030533, "n": 36}
- StressResidualNMO/Holdout_HighContrastReversal/AVOGradientError：{"baseline": 0.25217335049838163, "candidate": 0.3472116394489483, "allowance": 0.012608667524919083, "n": 36}
- StressResidualNMO/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 0.5673041493950052, "candidate": 0.6226178423708926, "allowance": 0.02836520746975026, "n": 36}
- StressResidualNMO/Holdout_StrongAtten/AVOInterceptError：{"baseline": 0.5448607146785555, "candidate": 0.5987836952855207, "allowance": 0.027243035733927775, "n": 36}
- StressResidualNMO/Holdout_StrongAtten/AVOGradientError：{"baseline": 0.6623827722504806, "candidate": 0.7149906823181058, "allowance": 0.03311913861252403, "n": 36}
- StressResidualNMO/Holdout_ThinInterbed/mean_noninferiority：{"Suite": "StressResidualNMO", "Model": "Holdout_ThinInterbed", "n": 36, "baseline": 0.656751460763643, "candidate": 0.688127540610477, "improvement_percent": -4.777466320417645, "baseline_median": 0.660611402299767, "candidate_median": 0.6898281499058195, "baseline_p90": 0.6806447347156039, "candidate_p90": 0.7086627982855189, "baseline_p95": 0.6844281730243426, "candidate_p95": 0.7139114044361939, "baseline_worst": 0.693024738977975, "candidate_worst": 0.720419802278013, "baseline_failures": 0, "candidate_failures": 0}
- StressResidualNMO/Holdout_ThinInterbed/EventAmplitudeError：{"baseline": 0.577344425631144, "candidate": 0.6180193261338067, "allowance": 0.0288672212815572, "n": 36}
- StressResidualNMO/Holdout_ThinInterbed/AVOInterceptError：{"baseline": 0.5573648263471421, "candidate": 0.5968654385599649, "allowance": 0.027868241317357108, "n": 36}
- StressResidualNMO/Holdout_ThinInterbed/AVOGradientError：{"baseline": 0.5921429113066162, "candidate": 0.6316533947816214, "allowance": 0.02960714556533081, "n": 36}
- StressResidualNMO/SoftLowQ/EventAmplitudeError：{"baseline": 0.5892954419508643, "candidate": 0.6391678223342939, "allowance": 0.029464772097543213, "n": 36}
- StressResidualNMO/SoftLowQ/AVOInterceptError：{"baseline": 0.5689443056222102, "candidate": 0.6177068765305002, "allowance": 0.02844721528111051, "n": 36}
- StressResidualNMO/SoftLowQ/AVOGradientError：{"baseline": 0.6661307727476444, "candidate": 0.7144109205653986, "allowance": 0.03330653863738222, "n": 36}
- StressResidualNMO/ThinHighFrequency/mean_noninferiority：{"Suite": "StressResidualNMO", "Model": "ThinHighFrequency", "n": 36, "baseline": 0.7449199471141004, "candidate": 0.7677296347873995, "improvement_percent": -3.062032069575582, "baseline_median": 0.7469922199889616, "candidate_median": 0.7686235642197685, "baseline_p90": 0.7625871700552075, "candidate_p90": 0.7833112850401245, "baseline_p95": 0.7659352471121698, "candidate_p95": 0.7854551814180477, "baseline_worst": 0.770704190839552, "candidate_worst": 0.789392207075031, "baseline_failures": 0, "candidate_failures": 0}
- StressResidualNMO/VelocityReversal/EventAmplitudeError：{"baseline": 0.5137322394106348, "candidate": 0.5486841386164026, "allowance": 0.025686611970531743, "n": 36}
- StressResidualNMO/VelocityReversal/AVOInterceptError：{"baseline": 0.42701385923088014, "candidate": 0.4672770030689753, "allowance": 0.02135069296154401, "n": 36}
- StressResidualNMO/VelocityReversal/AVOGradientError：{"baseline": 0.21752711947141665, "candidate": 0.30433950921676933, "allowance": 0.010876355973570832, "n": 36}
- StressVelocityError/Holdout_ThinInterbed/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Holdout_ThinInterbed", "n": 36, "baseline": 0.7109143913922318, "candidate": 0.7354667002933676, "improvement_percent": -3.4536238397218755, "baseline_median": 0.7103034223968014, "candidate_median": 0.7346131343780604, "baseline_p90": 0.721602435128547, "candidate_p90": 0.744875226517702, "baseline_p95": 0.7277316933351152, "candidate_p95": 0.7508604902994288, "baseline_worst": 0.730802703706062, "candidate_worst": 0.753767313449789, "baseline_failures": 0, "candidate_failures": 0}
- main_equal_model_improvement：{"baseline": 0.6565302490183482, "candidate": 0.6830441250575232, "relative_improvement": -0.040384850627703583}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
