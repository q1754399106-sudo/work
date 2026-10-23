# 校准参数收缩：开发诊断

完整传播算子与求解目标不变，仅在有效Q校准与分层校准之间对Q和主频取几何收缩估计。未使用评分真值选择单个观测的分支。

| 条件 | 候选 | 波形误差 | 弱事件误差 |
|---|---|---:|---:|
|Main|Shrink0.00_OriginalPrior|0.774361232400815|0.675405855202219|
|Main|Shrink0.00_LightPrior|0.782391619604876|0.636261277373668|
|LowerSNR|Shrink0.00_OriginalPrior|0.83486905520105|0.77438121997268|
|LowerSNR|Shrink0.00_LightPrior|0.846054250891867|0.772696020573906|
|Main|Shrink0.25_OriginalPrior|0.756873139374539|0.68875450705971|
|Main|Shrink0.25_LightPrior|0.766288613777906|0.648495073632778|
|LowerSNR|Shrink0.25_OriginalPrior|0.822222224495056|0.782674923048737|
|LowerSNR|Shrink0.25_LightPrior|0.83553304767598|0.784974074907264|
|Main|Shrink0.50_OriginalPrior|0.731261545988494|0.688113410255141|
|Main|Shrink0.50_LightPrior|0.741584244290677|0.63386254219885|
|LowerSNR|Shrink0.50_OriginalPrior|0.804267449096502|0.770573756932717|
|LowerSNR|Shrink0.50_LightPrior|0.819438290895482|0.769513341481575|
|Main|Shrink1.00_OriginalPrior|0.657489293471086|0.940218569324604|
|Main|Shrink1.00_LightPrior|0.67233852495385|0.826330935724644|
|LowerSNR|Shrink1.00_OriginalPrior|0.754046335422527|0.932892438039747|
|LowerSNR|Shrink1.00_LightPrior|0.77213347346048|0.908948247360999|

在 Dev_ThinPairWideContrast 的本次单校准/单观测开发案例中，50%收缩＋LightPrior 的 Main 波形误差从原 Combined 的0.77436降至0.74158，弱事件误差从0.67541降至0.63386。LowerSNR 对应波形0.81944、弱事件0.76951。此结果不能替代四模型筛选或独立确认。

0%原先验、100%原先验及100%轻先验的十二个端点/条件/路线组合已与前轮实际结果对照；具体位级CSV字段比较见 ENDPOINT_REPLAY_AUDIT.json。新的四模型筛选复用已完成校准，只重新解算受影响的Section观测，并为未变化的Offset候选保留精确别名与来源。
