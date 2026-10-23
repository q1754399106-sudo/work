# 实际 MATLAB 验证，2026-10-17

已连接用户启动的 ResearchMATLAB，MATLAB R2024b Update 6，PID 32404。
run_rebuild、validate_extra、run_validation_suite、run_angle_ablation 均实际执行成功。
完整记录：../outputs/attenuation_matlab_review_20261017/matlab_execution.log。
完整中文理论与实验报告：../outputs/attenuation_matlab_review_20261017/REPORT_ZH.md。

- 未启动新 MATLAB、未清空工作区、未关闭用户会话。函数自动恢复 RNG 状态；临时路径已恢复。
- config_model.m 与运行前 SHA256 一致；历史工程/结果完整保存在 ../work/backups/matlab_review_20261017。
- 新主实验：864 组观测 × 5 种方法 = 4320 行结果，另有 864 行射线近似消融。
- 所有网格和种子保存于 protocol.json；不使用干净信号选择正则化参数。
- 分步误差预算最大相对残差约 1.35e-14；精确白化等价性误差约 1.81e-13。
- 同一 MATLAB AVO 真值及含噪观测的 NumPy 独立核验通过，最大相对差约 1.45e-13。
- 历史 verification/numpy_checks.json 是另一套固定振幅核验，保留原样，不能混同。

原始主程序没有发生运行错误。加入可选输出目录后首次调用仍命中 MATLAB 旧函数签名缓存；
rehash 后解决，失败记录保留于 first_suite_attempt.log。最终完整运行 stderr 为空。

结果并非分步方法全面优于一步：正确 Q 的非零偏移下 Direct 更好；零偏移固定分步有较低误差；
Q 偏差下部分条件原协方差分步更好。一致算子+精确白化与同正则化的一步反演代数等价，
不作为创新宣称。实测资料、真实 NMO、精确 Zoeppritz 及未知噪声估计仍未验证。

修改：原 run_rebuild 仅增加可选输出目录与 RNG 恢复，validate_extra 增加可选输出目录；
其余新增为独立审计、对照、核验和报告程序，原配置保持。
