# Q 模型不确定时的稳健 Q 补偿偏移（UQ 逆时偏移）

**状态：阶段 1 进行中（开发）。** 方案见 [PLAN_STAGE1_ZH.md](PLAN_STAGE1_ZH.md)，评分前提交。

## 任务

把叠后已确认的 UQ-W 期望算子（`qunknown_comp_20261109`，确认 K1–K6 通过）推广到叠前两程（逆时）偏移。整体分五个阶段：

1. 逆时偏移中的 UQ 式增益，两条候选路线择一；
2. 扩展开发集：至少 4 个新模型、4 种 Q 误差形式、速度误差、相干噪声；
3. 冻结方案并确认；
4. 实测二维测线（只描述）；
5. 整理并正式查新。

判定规则 M1–M4 是草案，阈值在阶段 2 结束后定稿。基线的参数只用可行做法确定，不看参照成像。

## 内容

| 文件 | 作用 |
|---|---|
| `scripts/rtm_setup.m` | 网格、观测系统、频率、UQ 常数 |
| `scripts/make_model.m` | 开发模型 P1（倾斜层、细层、低 Q 透镜）和均匀模型 T0 |
| `scripts/fd2d_visco.m`、`gsls_design.m` | 时间域二维黏声有限差分，从 `fd2d_20261025` 原样复制 |
| `scripts/fd_shots.m` | 逐炮正演（黏声与弹性），可续跑 |
| `scripts/helm_op.m` | 频率域四阶 Helmholtz 算子，吸收带用拉伸坐标 |
| `scripts/eikonal_fs.m`、`tstar_fields.m`、`two_way.m` | 快速扫描法程函方程；按 Fermat 原理用走时对慢度扰动的导数求 t*；两程量用镜像检波点 |
| `scripts/uq_tab.m`、`gain_s.m` | UQ 期望（对数正态 Gauss–Hermite），稳定项 s |
| `scripts/prep_data.m` | 去直达波，加噪声，转到频率域 |
| `scripts/routec.m` | 路线 C：叠前逐道 UQ-W |
| `scripts/migrate.m` | 路线 A：频率域逆时偏移，成像条件中逐点施加期望算子，成像域奇偶炮信噪比 |
| `scripts/score_stage1.m` | 阶段 1 指标 |
| `scripts/test_engine.m` | 引擎检验 E1–E4 |
| `scripts/run_stage1.m` | 阶段 1 整条链 |
| `scripts/bench_*.m` | 测速 |
| `cache/`（不上传） | 炮记录、程函方程场、数据集、成像 |

## 运行

需要 MATLAB R2024b，只开一个进程。

```
matlab -batch "cd scripts; test_engine"
matlab -batch "cd scripts; run_stage1('P1')"
matlab -batch "cd scripts; score_stage1('P1')"
```
