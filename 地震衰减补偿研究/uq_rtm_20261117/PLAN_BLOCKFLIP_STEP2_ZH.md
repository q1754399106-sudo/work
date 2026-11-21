# 补充方案：分块估噪第二步（成像）的实现与判据写法（2026-11-21，第二步运行前提交）

第一步已通过，B\* = 300（`RESULTS_BLOCKFLIP_ZH.md`）。本文件只把 `PLAN_BLOCKFLIP_ZH.md` 第 4 节落到代码上，不改那里的数据、方法和门槛。

## 1. 代码

- `scripts/run_block2.m <name> <B>`：B = 1 和 B = 300 各跑一次。1/snr 取自 `cache/invblk_<name>.mat` 的 `<set>_B<B>`，存为 `cache/inv_<name>_bf<B>.mat`；B = 1 时打印它与阶段 2 `inv_<name>.mat` 的最大相对差（复现检查，应为 0）。
- 数据 m30、m40（`data_<name>.mat`）、l30（第一步生成的 `data_<name>_l30.mat`），E1，m = 0.5、0.7、1、1.4、2，真速度。
- 方法 PH、IQ30、IQ70、LP、W0、UQW0、W1、UQW1（UQW0 只为让 `score_stage2.m` 原样可用，不进判据）。LP 截止频率、W0、W1 用同一个 B 的估计。
- `migrate2.m` 只加两个可选字段（`C.dfiles` 多个数据文件、`C.invfile` 指定 1/snr 文件），不给时行为不变；`score_stage2.m` 只加可选第三个参数 C。偏移、走时场（沿用阶段 2 的 E1 场）、评分都不变。
- 输出 `cache/images2_<name>_bf<B>.mat`（不上传）、`analysis/stage2_<name>_bf<B>.json`；判定 `scripts/judge_block2.py` → `analysis/blockflip_B.json`。

## 2. 判据写法（W_deep 主窗）

d = UQW1 − B⁺，B⁺ = max(R0, PH, IQ30, IQ70, LP, W0, W1)，两者用同一个 B 的估计。M1 = m < 1 的平均，M2 = m = 1，M3 = m > 1 的最差（`score_stage2.m` 的 `M1_u1`、`M2_u1`、`M3_u1`）。“比 B = 1 高”指 d 的指标之差 [B = 300] − [B = 1]。

- **B1**：l30 上 4 个模型 ΔM2 的平均 ≥ 0.03，且 ΔM1 的平均 ≥ 0.03。
- **B2**：m30、m40 分别算，4 个模型 ΔM1 的平均 ≥ −0.02；B = 300 下每个（模型 × m30/m40）单元 M2 ≥ −0.04、M3 ≥ −0.06。
- 两条都成立才算第二步通过。只描述：UQW1 与 B⁺ 的绝对分数、LP 截止频率、m30/m40 上 B = 1 的结果与阶段 2 是否一致。

通过也只是开发结论；确认要新建 C5 起的模型，另行冻结方案。不通过就如实报告，不换 B、不调方框。

## 3. 运行

```
matlab -batch "cd scripts; run_block2('D1',1); run_block2('D1',300); run_block2('D2',1); run_block2('D2',300)"   % 进程 1
matlab -batch "cd scripts; run_block2('D3',1); run_block2('D3',300); run_block2('D4',1); run_block2('D4',300)"   % 进程 2
python scripts/judge_block2.py
```

内存约 14 GB，与其他会话共用，开 2 个进程。在本工作树的 `scripts/` 下运行，`cache/` 是指向主目录 `uq_rtm_20261117/cache` 的目录联接（不上传，清理时只删联接本身）。
