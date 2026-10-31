# 实验 17 方案：在开发数据上选 SC-a 的相干项权重 mu_relative，再在封存数据上确认

2026-10-31 写成，选择规则和检验方案在运行前一次写定：
- 开发阶段的选择规则和第 5 节的全部判定，在看到任何开发结果之前提交；
- 选出的值写入 `CHOICE.json`，与哈希一起在检验运行前另行提交。

## 1. 背景

实验 16 在 12 个封存数据集上否定了"mu_relative 从 10 改成 0"：两者无法区分。事后看，中间值 1 在 SR 路线上 12/12 好于 0 和 10，平均 S 为 0.900，高于 70 dB 反 Q 的 0.873。

但 mu 1 是在检验数据上看到的，而且只试了一个中间值。正确的做法是：在开发数据上按网格选，冻结后在新数据上确认。

## 2. 数据

| 阶段 | 数据 | 个数 |
|---|---|---|
| 开发 | 第一批（`semisynth_fault_q_20261030`）Z 组 A、B、C，c₀ = 0.75、1.05；第二批（`semisynth_fault_q2_20261031`）D、E、F。都已用过 | 24 |
| 检验 | 第三批（`semisynth_fault_q3_20261031`）G、H、I，c₀ = 0.75、1.05，新建并封存 | 12 |

- 第一批不取 c₀ = 1.4：那一档噪声被放大 2.3–3.3 倍，会造成实测上不存在的带宽塌缩。
- 第三批构建时没有运行任何补偿或选 Q 方法。

## 3. 开发阶段（`run17.py dev`、`select17.py`）

**每个开发数据集上计算**：
- 用 SR 系数 c_SR 和 QScan 系数 c_QS，分别跑 SC-a（a = 2⁻¹³）；
- mu_relative 取网格 {0, 0.25, 0.5, 1, 2, 4, 10}；
- 另算 IQ70 作参照。

**选择规则**：每条路线（SR、QT）分别取 24 个数据集平均 S 最高的网格值；并列时取较大的 mu（更接近原值）。

**另报（描述）**：按批次、构造、c₀ 分开时选出的值，以及每个 mu 胜过 mu 10 的数据集个数。

选择结果写入 `CHOICE.json`。若某路线选中 10，该路线的确认判定记为"不适用"。

## 4. 检验阶段方法（`run17.py test`，参数来自实验 12，不再调）

| 名称 | 内容 |
|---|---|
| Input | 不补偿 |
| IQ70 | 70 dB 稳定化反 Q，c_SR |
| VDB | 相位反 Q＋时变 Wiener，c_SR（描述） |
| SR_mu10 / **SR_sel** | SC-a，a = 2⁻¹³，c_SR，mu_relative = 10 / 开发选出的值 |
| **QSCa** | 用户原流程：c_QS，a = 0.125，mu 10 |
| QT_mu10 / **QT_sel** | c_QS，a = 2⁻¹³，mu_relative = 10 / 开发选出的值 |

指标 S 与实验 12、15、16 相同。

## 5. 判定（12 个封存数据集）

配对差的 90% t 区间，自由度 11。"支持"要求同时满足：区间下限 > 0；G、H、I 三个构造的平均差都 > 0。"相反"为区间上限 < 0。其余为"无法判定"。

| 判定 | 比较 | 含义 |
|---|---|---|
| **C1** | S(SR_sel) − S(SR_mu10) | SR 路线上选出的 mu 胜过原值 |
| **C2** | S(QT_sel) − S(QT_mu10) | 用户 QScan 路线（调小 a）上选出的 mu 胜过原值 |
| **C3** | S(SR_sel) − S(IQ70) | 选好 mu 的 SC-a 胜过 70 dB 反 Q |
| **C4** | S(QT_sel) − S(QSCa) | 对用户的总建议（QScan＋调小 a＋选出的 mu）胜过原流程 |

**另报（描述）**：
- QT_mu10 − QSCa，即只调小 a 的收益；
- QT_sel − IQ70；
- SR_sel − VDB；
- SR_mu10 − IQ70；
- 按 c₀ 分开的平均差；
- 带宽塌缩次数；
- QScan 选在网格边上的次数。

## 6. 解释规则

- **C1、C2 支持**：SC-a 的默认 mu_relative 改为选出的值，SR 路线和用户的 QScan 路线都适用。改动前还要在实验 13 的实测测线上做井旁检查，那是另一个实验。
- **C3 支持**：SR 路线优先推荐 SC-a（选出的 mu），不再是 70 dB 反 Q。不支持则维持原推荐。
- **C4**：对用户总建议的收益幅度。

## 7. 运行

按用户要求只开 4 个进程。

```
python scripts/run17.py dev K 4        # K = 0..3，开发阶段
python scripts/select17.py             # CHOICE.json、analysis/dev_summary.json
python scripts/hash_prerun.py          # PRE_RUN_HASHES.txt，与 CHOICE.json 一起提交
python scripts/run17.py test K 4       # 检验阶段
python scripts/analyze17.py            # analysis/summary.json
```
