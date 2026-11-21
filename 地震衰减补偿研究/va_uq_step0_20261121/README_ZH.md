# 失败原因分析与 VA-UQ 第 0 步计算（2026-11-21）

正文在文档《Q 补偿：最近失败的原因与新方法（2026-11-21）》。正文副本见 REPORT_ZH.md。本目录只放脚本和输出，没有跑任何偏移。

| 脚本 | 内容 | 输出 |
|---|---|---|
| `scripts/t_r0band.py` | R0 自身 −10 dB 箱形带参照下，任何补偿的井旁相关上限 | `analysis/t_r0band.log`：ρ₀ 0.950–0.981，上限 +0.012~+0.033 |
| `scripts/t_hfocus.py` | 单频：地下偏移距聚焦估噪 H 与 B1、B300 对比（独立、横向相关、线性相干事件） | `analysis/t_hfocus.json/.log` |
| `scripts/t_hfocus2.py` | 单频：多次波型事件、缺口画弧 | `analysis/t_hfocus2.json/.log`（单频结论有误导，以 t_hfocus3 为准） |
| `scripts/t_hfocus3.py` | ±5 Hz 频带求和（与 est_noise2.m 同形）；泄漏、平稳性、多次波型偏差 | `analysis/t_hfocus3.json/.log`：H 对多次波型低估比 B300 更多 → 不实现 H |
| `scripts/t_convention.py` | 一维：合成约定（100 Hz 到时）与实测约定（25 Hz 到时）；R0、PH、UQ(100)、UQ(25)、锚定 UQ；E1、平移相关、箱形相关 | `analysis/t_convention.json/.log` |
| `scripts/t_convention_noise.py` | 实测约定下算子把噪声低估 1/10/100 倍时的箱形相关 | `analysis/t_convention_noise.json/.log` |

几何模型与 `direction_theory_20261121/scripts/t_blockflip.py` 相同（常速、二维远场格林函数、300 道 × 8 m、51 炮 × 40 m、W_deep 点），独立噪声下 B300 的数字与它逐位一致。

运行（Python 3.11 + NumPy + SciPy，每个在 2 核上不到 10 分钟）：

```
python3 -I scripts/t_r0band.py
python3 -I scripts/t_hfocus3.py analysis/t_hfocus3.json
python3 -I scripts/t_convention.py analysis/t_convention.json
python3 -I scripts/t_convention_noise.py scripts/t_convention.py analysis/t_convention_noise.json
```

相位占比（文档根因 2）由 `uq_rtm_20261117/analysis/stage2_C1..C4, D1..D4.json` 的 E1（m = 1）直接算出：(PH − R0)/(B⁺ − R0)，16 个单元中位数 0.963，范围 0.829–1.000；B⁺ − PH 最大 0.055。
