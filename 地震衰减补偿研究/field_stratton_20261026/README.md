# Stratton 3D 适用性核查（2026-10-26）

目的：判断 BEG 公开的 Stratton 3D（南德州陆上）能否作为实测资料，检验 QScan、SC 等补偿方法。只读了文件头、文本表和井曲线名，没有做任何 Q 估计或补偿。

## 结论：不适合作为方法判别资料，不再下载炮集

1. **没有独立 Q。** VSP 只有 Wardner GP 202 一口井。`NVSP1`/`WVSP1` 是两条非零偏移 VSP（源距 2113 ft、3660 ft），各 72 级，内容是“归一化、子波整形反褶积后的上行波”；`NVSP2`/`WVSP2` 是上行波的偏移结果。没有下行波，所以不能用谱比法或质心频移法估 Q。
2. **做不了多井合成记录。** `Stratton_wells.zip` 里 21 口井的曲线只有 GR、SP、电阻率、中子、密度（`RHOB`），**没有任何一口井有声波**，无法计算反射系数。`TABLE2.TXT` 中各井的层位双程时是地震解释结果，不是独立参照。
3. **井震参照只有一口井。** `VSP_9.TXT` 是 9 号井的单道走廊叠加（1.056–3.0 s，2 ms 采样），`TABLE3.TXT` 是同井的时深表。这与实验 4、5 的“单井、分辨力不足”属于同一问题。

炮集和道集仍在线（见下表），以后如果要做“无独立 Q 时的定性对比”可以再下载。

| 文件 | 大小 | 本地 |
|---|---:|---|
| `segy/navmerged/swath_1–4_geometry.sgy`（未处理炮集，道头含几何） | 共约 6.6 GB | 未下载 |
| `segy/cdp_gathers/*.sgy`（预处理 CMP 道集） | 约 0.9 GB/个 | 未下载 |
| `segy/processed/Stratton3D_32bit.sgy`（偏移成像） | 444 MB | 未下载 |
| `Stratton_wells.zip`、`Stratton_Achieve_Log.pdf`、`stratton_acqui_proc.pdf`（扫描件，无文字层）、`segy/VSP/*` | 共约 8.7 MB | `raw/`，校验值见 `DOWNLOAD_MANIFEST_SHA256.txt` |

## 来源与许可

- 地址：`https://s3.amazonaws.com/open.source.geoscience/open_data/stratton/`，说明页 https://wiki.seg.org/wiki/Stratton_3D_survey ，2026-10-26 下载。
- 条件（2014 年 Bob Hardage 与 Scott Tinker 授权）：使用者须声明数据由 Bureau of Economic Geology, The University of Texas at Austin 采集并提供，用于全球教育和培训。
