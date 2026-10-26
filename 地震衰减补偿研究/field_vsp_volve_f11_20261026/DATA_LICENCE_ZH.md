# 数据来源与许可

`raw/` 中的全部文件来自 **Equinor Volve Data Village**，版权归 Equinor ASA 及前 Volve 许可方（ExxonMobil Exploration & Production Norway AS、Bayerngas Norge AS）。

- 获取途径：Databricks Marketplace，"Volve Data Village"（提供方 Equinor ASA），卷 `equinor_asa_volve_data_village.public.volve`，2026-10-26 下载（`fetch_volve.sh`，清单 `fetch_list.txt`）。
- 许可：Equinor Open Data Licence（基于 CC BY 4.0，另有两点改动：不得出售；许可覆盖数据集中的全部数据，无论其是否受版权法保护）。条款全文见 `raw/Volve_licence_terms.PDF`，以及 https://www.equinor.com/energy/volve-data-sharing 。
- 本目录中的原始文件未作修改；校验值见 `DOWNLOAD_MANIFEST_SHA256.txt`。
- 本目录中由这些数据计算得到的结果属于衍生材料（Adapted Material），按同一条款共享。它们不代表 Equinor 或前 Volve 许可方的观点，也未经其认可。

| 本地路径 | Volve 卷中的原路径 |
|---|---|
| `raw/vsp_velocity/` | `Well_logs_pr_WELL/15_9-F-11 T2/08.VSP_VELOCITY/`（除 `VSPNI_COMPUTED_2–38`、`VSPRS_COMPUTED_1–41` 两组处理成果 SEGY 外的全部文件） |
| `raw/composite/` | `Well_logs_pr_WELL/15_9-F-11 T2/04.COMPOSITE/` |
| `raw/reports/` | `Well_logs_pr_WELL/15_9-F-11 T2/14.DIV.REPORTS/` 中的 VSP 采集方案和完井图 |
| `raw/checkshots/` | `Seismic/VSP/Checkshots/checkshot_15_9_F_11T2.ASC`、`Seismic/VSP/README.txt` |
| `raw/Volve_licence_terms.PDF` | `HRS and Terms and conditions for license to data - Volve.pdf` |

地面地震沿用实验 5 已归档的 ST0202R08 MIG_RAW（`field_vsp_volve_20261025/raw/seismic/`，Git LFS），不重复下载。
