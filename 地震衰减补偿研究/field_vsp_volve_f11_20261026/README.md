# 实验 6：Volve 15/9-F-11 T2 实测 VSP 估 Q 与 SC、QScan＋SC 检验

- `PLAN_ZH.md`：预登记方案，**草案，未冻结**。冻结前没有运行任何 Q 估计、频谱、补偿或井震标定。
- `DATA_LICENCE_ZH.md`：数据来源、许可与署名；`DOWNLOAD_MANIFEST_SHA256.txt`：原始文件校验值。
- `raw/`：Volve 原始文件（Git LFS），未作修改。`fetch_volve.sh` 按 `fetch_list.txt` 从 Databricks 卷下载。
- `probe_vsp.py` → `probe/vsp_geometry.json`：NIVSP、RSVSP 每一级的 MD、TVDSS、坐标、入射角、承包商初至，以及对应的原始炮和水听器炮（只读道头）。
- `probe_cube.py` → `probe/cube_coverage.json`：检波点到 ST0202 MIG_RAW 最近面元的距离（只读道头）。MIG_RAW 从 `field_vsp_volve_20261025/raw/seismic/` 的 LFS 对象恢复到被忽略的 `cache/seismic/`。

为什么选这口井：同一工区里，F-11 T2 的垂直入射 VSP（NIVSP，152 级）加上近垂直的平台震源 VSP 上段（RSVSP，98 级，入射角 ≤ 12.5°），初至单程跨度约 1.05 s，比实验 5 的 F-15A（0.79 s）长，而且可以分两段估 Q；井位与 F-15A 不重叠。承包商报告里没有 Q 分析。
