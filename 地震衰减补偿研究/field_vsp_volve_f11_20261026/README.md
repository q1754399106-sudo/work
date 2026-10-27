# 实验 6：Volve 15/9-F-11 T2 实测 VSP 估 Q 与 SC、QScan＋SC 检验

- `PLAN_ZH.md`：预登记方案（2026-10-27 冻结）；`PRE_RUN_HASHES.txt`：运行前的代码与方案哈希。
- `scripts/`：冻结代码，运行顺序见方案末尾。`snapshot_inputs.py` 把引用的冻结代码和实验 5 归档文件复制到被忽略的 `cache/inputs/` 并逐字节核对；`extract_logs.py` 读测井；`vsp_q.m`、`prep_line.m`、`run_field.m`、`score_field.m`、`score_f15a.m` 由 `run_stage.py` 调用 MATLAB 运行；`run_sc.py` 跑 SC 与 QScan＋SC。
- `DATA_LICENCE_ZH.md`：数据来源、许可与署名；`DOWNLOAD_MANIFEST_SHA256.txt`：原始文件校验值。
- `raw/`：Volve 原始文件（Git LFS），未作修改。`fetch_volve.sh` 按 `fetch_list.txt` 从 Databricks 卷下载。
- `probe_vsp.py` → `probe/vsp_geometry.json`：NIVSP、RSVSP 每一级的 MD、TVDSS、坐标、入射角、承包商初至，以及对应的原始炮和水听器炮（只读道头）。
- `probe_cube.py` → `probe/cube_coverage.json`：检波点到 ST0202 MIG_RAW 最近面元的距离（只读道头）。MIG_RAW 从 `field_vsp_volve_20261025/raw/seismic/` 的 LFS 对象恢复到被忽略的 `cache/seismic/`。

为什么选这口井：同一工区里，F-11 T2 的垂直入射 VSP（NIVSP，152 级）加上近垂直的平台震源 VSP 上段（RSVSP，98 级，入射角 ≤ 12.5°），初至单程跨度约 1.05 s，比实验 5 的 F-15A（0.79 s）长，而且可以分两段估 Q；井位与 F-15A 不重叠。承包商报告里没有 Q 分析。
