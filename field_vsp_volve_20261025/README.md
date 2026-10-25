# 实验 5：Volve 15/9-F-15A 实测 VSP 估 Q 与补偿方法比较

- `PLAN_ZH.md`：预登记方案（提交 4285d0d）；`PLAN_ADDENDUM_ZH.md`：运行前补充；`PRE_RUN_HASHES.txt`：运行前的代码与方案哈希（提交 28cea18）。
- `RESULTS_ZH.md`：结果与判定。
- `DATA_LICENCE_ZH.md`：数据来源、许可与署名；`DOWNLOAD_MANIFEST_SHA256.txt`：原始文件校验值。
- `raw/`：Volve 原始文件（Git LFS），未作修改。
- `probe_*.py`、`probe/`：登记前的几何探查（只读道头）。
- `scripts/`：
  - `snapshot_inputs.py`：把冻结代码从登记的提交复制到被忽略的 `cache/inputs/`，并逐字节核对与 HEAD 一致；
  - `run_stage.py`：以隐藏 MATLAB 进程运行一个阶段，日志写入 `logs/`；
  - `vsp_q.m`（E5-1）→ `analysis/vsp_q.json`、`analysis/vsp_spectra.mat`；
  - `prep_line.m` → `inputs/curtain.mat`、`inputs/curtain_meta.json`；
  - `run_field.m` → `results/<run>.mat`（V0、Vm50、Vm30、Vp30、Vp50、R0）；
  - `score_field.m` → `analysis/scores.csv`、`analysis/gates.json`、`analysis/references.mat`。

运行顺序（MATLAB R2024B，4 线程）：

```
python scripts/snapshot_inputs.py
python scripts/run_stage.py vsp_q
python scripts/run_stage.py prep_line
python scripts/run_stage.py run_field V0 Vm30 Vp30
python scripts/run_stage.py run_field Vm50 Vp50 R0
python scripts/run_stage.py score_field
```
