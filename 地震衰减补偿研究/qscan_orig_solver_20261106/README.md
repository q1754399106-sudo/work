# 在已发表的原始求解器上复现扫描 Q 的压缩

投稿缺口第 1 条。方案 `PLAN_ZH.md`（运行前提交）。

- `scripts/orig_wang2021/`：Wang et al. 2021 公开代码原样副本，来源见 `SOURCE.txt`。
- `scripts/qcdn_core.m`：`Qcdn.m` 的算子和求解部分，去掉读写。
- `scripts/scan_q.m`、`peak_freq.m`、`fw_anchor.m`、`gen_data.m`、`orig_lib.m`：外层扫描、取子波、数据。
- `scripts/run_orig.m`：工作进程，结果写入 `results/`；`scripts/analyze_orig.py`：评分，写入 `analysis/orig.json`。
- `scripts/dev_*.m`：开发用（种子 901–903），不进入判定。

运行：在 `scripts/` 下启动 4 个 `matlab -singleCompThread -batch "run_orig(i)"`，再运行 `python analyze_orig.py`。
