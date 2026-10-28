# 实验 10：用频谱衰减估计 Q 系数（替代 QScan）

- `PLAN_ZH.md`：预登记方案；`PRE_RUN_HASHES.txt`：运行前的代码、方案和开发结果哈希。
- `scripts/sq.py`：SR 估计（相邻道互谱、多窗口拟合 t*/c、自助＋留一窗口区间）。
- `scripts/gen.py`：按实验 8 的做法生成半合成数据，可选用真实井段反射系数；`scripts/wells.py`：测井转时间域反射系数与颜色统计。
- `scripts/snapshot_inputs.py`：从登记提交恢复输入到被忽略的 `cache/inputs/`。
- `scripts/dev.py`、`scripts/dev_seeds.py`：登记前的开发（已见数据）；`scripts/smoke.py`：冒烟测试（不看分数）。
- `scripts/run.py synth K N | real`：检验数据与实测；`scripts/analyze.py`：E1–E4 判定。
