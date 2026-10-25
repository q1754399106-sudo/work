# 封存后的改动记录

预测在 `predictions/SEALED_SHA256.json` 记录的时刻封存。之后只改动评分侧脚本，预测、清单、输入与预测侧代码未变（analyze.py 逐一核对哈希）。

1. `scripts/score_diagnostics.m`：每个场景的目标只重算一次（缓存），不再逐道集重算；数值不变。
2. 新增 `scripts/posthoc_b0.m`（事后分析，E1 失败之后编写）：核对自助噪声模型与真实噪声，并用真实噪声水平重算开发 Q0 道集的 B0。它不改变封存预测，也不参与预先判定。
3. `scripts/run_stage.py`：允许的阶段名加入 `posthoc_b0`。
