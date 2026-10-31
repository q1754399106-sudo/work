# 断层＋复杂 Q 半合成数据第三批

说明见 [DATA_ZH.md](DATA_ZH.md)。第一批生成器原样调用，新构造 G、H、I（构造线 10271、10325、10363），只保留 Z 组，c₀ 取 0.75、1.05，共 12 个。为实验 17（`mu_tune_20261031`）的确认阶段而建，封存。

运行：`scripts/build3.py` → `scripts/summarize.py > analysis/summary_table.md` → `scripts/figures3.py`。数据在 `cache/`，不上传。
