# 断层＋复杂 Q 半合成数据第二批

说明见 [DATA_ZH.md](DATA_ZH.md)。第一批生成器（`semisynth_fault_q_20261030/scripts/fq.py`）原样调用，只换三个新构造 D、E、F（构造线 10253、10305、10351）和种子，只保留 Z 组，c₀ 取 0.75、1.05，共 12 个数据集。为实验 16（`mu_confirm_20261031`）而建，冻结方案前全部封存。

运行：`scripts/build2.py` → `scripts/summarize.py > analysis/summary_table.md` → `scripts/figures2.py`。数据在 `cache/`，不上传。
