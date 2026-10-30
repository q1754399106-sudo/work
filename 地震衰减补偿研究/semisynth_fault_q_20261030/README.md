# 新半合成数据：带断层、Q 结构复杂（2026-10-30）

数据说明、构建细节、使用规则见 [DATA_ZH.md](DATA_ZH.md)。

- 取代 F3/F4、V5 以及实验 8–14 的构造线作为今后检验用的未见数据。
- 3 条从未用过的 Volve 构造线（10243、10291、10345），每条 3 条正断层（地层缺失、破碎带、生长断距）；Q 挂在断开的地质层位上：模型 Q 台阶被移动、随机分层扰动、断块 Q 乘子、横向趋势、Q = 18 的气藏透镜，所以真 t* 不是 c × 模型 t*。
- Z 组（复杂 Q）18 个，孪生 Y 组（同断层、同噪声，旧配方的横向不变 Q）18 个。
- 建议 A 构造作开发，B、C 封存到检验方案冻结之后。

| 路径 | 内容 |
|---|---|
| `scripts/envf.py` | 从实验 13 的输入快照导入实验 7–13 的程序（不改动） |
| `scripts/fq.py` | 断层几何、Q 结构、按道 t* 正演、构建一个数据集；`fq.load(编号)` 读取 |
| `scripts/test_fq.py` | 正演与断层映射自检 |
| `scripts/build.py` | 构建全部数据集 → `cache/data/`（本地，可复现）、`analysis/datasets.json`、`geology.json`、`realism.json` |
| `scripts/summarize.py` | `analysis/summary.csv`（数据层面统计） |
| `scripts/figures.py` | `figures/overview.png`、`figures/tstar.png` |
