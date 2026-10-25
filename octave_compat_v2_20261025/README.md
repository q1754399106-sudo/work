# Octave 兼容层 v2：MATLAB 的 rand 与 randn 共用同一随机流

在 `octave_compat_20261025` 的基础上，加入 MATLAB 的 `rand`。`rand` 与 `randn` 共用 `rng(seed,'twister')` 的同一条 mt19937ar 流，因此可以运行在循环中交替调用 `rand`/`randn` 的模型生成器（例如 `make_dev.m`）。

- `rng.m`：通过 `mlwords.py` 取出原始 32 位字，在 Octave 中缓存。
- `rand.m`：每个数用两个字（genrand_res53），与 MATLAB 相同。
- `randn.m`：256 段 Ziggurat，算法见 `mlrandn.py` 文件头，在 Octave 中逐个实现。
- `sum.m`、`string.m`、`mlrandn.py`：与 v1 相同。

检验：
- `randn` 与 v1 逐位相同（F1 噪声种子，44,649 个数）；
- `rng(0)` 时 `rand` 给出 0.8147、0.9058、0.1270，与 MATLAB 相同；
- 用 `make_dev.m` 的生成器在 Octave 中重新生成 G1–G8，与存档的 MATLAB 模型相比：细层数相同，上尺度 Q 相差 ≤1.1e-13，无噪声道集相对差 ≤2.0e-15（`qscan_abstain_20261025/analysis/generator_check.json`）。

v1 保持不变，因为它的文件哈希已封存在 `b0_world_truth_20261024` 中。
