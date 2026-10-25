# 在 GNU Octave 上运行本仓库 MATLAB 实验的兼容层

2026-10-25 建立。用途：在没有 MATLAB 的环境（例如 Linux 云端容器）中，原样运行本仓库的 MATLAB 数值代码，并且结果可以与 MATLAB R2024b 逐值对照。

## 内容

| 文件 | 作用 |
|---|---|
| `rng.m`、`randn.m` | 实现 MATLAB 的 `rng(seed,'twister'); randn(...)`，即 mt19937ar 发生器加默认的 Ziggurat 正态变换。实际计算由 `mlrandn.py` 完成 |
| `mlrandn.py` | 由仓库中 MATLAB 生成的噪声重建出的 `randn` 算法，算法说明见文件头 |
| `validate_randn.py` | 检验 `mlrandn.py`：从 `bae_q_20261024` 的 F1/F2 观测减去无噪声道集，反推出 4 组 MATLAB `randn`（各 44,649 个数），逐个比较 |
| `sum.m` | 支持 `sum(x,'all')` 与 `sum(x,[1 2])`。Octave 8 不支持前者；对后者会静默给出错误结果 |
| `string.m` | 只支持 `switch lower(string(...))` 这类用法 |

## 用法

```matlab
addpath('octave_compat_20261025','-begin');   % 必须在仓库代码之前
warning('off','all');                         % 屏蔽“遮蔽内置函数”的提示
```

需要 Python 3 与 NumPy（`randn.m` 通过 `python3 -I mlrandn.py` 生成随机数）。

## 已验证的范围

- `validate_randn.py`：4 组 MATLAB 噪声场，最大偏差 8.7e-15（存储精度）。
- `b0_world_truth_20261024/octave_check/DEV_COMPARISON.json`：44 个开发道集上重跑 `diagnose2`，并与 MATLAB 结果比较：
  - QScan 因子与 FullCalRefit 路径（B0、B0M）一致到 ≤5e-14；
  - 相干收缩路径（B0C）差到 4e-3，原因是它对 BLAS 求和顺序敏感（MATLAB 自身 1 线程与 4 线程之间也差约 3e-3）。

## 限制

- 只实现了 `rng(seed,'twister')` 之后的 `randn`。没有实现 `rand`、`randi`、`randperm`、`rng('default')` 与 `rng(state)`。如果代码用到它们，Octave 会使用自己的发生器，与 MATLAB 不同（生成新模型的 `make_*.m` 就用到 `rand`）。
- 每次调用 `randn` 都会启动一次 Python，适合整块生成噪声，不适合在循环里逐个取数。
- 涉及相干收缩的结果与 MATLAB 只在约 1e-3 的量级上一致，不是逐位一致。
