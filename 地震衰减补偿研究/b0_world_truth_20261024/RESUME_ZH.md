# 恢复记录：B0M 未见模型确认在 Octave 上按原方案运行

登记：2026-10-25（本地时间），在对 U1–U4 运行任何诊断之前。冻结哈希见 `logs/RESUME_ADDENDUM_SHA256.json`。

## 决定

用户在 2026-10-25 的会话中要求“把暂停的继续”。按 [PAUSED_ZH.md](PAUSED_ZH.md) 末节的恢复条件，B0M 确认按 [PLAN_ZH.md](PLAN_ZH.md) 原样运行，变量、阈值 τ = 0.17076479629367916、K、噪声窗口、种子和判定 C1/C2 都不变。

暂停时的两点顾虑仍然成立，结论中会照实说明：

- 诊断与真值标签比较的是 FullCalRefit 的 c=1 与 c=ĉ，不是“相干标准版是否还应叠加 QScan”；
- B0 系列的前提是所用 Q 正确。

**代价**：本次确认之后，U1–U4 不再是未见数据。将来做组合层面的确认，需要另外生成新模型。

## 冻结文件核对

- `logs/PRE_CONFIRM_HASHES.json` 中的 36 个文件：
  - 27 个逐字节一致；
  - `analysis/dev_selection.json` 与 6 个 `analysis/dev_w*.csv` 只差行尾（Windows 工作副本为 CRLF，Git 存 LF），换成 CRLF 后哈希一致；
  - `cache/SNAPSHOT.json` 按原脚本重新生成。它引用的 103 个输入文件在冻结提交 b350377 与当前 HEAD 之间没有任何改动。
- U1–U4 的 8 个 `.mat` 与 PAUSED_ZH.md 中的 SHA-256 一致。

## 唯一的偏离：运行平台

这个云端会话里没有 MATLAB。确认改在 **GNU Octave 8.4.0（Linux，OpenBLAS，单线程，与 `bw_setup` 中 `P.threads=1` 相同）** 上运行。

- **数值代码不变**：`scripts/` 与快照中 `bae_q_20261024/core|source` 的数值代码一行未改。
- **兼容层** `octave_compat_20261025/`，只补 Octave 缺的功能：
  - `rng.m` / `randn.m`：重建 MATLAB R2024b 的 `rng(seed,'twister'); randn`（mt19937ar + 256 段 Ziggurat，见 `mlrandn.py` 文件头）。`validate_randn.py` 在仓库中 MATLAB 生成的 4 个 F1/F2 噪声场（各 44,649 个数）上检验，最大偏差 8.7e-15，属于存储精度。
  - `sum.m`：支持 `'all'` 与 `[1 2]`。
  - `string.m`：只用于 `switch lower(string(...))`。
- **驱动脚本** `scripts_octave/`：只替换输入输出部分（Octave 没有 `table`），计算调用与 `scripts/` 中对应文件逐句相同。
  - `run_confirm_oct.m` 对应 `run_confirm.m`。
  - `run_truth_oct.m` 对应 `run_truth.m`。差别只有一处：r=1–3 先重新生成观测，与 `inputs/` 中的诊断道集比较（相对偏差 ≤1e-12，不要求逐位相同），然后使用存档的道集。这样，真值所用的数据与预测时看到的数据完全相同。
- worker 数改为 4，方案允许，不影响结果。
- 原 6 个失败启动的日志已移到 `logs/failed_launch_20261024/`。

### 平台等价性检验（阶段 0，开发数据，已见）

在全部 44 个开发道集上，用 Octave 重跑 `diagnose2`，与 MATLAB 的 `analysis/dev_w*.csv` 比较，详见 `octave_check/DEV_COMPARISON.json`。检验标准在看到第一个道集（g065）后写定：

- QScan 因子逐个相同；
- B0 与 B0M 的最大绝对差 ≤1e-9；
- B0M 在 τ 下的 44 个判定全部相同；
- B0C 走相干收缩，对 BLAS 求和顺序敏感（`combo_bias_variance_20261024` 记录过 1 线程与 4 线程相差约 3e-3），只报告差值，不设门槛。

结果（冻结前填入）：**通过**。

| 项目 | Octave 与 MATLAB 的差别 |
|---|---|
| QScan 因子 | 44/44 相同（最大差 4.7e-15） |
| B0 | 最大绝对差 3.3e-14 |
| B0M | 最大绝对差 5.0e-14 |
| 噪声模型 φ̂、门限增益 G、m* 的能量比 | ≤1.9e-14 |
| B0M 在 τ 下的判定 | 44/44 相同 |
| B0 在 0 下的判定 | 44/44 相同 |
| B0C（只报告） | 最大绝对差 4.3e-3，在 0 下的判定 44/44 相同 |
| 平均耗时 | 每个道集 55 s（MATLAB 为 284 s） |

复制观测道集用重建的 `randn` 生成，与 MATLAB 一致到上述精度。这说明在本研究用到的数值路径上，重建的随机数与 MATLAB 相同。

## 封存

- 先运行未改动的 `scripts/seal_confirm.py`，再运行 `scripts_octave/seal_octave_extra.py`，补充封存 Octave 侧文件，生成 `predictions/SEALED_OCTAVE_SHA256.json`。两步都在读取任何真值之前完成。
- 评分用未改动的 `scripts/analyze_confirm.py`，通过 `scripts_octave/analyze_confirm_fonts.py` 调用，只为图中的中文补一个字体。

## 追加的描述性分析（不作门槛）

为回应暂停时的第一点顾虑，真值计算之后另跑 `scripts_octave/run_combo_truth_oct.m`：

- 在同样的 12 个实现上，计算相干标准版在 Q0 处（B）与在 QScan 因子处（D = QScan＋标准版）的误差；
- 报告冻结的 B0M 判定与组合层面期望符号的一致率。

这一项只作描述，不影响 C1/C2 的判定。
