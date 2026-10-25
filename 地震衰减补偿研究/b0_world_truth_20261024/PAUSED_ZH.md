# 暂停记录：B0M 未见模型确认在生成新模型后停止

记录：2026-10-24 20:05（本地时间）。用户决定。

## 状态

- 开发阶段完成，选定 B0M、τ = 0.17076479629367916，冻结方案见 [PLAN_ZH.md](PLAN_ZH.md)，冻结哈希见 `logs/PRE_CONFIRM_HASHES.json`。
- `make_unseen` 于 19:48–19:53 完成（exit 0），生成 U1–U4 的求解侧输入与真值。
- 19:53 另一会话按原计划启动了 6 个 `run_confirm` worker，**全部在启动阶段失败**：内存/页面文件不足（exit 1、3、0xC0000005、0xC00000FD、0xC0000417），见 `logs/run_confirm_*_6_*`。
  - `predictions/` 与 `cache/outputs/` 均为空：**没有产生任何确认预测，也没有读取或计算任何 U1–U4 的真值、期望 ΔMSE 或误差。**
- `run_truth`、`seal_confirm.py`、`analyze_confirm.py` 从未运行。

## 为什么停

确认的对象与当前要用的方法不一致：

- 诊断与真值标签都比较 FullCalRefit 的 c=1 与 c=ĉ（`scripts/diagnose2.m`、`scripts/run_truth.m`），不含相干收缩。当前推荐的方法是 QScan＋相干标准版，需要的判断是“相干标准版是否还应叠加 QScan”。相干收缩本身降低方差，FullCalRefit 层面的标签不能直接迁移。
- B0 系列的前提是所用 Q 正确；上一轮在 Q 误差 −30% 时 8/12 误报。即使 B0M 通过 C1/C2，也不能单独作为组合方法的开关。
- U1–U4 一旦用于确认就不再是未见数据。保留它们给完整组合层面的诊断更有价值。

下一步改做完整组合的偏差—方差实验（开发数据 F1/F2），见 `combo_bias_variance_20261024/`。

## U1–U4 保留为未见数据

自本记录起，U1–U4 只能用于将来一个**在使用前另行冻结**的确认方案。在此之前：

- 不对 U1–U4 运行任何求解、诊断或打分；不打开 `truth/`。
- 不重新运行 `make_unseen`（脚本会拒绝覆盖）。
- 已经公开的信息只有 `logs/make_unseen.log` 中的细层数、上尺度 Q 与无噪能量比，以及本文件的哈希。

使用前用下列 SHA-256 核对文件未改动：

```
4c4d861a1941bc43055e686d8d40741019aa69609657451d0a095210010d1fb4  inputs/U1.mat
ee07dfe364fafa2a0f6d1d86cc610750b563b8973db85282114b06fa1d56130c  inputs/U2.mat
1351ad61dd801a1683da81f44aa3292bbcbcbd482fb36bcde18872220a149fb5  inputs/U3.mat
85f813dfd9ae94228d4a185d3ba91dadf0294e7430a97b985bb643143be1d297  inputs/U4.mat
e09d2156292e69a8d81f0b898958cdb041b604f061a2eab070d508a94fe7e606  truth/U1_truth.mat
f67dd068a1c0efacaab27f45ee01990a8beb899440c260952e99a1a8538c4f7f  truth/U2_truth.mat
deb180ea8818a834ec896a5eb150c6d297df112145ff09a467338f523d4405f4  truth/U3_truth.mat
1c3becbebbc7bac22a42b525f393b65701e1ac7cc7550f8611dca7523930b9b2  truth/U4_truth.mat
```

## 恢复 B0M 确认的条件

若以后仍要按原方案确认 B0M：先把失败的 `logs/run_confirm_*` 移到 `logs/failed_launch_20261024/`（`run_stage.py` 会用同名文件覆盖），在内存充足时按 `PLAN_ZH.md` 原样运行（worker 数可以减少，不影响结果），并在结果中注明本次暂停与失败的启动。U1–U4 若已被组合诊断用掉，则 B0M 需要另生成新模型。

## 后续（2026-10-25）

用户决定恢复，确认已按原方案完成（在 Octave 上运行），见 [RESUME_ZH.md](RESUME_ZH.md) 与 [RESULTS_ZH.md](RESULTS_ZH.md)。上面“U1–U4 保留为未见数据”一节自此失效。
