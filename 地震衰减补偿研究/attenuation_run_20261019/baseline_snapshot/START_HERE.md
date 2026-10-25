# 稳定吸收补偿：精简主版与下一阶段任务包

日期：2026-10-19。用途：交给 AI 读取、执行、比较、报告。

**当前默认：Frozen 分步；必须同时保留相同先验的一步法。** 保留的是有条件收益和明确理论作用的模块，不宣称每项在所有模型都有效。

## 用户先看什么

1. `TASKBOOK_ZH.md`：做什么、为何做、怎样判断有效。
2. `reports/current_effects.html`：已有结果和局部代价，可浏览器打开。
3. `reports/target08_vs_frozen/EFFECTS_ZH.md`：自动评价示例，演示为何不能只看部分模型改善。
4. `STATUS.json`：本次实际完成和未完成的验证范围。

## AI 执行顺序

阅读 `AI_EXECUTION.md`，先运行 `python scripts/verify_package.py` 和 `python scripts/summarize_existing.py`。
然后使用 `scripts/run_matlab.py` 检查共享会话，再验证、回归与复跑；无共享会话时不能启动新实例。
下一阶段候选算法尚未实现，按任务书 P1 开发并按独立测试协议验收。

## 内容

- `core/`：原有已验证核心依赖的精简快照；只有 Baseline/Frozen 默认入口。
- `scripts/`：共享会话连接、已有结果汇总、候选配对评价、包完整性检查。
- `evidence/`：历史 MATLAB 结果、四方法报告/图、失败分支记录和文献定位；只读参考。
- `reports/`：本包生成的结果与验证记录；不冒充新增反演。
- `METHOD_POLICY.json` / `NEXT_EXPERIMENT.json`：机器可读的方法保留策略与实验约定。

基础工具脚本只依赖 Python 标准库。实际求解还需用户已安装的 MATLAB R2024b 与兼容官方 Engine。
高级四方法重绘脚本在 evidence/fourway/，依赖 NumPy、pandas、Matplotlib、Pillow；已附历史图，不必重绘才能看结果。
MATLAB/Engine 安装本体不打包。所有运行输出应写入新的目录。

`core/inputs/`保留原始配置快照，可能包含旧候选名称；新公开入口只允许Baseline/Frozen，且未打包调参/Combined等旧实验执行器。不要从配置中看到候选就自动启用。
`scripts/validate_offline.py --output 新检查目录`可复核24组保存波形及验收器。其人工构造的正/负测试仅测试报告逻辑，不能当作科学实验。

历史报告中的外部相对链接可能指向旧工程；本包导航以本文件、任务书和 SOURCE_MANIFEST.json 为准。历史文件原样保留以便校验。
