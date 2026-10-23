# 四方法同条件对比

先读[结论](CONCLUSIONS_ZH.md)，完整结果在[双语报告](REPORT_ZH_EN.md)。

本次直接重分析同一批已经实际运行的MATLAB结果，没有新增MATLAB运行，也没有重新调参。原始756行中选取四组共504行；另252行用于同结构先验的一步/分步对照。原工程未修改。

- `inputs/`：只读使用的源结果与有关MATLAB代码快照；`SOURCE_MANIFEST.json`记录来源和SHA256。
- `data/`：全部逐观测结果、六种子均值、逐情景汇总、所有六种方法配对差、参数触边和波形重算记录。
- `figures/`：9组中英图，每组PNG和PDF，均有图例。
- `ANALYSIS_CHECKS.json`：数值及配对核验。
- `VISUAL_QA_ZH.md`：人工图件检查及统计解释检查。
- `PLAN_ZH.md`：本次分析约定；不是原始实验预注册。

可直接从本包的源快照复核，不需要MATLAB，也不读取真实Q选参：

```text
python compare_fourway.py --source inputs --output 新的输出目录
```

依赖：Python 3.12、numpy 2.3.5、pandas 3.0.1、matplotlib 3.11.2、Pillow；中文字体使用Windows微软雅黑。输出目录不能已存在。本包重算的是分析、绘图及已存波形误差，不是从头重新求解MATLAB反演。

常规分步特指无多道结构、噪声自适应第二步且第一步lambda_x固定的版本，不是最初两步均固定参数的历史程序。Frozen和Target08均取分步路线。四方法共用估计有效Q；没有为一步法提供真实Q。
