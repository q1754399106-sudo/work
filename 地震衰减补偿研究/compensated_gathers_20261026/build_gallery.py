from pathlib import Path
import json,csv

root=Path(__file__).resolve().parent
cases=json.loads((root/'CASE_INDEX.json').read_text(encoding='utf-8'))
figures=json.loads((root/'FIGURE_INDEX.json').read_text(encoding='utf-8'))
figure_names={r['stem'] for r in figures}
main=[r for r in cases if not r['case'].endswith('_Qminus30')]
extra=[r for r in cases if r['case'].endswith('_Qminus30')]
assert len(main)==29 and len(extra)==6
assert len({r['model'] for r in main})==29
assert all(r['case'] in figure_names for r in cases)
lines=['# MATLAB 补偿道集图册（手机可直接看）','','覆盖上一册的 **29 个主要模型配置**，另给出 H1–H6 的 Q 低估 30% 场景；共 35 个场景、38 张图。', '',
       '每张主图：**无衰减目标 → 含衰减观测 → 补偿后 → 误差道集**，并附两道原始幅度波形叠加。点击图片或“原图”可放大，无需 MATLAB。', '',
       '使用项目当前推荐规则：给定参考 Q 用 **SC**，四层标定 Q 不确定或 Q 低估时用 **QScan＋SC**。这些是当前推荐方案，不声称每个模型上都是所有方法中的最低误差；V3/V4/V5 另给出自研方法对照。', '',
       '[先看 F3](figures/F3.png) · [先看复杂模型 V5](figures/V5.png) · [先看二维模型 H1](figures/H1.png) · [原模型图册](../model_visualization_20261026/README.md)', '',
       '[模型、数据与复现说明](README_ZH.md) · [全部指标 CSV](metrics.csv) · [MATLAB 绘图源码](plot_gathers.m)', '',
       '全部场景采用一个固定的名义 14 dB 噪声实现。V3/V5 还含附加相干噪声、多次波等干扰，14 dB 指原始随机噪声设定，并非总干扰信噪比。四层十配置是新生成的 FaultSection 展示观测；其余沿用原观测。F1/F2 为 SC 选参数据。图中相对波形误差为 ‖补偿−目标‖F / ‖目标‖F，越低越好；复杂模型死道保留，全道指标包含死道，另在 CSV 给出仅活道指标。', '',
       '统一原始幅度色标，无自动增益（AGC）、无逐道归一化、无分位数裁幅；残差也使用同一色标。SC 参数全局固定，没有对本册新模型调参。', '',
       '## 快速索引', '', '| 模型 / 场景 | 方法 | 观测误差 | 补偿误差 | 图 |','|---|---|---:|---:|---|']
for r in main+extra:
    name=r['case'];lines.append(f'| {name} | {r["method"]} | {r["input_error"]:.3f} | {r["output_error"]:.3f} | [查看](figures/{name}.png) |')
groups=[('一、四层模型',[r for r in main if r['family']=='四层模型']),
        ('二、Kennett 细层 F1–F4 / U1–U4',[r for r in main if r['family']=='Kennett 细层']),
        ('三、V1–V5 反射与复杂地层',[r for r in main if r['model'].startswith('V')]),
        ('四、二维模型 H1–H6：参考 Q',[r for r in main if r['family']=='二维有限差分']),
        ('五、H1–H6：Q 低估 30%，使用 QScan＋SC',extra)]
for heading,rows in groups:
    lines+=['',f'## {heading}','']
    for r in rows:
        stem=r['case'];flag='；新生成展示观测' if r['family']=='四层模型' else ''
        lines += [f'### {stem} · {r["method"]}', '',
                  f'相对波形误差：**{r["input_error"]:.3f} → {r["output_error"]:.3f}**{flag}。', '',
                  f'![{stem} 目标、观测、补偿及误差道集](figures/{stem}.png)', '',
                  f'[原图](figures/{stem}.png) · [PDF](figures/{stem}.pdf) · [可编辑 FIG](figures/{stem}.fig)', '']
        comp='compare_'+stem
        if comp in figure_names:
            lines += ['**同一观测上的自研方法对照（包括表现更好的候选）**', '',
                      f'![{stem} 固定候选方法对照](figures/{comp}.png)', '',
                      f'[对照原图](figures/{comp}.png) · [PDF](figures/{comp}.pdf)', '']
(root/'README.md').write_text('\n'.join(lines).rstrip()+'\n',encoding='utf-8')
print('Mobile gallery built:',len(cases),'cases,',len(figures),'figures')
