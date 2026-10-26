"""Build a GitHub-rendered gallery that also works in a phone browser."""
from pathlib import Path
import json

root = Path(__file__).resolve().parent
figures = json.loads((root / 'FIGURE_INDEX.json').read_text(encoding='utf-8'))
overview = next(x for x in figures if x['name'] == '00_research_models_overview')
ordered = [overview] + [x for x in figures if x is not overview]
lines = [
    '# 地震衰减补偿研究 · MATLAB 模型图册', '',
    '手机直接在本页向下查看图片，点图片可放大。无需下载 MATLAB，也无需打开 HTML 文件。', '',
    '**26 张图，覆盖 29 个主要模型配置。** 展示速度、Q、反射系数及配套地震道集。', '',
    '[模型来源与复现说明](README_ZH.md) · [MATLAB 源码](plot_research_models.m)', '',
]
for i, fig in enumerate(ordered, 1):
    stem = 'figures/' + fig['name']
    assert (root / (stem + '.png')).is_file()
    lines += [f'## {i}. {fig["family"]} · {fig["name"]}', '', fig['caption'], '',
              f'![{fig["family"]}：{fig["caption"]}]({stem}.png)', '',
              f'[查看原图]({stem}.png) · [PDF]({stem}.pdf) · [MATLAB FIG]({stem}.fig)', '']
(root / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
print(f'Built mobile GitHub gallery: {len(ordered)} figures')
