"""Recompute the existing MATLAB metric summary; no inverse solve is run."""
from pathlib import Path
import csv, html, json, statistics

ROOT=Path(__file__).resolve().parents[1]
def main():
    with (ROOT/'evidence/historical/metrics.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    methods=[('Baseline','Direct','一步 / Direct'),('Baseline','Split','常规分步 / Split'),
             ('Frozen','Split','Frozen'),('Tuned','Split','Target08')]
    models=sorted({r['Model'] for r in rows})
    groups=[]
    for model in models:
        g={'Model':model}
        for m,route,label in methods:
            rr=[r for r in rows if r['Model']==model and r['Method']==m and r['Route']==route and r['Suite']=='Main']
            assert len(rr)==30
            g[label]=statistics.mean(float(r['RelativeError']) for r in rr)
        groups.append(g)
    aggregate={label:statistics.mean(g[label] for g in groups) for _,_,label in methods}
    gain=1-aggregate['Frozen']/aggregate['常规分步 / Split']
    assert abs(gain-.064762)<.0001
    out=ROOT/'reports'; out.mkdir(exist_ok=True)
    (out/'current_summary.json').write_text(json.dumps({'scope':'existing MATLAB results reanalysis only','rows':len(rows),
        'groups':groups,'equal_model_mean':aggregate,'frozen_vs_split_relative_improvement':gain},ensure_ascii=False,indent=2),encoding='utf-8')
    body='<h1>当前方法效果 / Current evidence</h1><p><b>推荐研究主版：Frozen；保留相同先验的一步对照。</b></p>'
    body+='<p>来源为既有 MATLAB 实验。本页只重算已存指标，没有新增反演。18 dB、三模型、五场景、六观测噪声种子；每模型只使用一份校准估计。下表相对误差越小越好，不是恢复正确率。</p>'
    body+='<table><tr><th>模型 / Model</th>'+''.join('<th>'+label+'</th>' for _,_,label in methods)+'</tr>'
    for g in groups+[{'Model':'等模型平均 / Equal mean',**aggregate}]:
        body+='<tr><td>'+html.escape(g['Model'])+'</td>'+''.join('<td>'+f'{g[label]:.6f}'+'</td>' for _,_,label in methods)+'</tr>'
    body+='</table><p>Frozen 相对常规分步平均误差降低 '+f'{100*gain:.2f}%'+'。这属于有限合成证据，不能推断普遍改善。</p>'
    body+='<p><b>保幅代价：</b>薄层弱反射误差，常规分步1.079361，Frozen 1.178977，Target08 1.355316。不能只用整体平均判断效果。给一步法相同先验后，两条路线的模型平均差异最大约0.000618。</p>'
    body+='<p><a href="target08_vs_frozen/effects.html">自动验收示例：Target08 对比 Frozen</a> · <a href="../TASKBOOK_ZH.md">下一阶段任务书</a> · <a href="../evidence/fourway/REPORT_ZH_EN.md">完整历史报告</a></p>'
    for name,caption in [('01_overall_four_methods.png','整体误差 / Overall error'),('03_amplitude_and_local_risks.png','振幅与局部代价 / Amplitude and local risks'),('09_same_prior_direct_split.png','同先验对照 / Matched-prior control')]:
        body+='<h2>'+caption+'</h2><p>以下为已归档历史图 / Archived historical figure</p><img alt="'+caption+'" src="../evidence/fourway/figures/'+name+'">'
    (out/'current_effects.html').write_text('<!doctype html><meta charset="utf-8"><title>当前补偿效果</title><style>body{max-width:1100px;margin:32px auto;padding:18px;font:16px/1.7 system-ui;color:#203040}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd;padding:10px}img{max-width:100%;height:auto}a{color:#075aa1}</style>'+body,encoding='utf-8')
    print(json.dumps({'status':'PASS','historical_rows':len(rows),'relative_improvement':gain}))

if __name__=='__main__': main()
