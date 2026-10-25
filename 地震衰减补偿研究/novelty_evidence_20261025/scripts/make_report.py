from common import *
import csv,datetime
rows=json.loads((ROOT/'analysis/test_rows.json').read_text(encoding='utf-8'))
sel=json.loads((ROOT/'analysis/selection.json').read_text(encoding='utf-8'))
checks=json.loads((ROOT/'analysis/test_checks.json').read_text(encoding='utf-8'))
theory=json.loads((ROOT/'analysis/theory_checks.json').read_text(encoding='utf-8'))
audit=json.loads((ROOT/'analysis/archive_audit.json').read_text(encoding='utf-8'))
base=f"SC_a{sel['a']:g}_m{sel['mu_relative']:g}"
ablation=f"SC_a{sel['no_structure_a']:g}_m0"
cases=list(dict.fromkeys(r['case'] for r in rows));methods=list(dict.fromkeys(r['method'] for r in rows))
values={c:{r['method']:r['error'] for r in rows if r['case']==c} for c in cases}
summary={m:dict(mean=float(np.mean([values[c][m] for c in cases])),worst=float(max(values[c][m] for c in cases)),wins_over_SC=sum(values[c][m]<values[c][base] for c in cases),
               worst_ratio_to_SC=float(max(values[c][m]/values[c][base] for c in cases))) for m in methods}
ratio=summary['Standard']['mean']/summary[base]['mean']
gates=dict(mean_at_least_5pct=ratio<=.95,at_least_6_of_8=summary['Standard']['wins_over_SC']>=6,worst_no_more_than_5pct=summary['Standard']['worst_ratio_to_SC']<=1.05)
gates['all_pass']=all(gates.values())
writejson(ROOT/'analysis/summary.json',dict(selection=sel,summary=summary,gates=gates,mean_ratio_standard_to_SC=ratio,scope='Previously seen models; descriptive, not fresh confirmation.'))
with (ROOT/'analysis/comparison.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=['case']+methods);w.writeheader()
    for c in cases:w.writerow(dict(case=c,**values[c]))
lines=['# 创新性证据：数学边界、文献基线与机制核查','',
'## Material Passport','',
'- 日期：2026-10-25；上游固定提交198f04b。',
'- Verification Status：ANALYZED。实际执行理论数值核验、算子/梯度检查、归档算术复算，以及文献目标的独立实现对照。',
'- 范围：已见模型描述性比较；没有独立研究者复现，没有实测确认，也不是全球首创证明。',
'- 方法：基线SC为马雄等2023式17结构的加权独立实现（FISTA），不是作者原版。详见METHOD_COMPARISON_ZH.md。','',
'## 本轮最重要的判断','',
('本轮描述性对照支持Standard相对所实现的SC基线具有稳定优势，达到评估前固定的三条门槛。' if gates['all_pass'] else '本轮没有同时达到预先固定的三条优势门槛，不能据此声称已经证明Standard优于最接近的空间预测约束基线。'),
f"Standard平均误差/SC平均误差={ratio:.6f}，Standard优于SC {summary['Standard']['wins_over_SC']}/8场景，最坏逐场景比值={summary['Standard']['worst_ratio_to_SC']:.6f}。",'',
'这是实际计算结果；不因结果更换主候选、挑选模型或更改门槛。','',
'## 文献基线与公平性','',
f"开发使用F1/F2的4场景，最终SC全局参数a={sel['a']}, mu_relative={sel['mu_relative']}；无结构消融a={sel['no_structure_a']}。先有开发选择文件，再评分F3/F4/V4/V5。",'',
('所选参数仍处搜索边界；本轮只比较明示开发预算内的基线，未证明其已达全部参数空间的最好性能。' if sel['a'] in [.125,1.] or sel['mu_relative'] in [.01,10.] else '所选参数位于本轮搜索网格内部。'),'',
'初始搜索最优位于边界，故仅在开发阶段扩大网格；记录保存在DEV_ADDENDUM_ZH.md与初始源码/日志中。基线始终使用相同Q、子波、几何、噪声尺度和同一弹性目标。所有SC结果均报告KKT残差，不以未收敛基线制造优势。', '',
'## 8场景结果（相对波形误差，越低越好）','',
'| 场景 | SC文献目标 | 无结构项 | FullCalRefit | Soft | v1 | Standard | QScan+Standard归档 |',
'|---|---:|---:|---:|---:|---:|---:|---:|']
for c in cases:
    x=values[c];lines.append('| '+c+' | '+' | '.join(f'{x[m]:.6f}' for m in [base,ablation,'FullCalRefit','Soft','v1','Standard','QScan_StandardArchived'])+' |')
lines+=['','| 方法 | 平均误差 | 最坏误差 | 优于SC场景 |','|---|---:|---:|---:|']
for m in methods:
    s=summary[m];lines.append(f"| {m} | {s['mean']:.6f} | {s['worst']:.6f} | {s['wins_over_SC']}/8 |")
maxplat=max(checks[c]['standard_platform_difference'] for c in cases)
allsolvers=[v for k,v in checks.items() if '/' in k]
lines+=['','## 数值检查','',
f"- Python对归档FullCalRefit与Soft的重建误差最大为{max(max(checks[c]['full_reconstruction'],checks[c]['soft_reconstruction']) for c in cases):.3g}。",
f"- Standard与归档MATLAB结果的最大相对输出差为{maxplat:.3g}；两套Standard指标同时保留，避免把平台差异藏入比较。",
f"- 评估基线求解{len(allsolvers)}次，全部收敛={all(x['converged'] for x in allsolvers)}；最大KKT/惩罚残差={max(x['kkt_ratio'] for x in allsolvers):.6g}（预定容差0.02）。",
'- 预测算子伴随、目标梯度有限差分、物理反射率坐标变换检查见analysis/operator_checks.json。','',
'## 可证明的内容','',
'当前增益g=E/(E+lambda V)逐样本精确最小化可观测二次代理目标(1-u)^2 E+u^2 lambda V；这给出明确的优化解释，但该代理目标不是无偏MSE估计，不能冒充新的Wiener原理。',
'THEORY_ZH.md给出完整条件、证明和反例。确定性结论包括每组增益在[0,1]、不翻转系数符号、固定候选支撑的零系数保持为零，以及固定E/V时相干越高增益越大。',
'固定外生增益下风险为R(g)=(1-g)^2 A+g^2 V，最优固定增益为A/(A+V)。观测能量增益仅在固定支持、各向同性Gaussian噪声、常lambda等条件下有Stein型统一改进界。当前数据依赖相干、支持选择、相关噪声及Q/子波偏差不满足全部条件，不能外推为全流程保证。',
f"理论数值核对全部通过={theory['all_checks_pass']}，每个Monte Carlo案例50,000次抽样，包含改善例与失败例。",
'两个反例也已验证：不同组增益可破坏原子抵消、使合成波形能量增加；固定a的观测能量映射最大局部导数为9/8，所以g<=1不等于非扩张。','',
'## 组合机制证据','',
'从固定Git提交的decomposition.csv和factorial.csv复算全部10个region=all场景的raw/coh/interaction与rhoV/rhoB。三种MSE/偏差平方/方差对比的算术一致性见analysis/archive_audit.json。该检查是对归档表的独立算术复算，不是重跑680次噪声实验。',
'该证据支持既有机制判断：相干收缩消除了QScan的一部分方差收益，QScan的偏差代价仍保留；不支持把两个模块的收益简单相加。','',
'## 创新性可以怎样主张','',
'可以主张具体构造的差异：由未收缩反射事件组的邻道相干调节条件传播噪声收缩，邻道不直接加权加入目标道波形；以及对Q选择和相干收缩替代关系的机制研究。',
'不能把通用Wiener公式、稀疏求解后去偏、Q扫描、源域结构约束或Stein恒等式本身当首创。没有直接波形相加也不等于没有跨道增益影响。',
'是否已经证明性能优越，必须以本报告的SC对照门槛结果为准。即使通过，也只涉及这一公开声明的独立实现、开发预算和已见模型，不能推广成优于全部文献或实测资料。','',
'## 仍未完成的证据','',
'- 2018多道bSBL、2025扫描Q和RE-BM3D的完整作者代码/全部细节对照；尚未取得，不以简化替代冒充。',
'- 新地质模型的一次性确认和二维/实测优越性验证。本轮未改动或接管正在进行的二维正演。','',
'## 可复现入口','',
'运行环境为捆绑Python 3.12、NumPy与任务cache/deps中的SciPy；cache不上传。依次运行scripts/run_stage.py theory、scripts/check_operators.py、scripts/audit_archives.py、scripts/run_stage.py dev、scripts/run_stage.py test、scripts/make_report.py。',
'日志记录实际命令、开始/结束、退出码和源码SHA256。首次理论JSON序列化故障已修复，计算未改；首次开发的无关脚本变化与后续扩大参数搜索均已明示。','']
(ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
print(json.dumps(dict(gates=gates,mean_ratio=ratio,baseline=summary[base],standard=summary['Standard']),ensure_ascii=False))
