from pathlib import Path
import csv,hashlib,io,json,zipfile
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4,landscape
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from pypdf import PdfReader,PdfWriter,Transformation
import pypdfium2 as pdfium
from PIL import Image,ImageDraw
R=Path(__file__).resolve().parent
summary=json.loads((R/'analysis/summary.json').read_text(encoding='utf-8'))
meta=json.loads((R/'inputs/geometry_and_Q.json').read_text(encoding='utf-8'))
cape=json.loads((R/'analysis/cape_results.json').read_text(encoding='utf-8'))
controlled=summary['controlled_summary'];sens=list(csv.DictReader((R/'analysis/sensitivity_metrics.csv').open(encoding='utf-8-sig')))
def result(noise,method):return next(a for a in controlled if a['split']=='evaluation' and a['noise_factor']==noise and a['method']==method)['mean_relative_error']
labels={'input':'不补偿输入','fixed':'两步岭回归适配版','refit':'稀疏重拟合适配版','optimized':'稀疏剪枝适配版'}
lines=['# 有外部 Q 参考的真实资料补偿试验','','## Material Passport','','类型：自然地震单路径补偿适配、受控额外衰减一致性试验、井中 DAS 条件频谱检查。状态：执行完成；原生产流程未替换。','','## 核心结果','',
'Sumatra：4 条真实 P 波记录，来自 2 个地震。作者 α=0 的 Q 与 t* 均已匹配，SAC P 到时除以 Q 与作者 t* 一致。LWLI 只用于确定一个聚类频率，其余 3 条为评估记录，但同属另一个地震的不同台站，不是 3 个独立地震。',
'完成 12 组全路径补偿（4 条 × t* 的 0.75/1/1.25 倍）和 24 组额外衰减回放（4 条 × 3 个额外 t* × 2 个噪声条件）。所有核心优化标志为 0。','',
'### 评估记录受控回放误差','',
'| 额外噪声 | 不补偿输入 | 两步岭回归适配版 | 未剪枝重拟合 | 剪枝适配版 |','|---|---:|---:|---:|---:|']
for n in [0,1]:lines.append('| '+str(n)+'σ | '+' | '.join(f'{result(n,m):.4f}' for m in ['input','fixed','refit','optimized'])+' |')
lines+=['',
'每格是 3 条评估记录 × 3 个额外 t* 的 9 个情况的误差算术平均，重复情况有相关性，不给独立样本显著性解释。误差参考是原始已衰减真实记录，不是真实未衰减信号。',
'本次适配中，两步岭回归的额外衰减恢复误差更低；稀疏与剪枝对连续、复杂 P 波引入较大失真。不能将先前合成反射模型上的优化结论直接推广到自然地震。也不能据此判定原版反射算法在其原目标场景中优劣。','',
'## 适配内容与不变部分','',
'原反射模型要求分层几何和 Ricker 子波。这些自然地震没有对应的反射层几何和独立源子波，因而没有伪造层厚、速度或 Ricker 源。改用各条记录作者的常 Q t* 构造平稳单路径算子，偏移距分步为单位阵。',
'7 个核心 MATLAB 文件逐字节保持不变：固定两步岭回归、岭回归求解、吸收矩阵、稀疏求解、支持集重拟合、分组剪枝和 JSON 输出。新增试验入口提供字典 W=I，将未知量改为带限源波形系数；这不是原版 Ricker 反射率字典，故名称全部标明“适配版”。',
'吸收算子使用原代码的傅里叶约定：正频率 H(f)=exp[-π f t* + i 2 f ln(f/fref) t*]，fref=1 Hz。作者表只约束衰减参考；这里的相位模型是计算假设，没有相位真值验证。固定岭参数仍为 1e-3；重拟合与剪枝使用原自适应噪声规则，没有按回放误差调参。',
'原始资料已做 0.2 Hz 高通。本轮统一使用 0.2-0.4 Hz 和 6.5-8 Hz 平滑边缘的带通，dt=0.05 s；取 P 前 4 s 至后 12 s，评分 P 前 1 s 至后 7 s。背景取 P 前 24 s 至前 8 s，用 MAD 估计尺度。求解时作标量数值缩放，绘图恢复原文件幅度；无 AGC 或逐道显示归一化。','',
'## 同算子回放的边界','',
'对原真实波形 x 再施加已知额外吸收算子 G，生成 Gx + 新增有色噪声，然后恢复到 x。额外 t*=0.03/0.06/0.12 s；额外噪声为 0 或背景噪声标准差 1 倍，固定随机种子 20261023。',
'生成额外衰减与逆解使用同一吸收模型，属于匹配模型的一致性与先验失真检查，不是独立传播模型验证。没有新增噪声也不代表无原有实测噪声；回放目标仍包含真实原记录中的自然衰减、仪器及噪声。','',
'## 全路径补偿：可计算，但不能以真值误差评分','',
'| 台站 | 作者 Q | 作者 t* / s | 岭回归幅度增益 | 剪枝幅度增益 | 岭回归重放残差 | 剪枝重放残差 |','|---|---:|---:|---:|---:|---:|---:|']
for a in meta:
    rr=[r for r in summary['nominal_field_metrics'] if r['station']==a['station']];b=next(r for r in rr if r['method']=='fixed');p=next(r for r in rr if r['method']=='optimized')
    lines.append(f"| {a['station']} | {a['Q']:.2f} | {a['tstar']:.6f} | {b['amplitude_gain']:.2f} | {p['amplitude_gain']:.2f} | {b['replay_relative_residual']:.3f} | {p['replay_relative_residual']:.3f} |")
lines+=['',
'幅度增益是输出 / 输入在评分窗内的 L2 范数比。重放残差是将补偿输出通过假设吸收算子再传播后与输入的差异；小残差仅说明模型拟合，不是补偿真值。',
'背景响应为将单独背景窗用同一求解器处理的范数变化。稀疏求解器是非线性的，单独噪声输出为 0 不能解释成真实 P 波中的噪声已全部去除或 SNR 无穷大。','',
'## Cape 井中 DAS：两事件条件频谱试验','',
'使用 med_Q_chs 与 MAD 和配套测井速度，取各事件有效测量通道深度的 25%/75% 分位附近各 ±40 m 的窗，对振幅谱取中位数。只用实测估计范围内的 Q，区间内插值，不使用包含外推假定的 Q_interp。',
'为了构造条件检验，假设事件在下方、波沿竖直路径向上传播，计算 Δt*=∫dz/(Vp Q) 并对较浅窗乘 exp(πfΔt*)，增益上限 10 倍，检查 20-200 Hz 两窗均 SNR≥3 的频点。未独立校正 DAS 指向性 / gauge-length 响应、未核实精确射线路径，因此不能把谱比变平作为 Q 真值验证。Q±MAD 只用于压力测试，不是置信区间。',
f"事件一：{cape[0]['usable_bins']} 个可用频点，原始谱比斜率 {cape[0]['raw_ratio_slope_per_hz']:.6f} /Hz，重加权后 {cape[0]['corrected_ratio_slope_per_hz']:.6f} /Hz，反向变斜。简化路径与响应条件不充分，不能认定成功。", 
f"事件二：{cape[1]['usable_bins']} 个频点同时达到两窗 SNR≥3，没有可靠的谱比拟合。只展示条件增益，不报告恢复成功。",'',
'## 数值检查与复现','',
'MATLAB 退出码 0，36 组完成，优化标志 0；7 个核心源文件哈希与原快照一致。首次试验入口的空结构体追加发生容器问题，后改为 cell 容器，原错误日志保留，科学参数未改。',
'运行顺序：prepare_trial.py → MATLAB run_trial → cape_spectral_trial.py → analyze_trial.py → build_report.py。运行 MATLAB 时在本目录调用 run_trial。依赖 NumPy、Matplotlib、MATLAB；报告依赖 ReportLab、pypdf、pypdfium2、Pillow。',
'所有实际结果、背景响应、参数敏感性和回放目标在 matlab/*.json；CSV 在 analysis/；原始 SAC、作者结果表、Cape Q 及两个频谱样本在 raw/。报告与代码没有覆盖之前的试验。','',
'## 数据来源','',
'Sumatra 波形与 Q/t* 表：https://zenodo.org/records/6932666',
'Cape Q、速度、密度及频谱：https://gdr.openei.org/submissions/1823',
'Cape 方法与 DAS 响应讨论：https://doi.org/10.1029/2026JB034276',
'数据许可：两套原数据均为 CC BY 4.0；保留原作者、来源和修改说明。无真实无衰减参考，无生产流程推广。']
(R/'RESULTS_ZH.md').write_text('\n'.join(lines),encoding='utf-8')

figs=[('真实 P 波：LWLI','LWLI_field_trial'),('真实 P 波：MNSI','MNSI_field_trial'),('真实 P 波：SISI','SISI_field_trial'),('真实 P 波：MLSI','MLSI_field_trial'),('受控额外衰减：全部结果','controlled_recovery'),('受控回放：预先指定示例','controlled_example')]
figs += [('Cape 事件一：条件频谱检查','spec_D-1OB_P_arr_20240218145151_891022_cape_trial'),('Cape 事件二：SNR 不足','spec_D-1OB_P_arr_20240228125612_509795_cape_trial')]
html='<!doctype html><meta charset="utf-8"><title>有 Q 参考的真实补偿试验</title><style>body{max-width:1400px;margin:30px auto;font:18px sans-serif;line-height:1.75;color:#16354b}img{width:100%}a{color:#087f8c}</style><h1>有外部 Q 参考的真实补偿试验</h1><p>36 组完成。自然地震适配试验中，两步岭回归回放误差低于稀疏剪枝；Cape 简化条件未支持成功验证。实测无衰减真值缺失，回放只逆转额外施加的衰减。</p><p><a href="q_supported_trial_results.pdf">PDF 报告</a> | <a href="RESULTS_ZH.md">完整方法与限制</a></p>'
for title,name in figs:html+=f'<h2>{title}</h2><img src="figures/{name}.png" alt="{title}">'
(R/'index.html').write_text(html,encoding='utf-8')

pdfmetrics.registerFont(TTFont('CN','C:/Windows/Fonts/msyh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('CNB','C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0))
W,H=landscape(A4);M=38;writer=PdfWriter()
style=ParagraphStyle('body',fontName='CN',fontSize=11,leading=19,wordWrap='CJK',textColor=colors.HexColor('#16354b'));small=ParagraphStyle('small',parent=style,fontSize=9,leading=15)
def page(title):
    b=io.BytesIO();c=canvas.Canvas(b,pagesize=(W,H));c.setFillColor(colors.HexColor('#087f8c'));c.rect(0,H-8,W,8,fill=1,stroke=0);c.setFont('CN',9);c.drawString(M,H-32,'2026-10-23 / 有外部 Q 参考的真实资料试验');c.setFillColor(colors.HexColor('#16354b'));c.setFont('CNB',22);c.drawString(M,H-66,title);c.setFont('CN',8);c.drawString(M,23,'真实资料适配版 | 作者 Q 为参考估计 | 回放仅逆转新增衰减');c.drawRightString(W-M,23,str(len(writer.pages)+1));return b,c
def para(c,txt,y,s=style):
    p=Paragraph(txt,s);_,h=p.wrap(W-2*M,H);p.drawOn(c,M,y-h);return y-h-12
def save(b,c):c.save();b.seek(0);writer.add_page(PdfReader(b).pages[0])
def table(c,rows,y,widths):
    t=Table([[Paragraph(str(a),small) for a in rr] for rr in rows],colWidths=widths);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e1eff2')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f4f7f8')]),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8)]));_,h=t.wrap(W-2*M,H);t.drawOn(c,M,y-h);return y-h-16
b,c=page('真实资料补偿试验：结果概览');y=H-100
y=para(c,'已运行 4 条真实 P 波记录的全路径补偿与参数敏感性，以及已知额外衰减的受控回放，共 36 组。采用作者发布的 Q/t*，没有重新拟合这些参考值。',y)
rr=[['评估记录平均相对误差','不补偿输入','两步岭回归适配','未剪枝重拟合','剪枝适配']]
for n in [0,1]:rr.append([f'新增噪声 {n}σ']+[f'{result(n,m):.4f}' for m in ['input','fixed','refit','optimized']])
y=table(c,rr,y,[230,125,145,135,130])
for txt in ['本次适配中，两步岭回归的受控恢复误差更低，稀疏与剪枝会损失连续 P 波中的部分信息。不能把合成反射模型中的优化优势直接推广到自然地震。',
'参考目标是原真实记录，本身仍有自然衰减和噪声。额外衰减与逆解用同一算子，这是匹配模型回放，不是恢复真实无衰减波形的证明。',
'Cape 只做振幅谱检查：一个事件重加权后谱比反向变斜，另一个事件无足够信噪比频点。现有简化路径和响应条件不支持成功结论。']:
    y=para(c,txt,y)
save(b,c)
b,c=page('方法适配、输入条件与解释范围');y=H-98
for txt in ['4 条记录来自 2 个地震。LWLI 用于确定统一聚类频率，其余 3 个台站属于另一个地震。参数未按评估误差调整；重复条件与同事件台站不作为独立统计样本。',
'7 个核心 MATLAB 文件保持原哈希：两步岭回归、吸收算子、稀疏求解、支持重拟合及分组剪枝等。新增入口使用单路径 t*；偏移阶段为单位阵；W=I 估计未知源波形，不再假设原分层 Ricker 反射率。',
'作者 α=0 的 Q、t* 与 SAC P 到时逐条匹配。dt=0.05 s，统一带通边缘 0.2/0.4/6.5/8 Hz；P 前 4 s 到后 12 s 求解，P 前 1 s 到后 7 s 评分。背景为 P 前 24 s 到前 8 s，MAD 估计噪声尺度。',
'固定岭参数仍为 1e-3，稀疏及剪枝保留原规则。吸收相位采用原代码模型，fref=1 Hz，属于计算假设；无相位真值标定。标量数值预条件化后恢复原幅度，没有 AGC 或逐道显示归一化。',
'额外 t*=0.03/0.06/0.12 s，新增有色噪声为 0 或 1σ。无新增噪声时，目标仍含原实测噪声。原始全路径补偿只报告幅度、重放拟合及敏感性，不报告未知真值误差。']:
    y=para(c,txt,y)
save(b,c)
b,c=page('全路径补偿与诊断限制');y=H-99
rr=[['台站','Q / t*（s）','岭回归增益','剪枝增益','岭回归重放残差','剪枝重放残差']]
for a in meta:
    ff=next(z for z in summary['nominal_field_metrics'] if z['station']==a['station'] and z['method']=='fixed');pp=next(z for z in summary['nominal_field_metrics'] if z['station']==a['station'] and z['method']=='optimized')
    rr.append([a['station'],f"{a['Q']:.2f} / {a['tstar']:.6f}",f"{ff['amplitude_gain']:.2f}",f"{pp['amplitude_gain']:.2f}",f"{ff['replay_relative_residual']:.3f}",f"{pp['replay_relative_residual']:.3f}"])
y=table(c,rr,y,[80,185,120,120,130,130])
for txt in ['增益是评分窗内输出 / 输入的范数比。重放残差是对补偿结果再次施加假设吸收算子后与观测的差异，不能当作无衰减恢复误差。',
'背景窗也独立通过每个求解器。稀疏求解是非线性的，背景单独输出为 0 不表示真实信号中噪声全部消失，不能据此计算无穷大 SNR。',
'Cape 使用有实测估计的 Q、MAD 与速度，区间内插值，假设竖直向上传播；未独立校正 DAS 指向性和 gauge-length 响应。事件一谱比斜率 -0.000895 → +0.009657 /Hz；事件二两窗均 SNR≥3 的频点数为 0。均不能认定验证成功。']:
    y=para(c,txt,y)
save(b,c)
for title,name in figs:
    b,c=page(title);para(c,'自然地震单路径适配与条件频谱试验；原文件幅度标定。没有真实无衰减参考，不能根据高频增强或振幅增加判定恢复成功。',H-84,small);save(b,c)
    src=PdfReader(R/'figures'/f'{name}.pdf').pages[0];sw,sh=float(src.mediabox.width),float(src.mediabox.height);s=min((W-2*M)/sw,(H-175)/sh)
    writer.pages[-1].merge_transformed_page(src,Transformation().scale(s).translate((W-sw*s)/2,46+(H-175-sh*s)/2))
b,c=page('数据来源、检查与后续判断');y=H-99
for txt in ['MATLAB 退出码 0，36 组全部完成，优化标志 0；7 个核心源文件哈希一致。首次入口的结构体容器问题已修正为 cell 容器，原失败日志保留，科学参数未变。',
'当前建议：自然地震单路径试验保留较简单的两步岭回归作为基线；稀疏剪枝不直接推广。反射炮集的毕设主验证仍需要配套 VSP / 井资料和真实传播条件。',
'完整包包含原始 SAC、作者 Q 表、Cape Q / 频谱样本、冻结参数、MATLAB 核心快照、所有输出和表格。RESULTS_ZH.md 给出完整定义，index.html 为图集。']:
    y=para(c,txt,y)
for label,url in [('Sumatra 波形与 Q/t*','https://zenodo.org/records/6932666'),('Cape Q 及频谱','https://gdr.openei.org/submissions/1823'),('Cape 衰减与 DAS 响应研究','https://doi.org/10.1029/2026JB034276')]:y=para(c,f'{label}：<link href="{url}" color="#087f8c">{url}</link>',y,small)
y=para(c,'原始数据为 CC BY 4.0；引用原作者及数据集。作者 Q 是参考估计，不是介质精确真值。',y,small);save(b,c)
pdf=R/'q_supported_trial_results.pdf';writer.add_metadata({'/Title':'有外部Q参考的真实资料补偿试验','/Author':'研究实验记录'})
with pdf.open('wb') as f:writer.write(f)
tmp=R/'tmp/pdfs';tmp.mkdir(parents=True,exist_ok=True);doc=pdfium.PdfDocument(str(pdf));tiles=[]
for i in range(len(doc)):
    im=doc[i].render(scale=1.2).to_pil().convert('RGB');im.save(tmp/f'page_{i+1:02d}.png');im.thumbnail((560,400));tile=Image.new('RGB',(580,430),'#e6ebee');tile.paste(im,((580-im.width)//2,8));ImageDraw.Draw(tile).text((10,410),str(i+1),fill='black');tiles.append(tile)
sheet=Image.new('RGB',(1740,430*((len(tiles)+2)//3)),'white')
for i,tile in enumerate(tiles):sheet.paste(tile,((i%3)*580,(i//3)*430))
sheet.save(tmp/'contact_sheet.png');assert len(doc)==12
assert all(len(p.extract_text() or '')>80 for p in PdfReader(pdf).pages)
manifest={str(p.relative_to(R)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in R.rglob('*') if p.is_file() and not any(x in p.parts for x in ['tmp','__pycache__']) and p.suffix!='.log' and p.name not in ['MANIFEST.json','DELIVERY_CHECK.json']}
(R/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
dest=R.parent/'q_supported_trial_20261023.zip'
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(R.rglob('*')):
        if p.is_file() and not any(x in p.parts for x in ['tmp','__pycache__']) and p.name not in ['report_build.log','DELIVERY_CHECK.json']:z.write(p,'q_supported_trial_20261023/'+str(p.relative_to(R)).replace('\\','/'))
with zipfile.ZipFile(dest) as z:
    assert z.testzip() is None
    for name,h in manifest.items():assert hashlib.sha256(z.read('q_supported_trial_20261023/'+name)).hexdigest()==h
check=dict(pdf=str(pdf),pages=len(doc),zip=str(dest),bytes=dest.stat().st_size,all_members_verified=True,matlab_exit_code=0,completed_cases=36)
(R/'DELIVERY_CHECK.json').write_text(json.dumps(check,indent=2),encoding='utf-8');print(json.dumps(check,ensure_ascii=True))
