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
s=json.loads((R/'analysis/summary.json').read_text(encoding='utf-8'));rows=s['metrics'];plan=json.loads((R/'PLAN.json').read_text(encoding='utf-8'))
names={'fixed':'最初分步法（共享 Q / 子波）','unpruned':'未剪枝重拟合','optimized':'优化剪枝法'}
def row(protocol,condition,method,error):return next(a for a in rows if a['protocol']==protocol and a['condition']==condition and a['method']==method and a['uniform'] and a['Q_error_percent']==error)
def fmt(v):return f'{v:.4f}'
lines=['# 合成资料 Q 误差稳定性试验','','## Material Passport','','任务：停止真实资料实验，回到原反射模型，检验补偿所用 Q 发生偏差时的稳定性。执行完成，旧生产入口保持不变。','','## 设计与误差定义','',
'主试验为已知四层真 Q=[40,60,85,110] 的反射模型，32 Hz 零相位 Ricker 子波，dt=4 ms、451 点、11 道、450 m 固定偏移距。九组反射事件含四组弱事件。两组观测分别有 14 dB 和 8 dB 的时空相关随机噪声，所有 Q 情况复用同一观测。',
'统一相对误差定义为 δ=(Q使用/Q真值)-1，所有层同乘 1+δ，测试 -50% 到 +50%、步长 10%。另测两种深度变化：[0.7,0.9,1.1,1.3] 与 [1.3,1.1,0.9,0.7]。没有重新估计 Q、子波、几何、噪声尺度或正则化。',
'复杂模型补充试验复用之前 Moderate / Strong 观测，预先取原第 1、11、21、30、40、50、60、70、79、89、99 道，测试 -50/-30/-10/0/+10/+30/+50%。其 Q 基准来自既有四层标定；这部分是相对标定 Q 的扰动，不是相对二维空间真 Q 的百分比误差。',
'所有方法使用相同的所施加 Q 与子波。最初分步法保留原两步固定岭回归公式；优化方法保留原 Ricker 反射率字典、支持集重拟合与分组剪枝，未使用上一轮自然地震 W=I 的适配。65 个核心 MATLAB 文件保持原哈希。','',
'## 已知真 Q 模型的整体误差','',
'| 条件 | 方法 | Q -50% | Q -20% | Q 正确 | Q +20% | Q +50% |','|---|---|---:|---:|---:|---:|---:|']
for cond in ['Moderate','Strong']:
    for method in names:lines.append('| '+cond+' | '+names[method]+' | '+' | '.join(fmt(row('primary',cond,method,e)['relative_error']) for e in [-50,-20,0,20,50])+' |')
lines+=['',
'相同百分比的 Q 低估与高估并不对应相同的吸收误差。传播时间不变时 t* 与 1/Q 成正比：Q 低估 50% 使吸收积分翻倍，Q 高估 50% 使吸收积分降为原来的 2/3。低估 Q 往往造成更强的逆补偿和噪声响应，需要结合实测曲线判断。',
'复杂模型显示更严重的失配风险：Strong 条件下，相对标定 Q 低估 50%，优化剪枝误差为 '+fmt(row('complex','Strong','optimized',-50)['relative_error'])+'，原分步法为 '+fmt(row('complex','Strong','fixed',-50)['relative_error'])+'；优化法在这一大失配下反而更差。',
'弱反射也不能仅归因于 Q：主模型 Strong、Q 正确时，剪枝版弱振幅误差为 '+fmt(row('primary','Strong','optimized',0)['weak_error'])+'，未剪枝重拟合为 '+fmt(row('primary','Strong','unpruned',0)['weak_error'])+'。剪枝降低整体噪声响应，同时仍会丢失部分弱事件。',
'相对敏感性采用每个方法自身零偏差结果作基准。某方法在零偏差时误差较高，曲线可能显得更平；因此绝对恢复误差、相对变化及弱事件指标同时展示。噪声和正则化也可能使最低误差点偏离 Q 真值，不能据此反推应该故意偏移 Q。','',
'## 弱反射与噪声响应定义','',
'主模型弱事件振幅：在真事件位置 ±40 ms 窗内，将输出投影到该事件 Ricker 模板，比较估计振幅与生成振幅。真事件时刻和幅度只用于合成设计与评分，不传入求解器。报告相对振幅误差与保留符号的振幅比中位数。复杂模型使用旧试验冻结的低能量区掩膜，不能把两类弱指标混为同一量。',
'噪声响应增益仅在主模型定义：||F_Q(D含噪)-F_Q(D纯净)|| / ||新增真实噪声||。两次求解都使用同一组从含噪观测估计并冻结的 sigma；非线性支持变化也计入这一响应，它不是假定的线性频率增益。另保存纯净输入失真以区分模型偏差与噪声影响。',
'主试验生成与字典使用同一分层吸收模型，属于受控模型误差与先验稳定性诊断，不是独立传播模拟器验证；复杂模型包含子波失配、横向变化、多次波、相干干扰及静校正等额外偏差，且只评价预先选定的 11 道。','',
'## 注册的实用容差检查','',
'相对自身零偏差结果，要求整体误差不超过 110%、弱振幅误差不超过 115%、噪声响应增益不超过 200%。三者同时通过，记录包含 0 的连续已测试区间。该容差是本轮预先选定的描述性规则，不是物理通用阈值或置信区间。到达 ±50% 边界仅表示没有在本测试网格中触及失败，不能外推。','',
'| 条件 | 方法 | 连续通过的已测试 Q 误差区间 |','|---|---|---|']
for g in s['practical_intervals']:
    lo,hi=g['contiguous_tested_interval_percent'];lines.append(f"| {g['condition']} | {names[g['method']]} | [{lo:+.0f}%, {hi:+.0f}%] |")
lines+=['','## 随深度变化的误差','',
'| 条件 | Q 因子（浅→深） | 方法 | 整体误差 | 弱振幅误差 |','|---|---|---|---:|---:|']
for a in rows:
    if a['protocol']=='primary' and not a['uniform']:
        pat='0.7/0.9/1.1/1.3' if a['scenario']=='shallow_low_deep_high' else '1.3/1.1/0.9/0.7'
        lines.append(f"| {a['condition']} | {pat} | {names[a['method']]} | {a['relative_error']:.4f} | {a['weak_error']:.4f} |")
lines+=['','## 复杂模型：相对标定 Q 的扰动','',
'| 条件 | 方法 | 标定 Q -50% | 标定 Q 正常 | 标定 Q +50% |','|---|---|---:|---:|---:|']
for cond in ['Moderate','Strong']:
    for method in names:lines.append('| '+cond+' | '+names[method]+' | '+' | '.join(fmt(row('complex',cond,method,e)['relative_error']) for e in [-50,0,50])+' |')
lines+=['',
'复杂模型没有重新生成或去噪，不让 Q 变化同时改变观测；其低能量掩膜也保持不变。这里的 0% 只是沿用原标定，并不表示所有横向位置的真实 Q 被准确恢复。','',
'## 完成检查与复现','',
f"主试验完成 {s['checks']['primary']['scenarios']} 组 Q 设置 / {s['checks']['primary']['dataset_cases']} 个含噪条件；复杂试验完成 {s['checks']['complex']['scenarios']} 组 / {s['checks']['complex']['dataset_cases']} 个条件。两次 MATLAB 退出码均为 0，优化标志均为 0。数值收敛与恢复精度分别报告，不能因为程序正常就认定大 Q 误差安全。", 
'首次主模型入口出现空结构体初始化错误，修正为预声明字段后复算；原失败日志保留。未改变科学参数或算法核心。',
'inputs/ 为原始观测生成要素、冻结标定和复杂子集；matlab_primary/ 与 matlab_complex/ 保存所有输出、纯净响应、诊断及进度；analysis/metrics.csv 为完整指标。PLAN.json 注册设计，FROZEN_CORE_HASHES.json 校验 65 个核心源文件。',
"MATLAB 入口为 run_q_stability('primary') 和 run_q_stability('complex')；已有同名输出目录时拒绝覆盖，复算前在副本中将目录改名存档。Python analyze.py 评分作图，build_report.py 生成报告与交付包。prepare.py 是首次准备脚本，引用工作区此前复杂试验，复算已有输入不需要再次运行。",
'本轮未读取或处理真实资料。结果覆盖一个已知四层模型及一个复杂模型子集、两个噪声条件，不是跨地质模型的通用稳定性证明。']
(R/'RESULTS_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
figs=[('真 Q 误差：整体、弱反射与噪声','primary_stability'),('相对于自身基线的敏感性','relative_sensitivity'),('弱反射振幅恢复','weak_amplitudes'),('四层模型：低估 / 正确 / 高估 Q','primary_gathers'),('复杂模型：相对标定 Q 的扰动','complex_stability'),('复杂模型：波形对比','complex_gathers')]
html='<!doctype html><meta charset="utf-8"><title>Q 误差稳定性</title><style>body{max-width:1450px;margin:30px auto;font:18px sans-serif;line-height:1.8;color:#16354b}img{width:100%}a{color:#087f8c}</style><h1>合成资料 Q 误差稳定性试验</h1><p>停止真实资料试验。已知真 Q 模型测试 −50% 到 +50% 的误差与两种深度误差；复杂模型仅扰动既有标定 Q，二者定义不同。</p><p><a href="q_stability_results.pdf">PDF 报告</a> | <a href="RESULTS_ZH.md">完整数值与方法</a></p>'
for title,name in figs:html+=f'<h2>{title}</h2><img src="figures/{name}.png" alt="{title}">'
(R/'index.html').write_text(html,encoding='utf-8')

pdfmetrics.registerFont(TTFont('CN','C:/Windows/Fonts/msyh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('CNB','C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0))
W,H=landscape(A4);M=38;writer=PdfWriter();style=ParagraphStyle('body',fontName='CN',fontSize=11,leading=19,wordWrap='CJK',textColor=colors.HexColor('#16354b'));small=ParagraphStyle('small',parent=style,fontSize=9,leading=15)
def page(title):
    b=io.BytesIO();c=canvas.Canvas(b,pagesize=(W,H));c.setFillColor(colors.HexColor('#087f8c'));c.rect(0,H-8,W,8,fill=1,stroke=0);c.setFont('CN',9);c.drawString(M,H-32,'2026-10-23 / 合成反射模型 Q 误差稳定性');c.setFillColor(colors.HexColor('#16354b'));c.setFont('CNB',22);c.drawString(M,H-66,title);c.setFont('CN',8);c.drawString(M,23,'同一观测与参数 | 65 个核心源文件不变 | 本轮无真实资料');c.drawRightString(W-M,23,str(len(writer.pages)+1));return b,c
def para(c,text,y,st=style):
    p=Paragraph(text,st);_,h=p.wrap(W-2*M,H);p.drawOn(c,M,y-h);return y-h-12
def save(b,c):c.save();b.seek(0);writer.add_page(PdfReader(b).pages[0])
def table(c,rr,y,widths):
    tab=Table([[Paragraph(str(v),small) for v in r] for r in rr],colWidths=widths);tab.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e1eff2')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f3f6f8')]),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7)]));_,h=tab.wrap(W-2*M,H);tab.drawOn(c,M,y-h);return y-h-15
b,c=page('Q 误差稳定性：绝对精度与敏感性');y=H-100
y=para(c,'回到原反射模型与 Ricker 字典，保持观测、几何、子波、噪声尺度和正则化不变，只改变补偿 Q。主模型四层真 Q=[40,60,85,110]；两个噪声条件为 14 dB / 8 dB。',y)
rr=[['条件 / 方法','Q -50%','Q 正确','Q +50%']]
for cond in ['Moderate','Strong']:
    for method in ['fixed','optimized']:rr.append([cond+' / '+names[method]]+[fmt(row('primary',cond,method,e)['relative_error']) for e in [-50,0,50]])
y=table(c,rr,y,[330,145,145,145])
y=para(c,'优化方法在 Q 正确时可有更低误差，但大幅低估 Q 时仍会严重恶化。相同的 Q 百分比高估 / 低估对应不对称的吸收误差，不能只看求解器是否收敛。',y)
y=para(c,'复杂模型 Strong、标定 Q 低估 50% 时：优化法误差 '+fmt(row('complex','Strong','optimized',-50)['relative_error'])+'，分步法 '+fmt(row('complex','Strong','fixed',-50)['relative_error'])+'，优化法反而更差。主模型 Q 正确时，剪枝弱振幅误差 '+fmt(row('primary','Strong','optimized',0)['weak_error'])+' 也高于未剪枝的 '+fmt(row('primary','Strong','unpruned',0)['weak_error'])+'；弱反射不足还受剪枝影响。',y)
save(b,c)
b,c=page('试验设计、评分与实用容差');y=H-99
for txt in ['主模型：11 道、dt=4 ms、451 点、450 m 偏移距、32 Hz 子波、九组反射事件及四组弱事件。统一 Q 误差 -50% 至 +50%，步长 10%；另测浅低深高与浅高深低两种模式。',
'最初分步法与优化方法共享相同 Q / 子波，保留原核心公式。弱振幅用真事件 ±40 ms 的 Ricker 模板投影评分；噪声响应为含噪输出减纯净输出，除以输入噪声范数，两次使用同一冻结 sigma。真事件信息不传给求解器。',
'容差同时要求：整体误差≤自身零偏差的 110%，弱振幅误差≤115%，噪声响应≤200%。下表只是包含 0 的连续通过网格区间，不是置信区间或通用安全范围。']:
    y=para(c,txt,y)
rr=[['条件','方法','连续通过的已测试 Q 误差区间']]
for g in s['practical_intervals']:
    lo,hi=g['contiguous_tested_interval_percent'];rr.append([g['condition'],names[g['method']],f'[{lo:+.0f}%, {hi:+.0f}%]'])
y=table(c,rr,y,[150,330,285]);save(b,c)
b,c=page('深度变化误差与复杂模型范围');y=H-97
rr=[['条件 / 模式（浅→深）','方法','整体误差','弱振幅误差']]
for a in rows:
    if a['protocol']=='primary' and not a['uniform']:
        pat='0.7/0.9/1.1/1.3' if a['scenario']=='shallow_low_deep_high' else '1.3/1.1/0.9/0.7'
        rr.append([a['condition']+' / '+pat,names[a['method']],fmt(a['relative_error']),fmt(a['weak_error'])])
y=table(c,rr,y,[250,275,120,120]);y=para(c,'复杂模型沿用原第 1、11、21、30、40、50、60、70、79、89、99 道。低能量区沿用旧冻结掩膜，不等同主模型的弱事件振幅。这里只改变原标定 Q，不宣称每道的空间真 Q 已知。',y,small);save(b,c)
for title,name in figs:
    b,c=page(title);para(c,'原反射模型；所有方法共享所施加 Q 和子波。真 Q 误差与相对标定 Q 扰动分别标注；复杂模型只评价预先固定的 11 道。',H-84,small);save(b,c)
    src=PdfReader(R/'figures'/f'{name}.pdf').pages[0];sw,sh=float(src.mediabox.width),float(src.mediabox.height);scale=min((W-2*M)/sw,(H-175)/sh)
    writer.pages[-1].merge_transformed_page(src,Transformation().scale(scale).translate((W-sw*scale)/2,46+(H-175-sh*scale)/2))
b,c=page('数值检查、复现与结论边界');y=H-100
for txt in [f"主模型 {s['checks']['primary']['dataset_cases']} 个含噪条件，复杂模型 {s['checks']['complex']['dataset_cases']} 个条件；两次 MATLAB 退出码均为 0，所有优化标志为 0；65 个核心源文件与原快照哈希一致。", 
'数值收敛不等于恢复稳定。Q 过低会增大吸收逆补偿；高噪声和先验失配下，最佳误差点也可能偏离 Q 真值，因此同时报告绝对误差与相对变化。',
'主试验生成与字典使用同一分层模型，评价受控 Q 偏差；复杂模型补充了多种模型失配，但仅为单一地质子集。结果不构成跨地质模型的通用稳定性证明。',
'首次主入口空结构体初始化错误已修正，原日志保留，科学参数未改。旧生产入口未修改，本轮未读取或处理真实资料。',
'交付包含 inputs/、65 个核心源文件、全部 MATLAB 输出、误差与弱振幅表、纯净输入响应、实用容差区间、图件和 PLAN.json。先看 index.html 或 RESULTS_ZH.md；复算已有输入无需再次运行准备脚本。']:
    y=para(c,txt,y)
save(b,c)
pdf=R/'q_stability_results.pdf';writer.add_metadata({'/Title':'合成资料Q误差稳定性试验','/Author':'研究实验记录'})
with pdf.open('wb') as f:writer.write(f)
tmp=R/'tmp/pdfs';tmp.mkdir(parents=True,exist_ok=True);doc=pdfium.PdfDocument(str(pdf));tiles=[]
for i in range(len(doc)):
    im=doc[i].render(scale=1.2).to_pil().convert('RGB');im.save(tmp/f'page_{i+1:02d}.png');im.thumbnail((560,400));tile=Image.new('RGB',(580,430),'#e6ebee');tile.paste(im,((580-im.width)//2,8));ImageDraw.Draw(tile).text((10,410),str(i+1),fill='black');tiles.append(tile)
sheet=Image.new('RGB',(1740,430*((len(tiles)+2)//3)),'white')
for i,tile in enumerate(tiles):sheet.paste(tile,((i%3)*580,(i//3)*430))
sheet.save(tmp/'contact_sheet.png');assert len(doc)==10
assert all(len(p.extract_text() or '')>80 for p in PdfReader(pdf).pages)
manifest={str(p.relative_to(R)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in R.rglob('*') if p.is_file() and not any(x in p.parts for x in ['tmp','__pycache__']) and p.suffix!='.log' and p.name not in ['MANIFEST.json','DELIVERY_CHECK.json']}
(R/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
dest=R.parent/'q_error_stability_20261023.zip'
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(R.rglob('*')):
        if p.is_file() and not any(x in p.parts for x in ['tmp','__pycache__']) and p.name not in ['build_report.log','DELIVERY_CHECK.json']:z.write(p,'q_error_stability_20261023/'+str(p.relative_to(R)).replace('\\','/'))
with zipfile.ZipFile(dest) as z:
    assert z.testzip() is None
    for name,h in manifest.items():assert hashlib.sha256(z.read('q_error_stability_20261023/'+name)).hexdigest()==h
check=dict(pdf=str(pdf),pages=len(doc),zip=str(dest),bytes=dest.stat().st_size,all_members_verified=True,completed_noisy_dataset_conditions=40,core_source_files_unchanged=65)
(R/'DELIVERY_CHECK.json').write_text(json.dumps(check,indent=2),encoding='utf-8');print(json.dumps(check,ensure_ascii=True))
