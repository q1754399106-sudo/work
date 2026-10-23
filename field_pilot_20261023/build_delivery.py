from pathlib import Path
import csv, hashlib, io, json, shutil, zipfile
from statistics import mean
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
res=json.loads((R/'analysis/results.json').read_text(encoding='utf-8'))
frozen=json.loads((R/'FROZEN_PARAMETERS.json').read_text(encoding='utf-8'))
rows=res['metrics'];held=[r for r in rows if r['split']=='held_out']
q=[r for r in res['apparent_Q_summary'] if r['stage']=='local' and r['record'] in [115,116,117,118]]
mx=max(c['relative_change'] for c in res['probe_controls'])
lines=['# Alaska 多炮局部去干扰与 Q 可估计性试验','','## Material Passport','','类型：真实资料方法试验；执行状态：已完成；评价类型：描述性处理指标和有限人工反射控制。无无噪声实测真值。','','## 核心结果','',
f"4 个开发炮（111-114）、4 个留出炮（115-118）。留出炮位于同一测线，不能等同独立地区泛化。冻结 f90={frozen['observed_f90_hz']:.3f} Hz，减除强度={frozen['subtraction_strength']:.2f}。参数选择后才解析留出波形；修正末尾数值检查后按完全相同参数复算，首次失败日志保留。",
f"留出炮减去输入范数范围 {min(r['removed_norm_fraction'] for r in held):.2%}-{max(r['removed_norm_fraction'] for r in held):.2%}；慢速低频谱锥功率平均降低 {mean(r['slow_cone_power_reduction'] for r in held):.2%}。这些是处理量指标，不是 SNR 或保幅改善。",
f"27 个人工反射探针（RMS 速度 2000/3000/4500 m/s，零偏移时 0.5/1.0/1.8 s，频率 20/35/55 Hz）最大相对变化 {mx:.2%}。探针数量不等于独立样本数，也不能替代真实反射保真检查。",'',
'## 逐炮结果','','| 记录号 | 划分 | 减去范数 | 慢速谱锥功率降低 | 晚时近偏移区变化 | 原全局基线减去范数 |','|---|---|---:|---:|---:|---:|']
for r in rows:lines.append(f"| {r['record']} | {r['split']} | {r['removed_norm_fraction']:.2%} | {r['slow_cone_power_reduction']:.2%} | {r['late_near_offset_relative_change']:.2%} | {r['baseline_removed_norm_fraction']:.2%} |")
lines+=['','晚时近偏移区定义为 1.5-3.5 s、|偏移距|<800 m；它混合反射与噪声，变化量只作可能泄漏提示。|偏移距|<250 m 的振幅保持不变是显式保护门控的结果，不是自动证明保真。','',
'## 局部方法','',
'左右排列分别计算，512 点 × 24 道局部窗，步长 128 点 × 6 道，平方根 Hann 窗加权重叠。慢速谱锥在 1000-2200 m/s 间平滑过渡，频率在开发炮 f90 的 0.3-0.7 倍间平滑过渡；排除最慢目标速度的空间混叠频率以上部分。仅在声明的慢速到时范围和 3 s 以内减除，中央偏移距保护。',
'强度从 0.20/0.35/0.50 中按开发炮慢速谱锥功率降低选择，并要求所有人工反射探针变化不超过 2%。使用与原全局软 f-k 基线相同观测对照。未降低原长孔径候选门槛冒充原方法成功。','',
'## 表观 Q 检查','',
'在开发炮平均双曲线相似性上自动选取两个相隔至少 0.35 s 的事件假设，随后固定时间与 RMS 速度，在留出炮原始幅度上直接截取事件谱。未对频谱输入做 AGC、道归一化或 NMO 拉伸。每个事件使用 0.192 s 窗及三个正弦窗功率平均，在 18-55 Hz 之间拟合 ln(A晚/A早)=a+b f；仅 b<0 时计算表观 Q=-πΔt/b。它依赖常 Q、可比反射子波及频率无关幅度比等假设。',
'两事件为自动假设，尚未确认对应同一地质层位；不同反射的反射系数、传播路径和耦合差异仍会污染结果。窗口前后偏移 16 ms 用于敏感性检查。有限频率点、同炮多道及相邻炮具有相关性，本轮不报告虚假的独立样本置信区间。','',
'| 留出记录 | 负谱斜率比例 | 负斜率且 R²≥0.70 | R² 中位数 |','|---|---:|---:|---:|']
for z in q:lines.append(f"| {z['record']} | {z['negative_slope_fraction']:.1%} | {z['high_r2_negative_fraction']:.1%} | {z['median_r2']:.3f} |")
lines+=['','Q 可估计性门槛：未通过。理由：']+['- '+s for s in res['q_gate']['reasons']]+['',
'因此没有填入任意 Q，也没有运行两种方法的实测 Q 补偿对比。完成的是多炮去干扰及 Q 可估计性检查。后续需要确认反射事件、核对几何 / 静校正及速度，或引入 VSP / 井资料，再验证估计与补偿。','',
'## 文件与复现','',
'run_field_pilot.py 执行实验；build_delivery.py 生成报告。所有原始与输出均保留。data/record_*.npz 包含 original、local、removed、baseline 等数组。analysis/ 保存逐道谱比、所有 Q 候选、参数敏感性、速度相似性扫描和保真控制表。PLAN.json 和 FROZEN_PARAMETERS.json 保留首次注册及冻结内容。',
'完整 ZIP 包额外携带原始 SEG-Y 和基线脚本，解压后可使用 NumPy、Matplotlib 运行 run_field_pilot.py。重新运行将校验已冻结参数，不接受留出结果驱动的参数改变。',
'','## 来源','',
'USGS 调查目录：https://walrus.wr.usgs.gov/namss/survey/n-pr-81-ak/',
'几何与处理脚本：https://ahay.org/RSF/book/data/alaska/line31-81-m8r.html',
'原始下载：https://ahay.org/data/alaska/L23535.SGY',
'光谱比衰减估计及路径问题：https://doi.org/10.1093/jge/gxab025',
'原始文件 SHA256：'+res['source_sha256'],
'几何声明：偏移距由发布者 -5225 ft 起点、110 ft 道距重建；原始源检坐标与偏移距头均为 0。未使用这些缺失的坐标假装建立可信 CMP 或射线路径模型。']
(R/'RESULTS_ZH.md').write_text('\n'.join(lines),encoding='utf-8')

figs=[('开发炮 111','record_111'),('开发炮 114','record_114'),('留出炮 115','record_115'),('留出炮 118','record_118'),('处理量与人工反射控制','metrics'),('Q 可估计性','q_feasibility')]
html='''<!doctype html><meta charset="utf-8"><title>真实多炮试验</title><style>body{max-width:1400px;margin:30px auto;font:18px sans-serif;color:#16354b}img{width:100%}p{line-height:1.8}a{color:#087f8c}</style><h1>Alaska 多炮局部去干扰与 Q 可估计性试验</h1>'''
html+=f'<p>开发 111-114，留出 115-118。留出炮慢速谱锥功率平均降低 {mean(r["slow_cone_power_reduction"] for r in held):.2%}，不是实测 SNR 改善。Q 门槛未通过，本轮未运行实测 Q 补偿。</p><p><a href="field_pilot_results.pdf">下载 PDF 报告</a> | <a href="RESULTS_ZH.md">完整文字与方法记录</a></p>'
for title,file in figs:html+=f'<h2>{title}</h2><img src="figures/{file}.png" alt="{title}对比图">'
(R/'index.html').write_text(html,encoding='utf-8')

pdfmetrics.registerFont(TTFont('CN','C:/Windows/Fonts/msyh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('CNB','C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0))
W,H=landscape(A4);M=38;writer=PdfWriter();NAVY=colors.HexColor('#16354b')
style=ParagraphStyle('body',fontName='CN',fontSize=11,leading=19,wordWrap='CJK',textColor=NAVY)
small=ParagraphStyle('small',parent=style,fontSize=9,leading=15)
def page(title):
    b=io.BytesIO();c=canvas.Canvas(b,pagesize=(W,H));c.setFillColor(colors.HexColor('#087f8c'));c.rect(0,H-8,W,8,fill=1,stroke=0);c.setFont('CN',9);c.drawString(M,H-32,'2026-10-23 / 真实资料方法试验');c.setFillColor(NAVY);c.setFont('CNB',22);c.drawString(M,H-66,title);c.setFont('CN',8);c.drawString(M,23,'Alaska 31-81 | 同一测线相邻炮，非独立地区验证');c.drawRightString(W-M,23,str(len(writer.pages)+1));return b,c
def para(c,txt,y,s=style):
    p=Paragraph(txt,s);_,h=p.wrap(W-2*M,H);p.drawOn(c,M,y-h);return y-h-12
def save(b,c):c.save();b.seek(0);writer.add_page(PdfReader(b).pages[0])
def table(c,rr,y,widths):
    z=Table([[Paragraph(str(a),small) for a in r] for r in rr],colWidths=widths)
    z.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e1eff2')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f3f6f8')]),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7)]));_,h=z.wrap(W-2*M,H);z.drawOn(c,M,y-h);return y-h-15

b,c=page('多炮局部去干扰与 Q 可估计性');y=H-99
for txt in [lines[8],lines[9],lines[10],
'局部候选确实改变了真实炮集，减除量明显小于之前的全局基线；这不等同于恢复更准确。减去分量仍可能包含有效波，人工探针只覆盖有限反射族。',
'Q 可估计性门槛未通过。本轮不使用任意 Q 做补偿演示，尚未形成两种 Q 补偿方法的实测效果对比。']:
    y=para(c,txt,y)
save(b,c)
b,c=page('逐炮处理量与保真提示');y=H-96
rr=[['记录 / 划分','减去输入范数','慢速谱锥功率降低','晚时近偏移变化','原全局基线减去范数']]
for r in rows:rr.append([str(r['record'])+(' / 开发' if r['split']=='development' else ' / 留出'),f"{r['removed_norm_fraction']:.2%}",f"{r['slow_cone_power_reduction']:.2%}",f"{r['late_near_offset_relative_change']:.2%}",f"{r['baseline_removed_norm_fraction']:.2%}"])
y=table(c,rr,y,[160,140,155,155,155]);y=para(c,'慢速谱锥功率变化只衡量目标谱区被压制多少，不是 SNR。晚时近偏移变化用于提示可能的有效波泄漏。中央 250 m 不变来自显式门控保护；没有实测无噪声参考。',y,small)
y=para(c,'参数在开发炮确定，留出炮波形随后才解析。首次执行末尾的绝对数值容差不适应原始幅度量级，改为相对容差后按原冻结参数复算，方法和参数未改变；失败记录完整保留。',y,small);save(b,c)
b,c=page('事件假设与表观 Q 诊断');y=H-98
for e in res['event_hypotheses']:y=para(c,f"事件假设：t0={e['t0_s']:.3f} s，RMS 速度={e['rms_velocity_m_s']:.0f} m/s，开发炮平均相似性={e['development_mean_semblance']:.3f}。",y)
y=para(c,'在每道的假设到时直接截取 0.192 s 波形，三正弦窗估计振幅谱，在 18-55 Hz 拟合晚 / 早事件对数谱比。未做 NMO 拉伸，未将显示增益或逐道归一化用于谱比。窗口偏移 ±16 ms 检验敏感性。',y,small)
rr=[['留出记录','负谱斜率比例','负斜率且 R²≥0.70','R² 中位数']]+[[a['record'],f"{a['negative_slope_fraction']:.1%}",f"{a['high_r2_negative_fraction']:.1%}",f"{a['median_r2']:.3f}"] for a in q]
y=table(c,rr,y,[170,190,220,185]);y=para(c,'自动事件不等于确认层位。不同反射系数、路径、耦合及噪声可导致频谱差异；拟合得出的表观 Q 不能直接作为真实介质 Q。所有候选值及拒绝原因均保存在 analysis/。',y,small)
y=para(c,'本轮诊断：负斜率比例、拟合质量、窗口扰动稳定性及跨炮离散度按注册门槛检查；独立事件解释与传播标定尚不满足，因此不进入实测 Q 补偿。',y,small);save(b,c)
for title,file in figs:
    b,c=page(title);para(c,'共同色标与统一显示增益；真实炮集没有纯净参考。图中减去分量用于检查泄漏，不能直接当作噪声真值。',H-83,small);save(b,c)
    src=PdfReader(R/'figures'/f'{file}.pdf').pages[0];sw,sh=float(src.mediabox.width),float(src.mediabox.height);s=min((W-2*M)/sw,(H-175)/sh)
    writer.pages[-1].merge_transformed_page(src,Transformation().scale(s).translate((W-sw*s)/2,46+(H-175-sh*s)/2))
b,c=page('方法、数据来源与交付文件');y=H-99
for txt in [
'左右排列分开做局部 f-k，512 点 × 24 道重叠窗，时间 / 空间步长 128 点 / 6 道。慢速掩膜在 1000-2200 m/s 间平滑过渡，频带与开发炮 f90 关联；3 s 以后和中央偏移距受保护。',
'原始源检坐标和偏移距头为零；按发布者脚本重建偏移距。未据此伪造 CMP 或精确射线路径。未覆盖原始文件，数据守恒采用按输入幅度归一化的最大残差检查。',
'后续工作：确认反射事件与几何 / 静校正，约束速度和传播模型，再引入 Q 估计及两种补偿方法的对照。当前记录仅支持有限场景下的处理量和探针响应结论。']:
    y=para(c,txt,y)
for label,url in [('USGS 原始调查','https://walrus.wr.usgs.gov/namss/survey/n-pr-81-ak/'),('几何与处理脚本','https://ahay.org/RSF/book/data/alaska/line31-81-m8r.html'),('原始 SEG-Y 下载','https://ahay.org/data/alaska/L23535.SGY'),('光谱比衰减估计及路径问题','https://doi.org/10.1093/jge/gxab025')]:y=para(c,f'{label}：<link href="{url}" color="#087f8c">{url}</link>',y,small)
y=para(c,'完整包包含原始 SEG-Y、八炮输入与输出、谱比与控制表、冻结参数、程序、首次校验失败记录、全部图件及报告。RESULTS_ZH.md 提供方法细节；index.html 为图集。',y,small);save(b,c)
pdf=R/'field_pilot_results.pdf';writer.add_metadata({'/Title':'真实多炮局部去干扰与Q可估计性试验','/Author':'研究实验记录'})
with pdf.open('wb') as f:writer.write(f)
tmp=R/'tmp/pdfs';tmp.mkdir(parents=True,exist_ok=True);doc=pdfium.PdfDocument(str(pdf));ims=[]
for i in range(len(doc)):
    im=doc[i].render(scale=1.2).to_pil().convert('RGB');im.save(tmp/f'page_{i+1:02d}.png');im.thumbnail((560,400));tile=Image.new('RGB',(580,430),'#e6ebee');tile.paste(im,((580-im.width)//2,8));ImageDraw.Draw(tile).text((10,410),str(i+1),fill='black');ims.append(tile)
contact=Image.new('RGB',(1740,430*((len(ims)+2)//3)),'white')
for i,im in enumerate(ims):contact.paste(im,((i%3)*580,(i//3)*430))
contact.save(tmp/'contact_sheet.png')
assert len(PdfReader(pdf).pages)==10
assert all(len(p.extract_text() or '')>80 for p in PdfReader(pdf).pages)

scripts=R/'scripts';scripts.mkdir(exist_ok=True)
original_baseline=R.parent/'denoise_and_field_20261023/scripts/fk_baseline.py'
if original_baseline.exists():shutil.copy2(original_baseline,scripts/'fk_baseline.py')
raw_source=R.parent/'denoise_and_field_20261023/field_data/L23535.SGY'
if not raw_source.exists():raw_source=R/'raw/L23535.SGY'
manifest={str(p.relative_to(R)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in R.rglob('*') if p.is_file() and not any(x in p.parts for x in ['tmp','__pycache__']) and p.name not in ['MANIFEST.json','DELIVERY_CHECK.json','build_delivery.log']}
(R/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
dest=R.parent/'field_pilot_20261023.zip'
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(R.rglob('*')):
        if p.is_file() and not any(x in p.parts for x in ['tmp','__pycache__']) and p.name not in ['DELIVERY_CHECK.json','build_delivery.log']:z.write(p,'field_pilot_20261023/'+str(p.relative_to(R)).replace('\\','/'))
    if raw_source!=R/'raw/L23535.SGY':z.write(raw_source,'field_pilot_20261023/raw/L23535.SGY')
with zipfile.ZipFile(dest) as z:
    assert z.testzip() is None
    for k,h in manifest.items():assert hashlib.sha256(z.read('field_pilot_20261023/'+k)).hexdigest()==h
    assert hashlib.sha256(z.read('field_pilot_20261023/raw/L23535.SGY')).hexdigest()==res['source_sha256']
checks=dict(pdf=str(pdf),pages=len(doc),zip=str(dest),bytes=dest.stat().st_size,all_members_verified=True,raw_unchanged=True)
(R/'DELIVERY_CHECK.json').write_text(json.dumps(checks,indent=2),encoding='utf-8');print(json.dumps(checks,ensure_ascii=True))
