from pathlib import Path
import csv, io, json
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.lib.styles import ParagraphStyle
from pypdf import PdfReader, PdfWriter, Transformation
import pypdfium2 as pdfium
from PIL import Image, ImageOps, ImageDraw

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output'/'pdf'; OUT.mkdir(parents=True,exist_ok=True)
TMP=ROOT/'tmp'/'pdfs'; TMP.mkdir(parents=True,exist_ok=True)
pdfmetrics.registerFont(TTFont('CN','C:/Windows/Fonts/msyh.ttc',subfontIndex=0))
pdfmetrics.registerFont(TTFont('CNB','C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0))
W,H=landscape(A4); M=38
NAVY=colors.HexColor('#16354b'); TEAL=colors.HexColor('#087f8c')
body=ParagraphStyle('body',fontName='CN',fontSize=11,leading=19,textColor=NAVY,wordWrap='CJK')
small=ParagraphStyle('small',parent=body,fontSize=9,leading=15)
writer=PdfWriter()
def page(title,kicker='相干干扰压制与真实地震资料试验'):
    buf=io.BytesIO(); c=canvas.Canvas(buf,pagesize=(W,H))
    c.setFillColor(TEAL); c.rect(0,H-8,W,8,fill=1,stroke=0)
    c.setFont('CN',9); c.setFillColor(TEAL); c.drawString(M,H-32,kicker)
    c.setFillColor(NAVY); c.setFont('CNB',22); c.drawString(M,H-66,title)
    c.setStrokeColor(colors.HexColor('#d8e3e8')); c.line(M,36,W-M,36)
    c.setFont('CN',8); c.drawString(M,22,'实验日期 2026-10-23  |  同一观测、冻结补偿与标定参数')
    c.drawRightString(W-M,22,f'{len(writer.pages)+1:02d}')
    return buf,c
def para(c,text,y,style=body,width=None):
    p=Paragraph(text,style); _,h=p.wrap(width or W-2*M,H); p.drawOn(c,M,y-h); return y-h-12
def save(buf,c):
    c.save();buf.seek(0);writer.add_page(PdfReader(buf).pages[0])
def table(c,rows,y,widths):
    rows=[[Paragraph(str(v),small) for v in r] for r in rows]
    t=Table(rows,colWidths=widths,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e1eff2')),('VALIGN',(0,0),(-1,-1),'TOP'),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f3f6f8')]),('LEFTPADDING',(0,0),(-1,-1),9),('RIGHTPADDING',(0,0),(-1,-1),9),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.6,TEAL)]))
    _,h=t.wrap(W-2*M,H);t.drawOn(c,M,y-h);return y-h-16

buf,c=page('相干干扰压制与真实地震资料试验','实验结果报告  /  2026-10-23')
y=H-101
y=para(c,'本轮目标：减弱复杂模型中可见的斜向相干干扰，并在公开真实炮集上检查方法是否能够直接工作。',y)
y=table(c,[['主要结果','Moderate','Strong'],['优化方法整体相对误差','0.568153 → 0.563300','0.701367 → 0.684622'],['相对误差降幅','0.85%','2.39%'],['纯净输入控制试验变化','0','0']],y,[310,225,230])
y=para(c,'合成资料：分频长孔径检测识别出约 1846 m/s 的低频波包，斜向干扰减弱，补偿误差小幅下降。较缓、与反射频带重叠的干扰仍有残留。',y)
y=para(c,'真实资料：已下载两套数据并处理四个炮集。主候选均未触发；补充软 f-k 基线减去了部分相干分量，但尚不能确认有效波保真或信噪比改善。',y)
y=para(c,'当前结论：保留为研究候选，暂不替换原补偿流程。实测实验限于去干扰，没有完成可信的真实资料 Q 保幅补偿验证。',y)
save(buf,c)

buf,c=page('01  试验设置与评价方式');y=H-97
for text in [
'复杂合成模型包含 60 个反射层、褶皱、断层、河道、楔状体、低 Q 透镜体、空间变化子波、鬼波、多次波、随机及相干噪声。剖面有 99 道，采样间隔 2 ms，时长 1.7 s；两个观测条件分别为 14 dB 与 8 dB 随机噪声。',
'生成参考采用局部一维射线与褶积近似，并非完整二维弹性波模拟。原始死道保留在全道误差中，另记录有效道误差。',
'第一版宽频长孔径检测未触发。第二版按观测 f90 分频，扫描 1201 个斜率、视速度不低于 400 m/s；至少覆盖 70% 有效孔径，相似性门槛 0.65，每频带最多 3 个波包，减除比例 0.60。预测时按道 RMS 预条件化，减除前恢复原幅度单位。',
'噪声轨迹、生成器速度和纯净参考没有传给预处理。原补偿算法与标定参数保持冻结，仅在预处理后重新估计噪声水平。两个合成场景属于开发资料，不能据此宣称独立泛化。',
'整体相对误差为 ||输出 - 匹配纯净参考||₂ / ||匹配纯净参考||₂，越小越好。低能量区及深部使用原固定评分掩膜；各项误差均在同一幅度标定下比较。']:
    y=para(c,text,y)
save(buf,c)

rows=list(csv.DictReader((ROOT/'analysis/compensation_metrics.csv').open(encoding='utf-8-sig')))
buf,c=page('02  完整量化对比');y=H-96
data=[['场景 / 方法','整体误差：前 → 后','低能量区误差：前 → 后']]
for scene in ['Moderate','Strong']:
    for method in ['FixedSplit','FixedSplitSharedCal','FullCalRefit','PrunedFullCalRefit']:
        b=next(r for r in rows if r['Dataset']==scene and r['Method']==method and r['Stage']=='Before')
        a=next(r for r in rows if r['Dataset']==scene and r['Method']==method and r['Stage']=='After')
        data.append([scene+' / '+method,f"{float(b['RelativeError']):.6f} → {float(a['RelativeError']):.6f}",f"{float(b['LowEnergyError']):.5f} → {float(a['LowEnergyError']):.5f}"])
y=table(c,data,y,[300,225,240])
y=para(c,'FixedSplit：最初分步法；FixedSplitSharedCal：采用共享标定的分步法；FullCalRefit：未剪枝的全标定重拟合；PrunedFullCalRefit：当前优化方法。表中“前 / 后”均指加入去干扰预处理前后。',y,small)
y=para(c,'本轮 FullCalRefit 的整体误差低于当前剪枝版本，说明剪枝带来的稀疏性并未转化为这两个场景下的最低误差。弱反射恢复仍需继续改进。',y,small)
save(buf,c)

def figure_page(title,name,caption):
    buf,c=page(title);para(c,caption,H-83,small);save(buf,c)
    target=writer.pages[-1];src=PdfReader(ROOT/'figures'/f'{name}.pdf').pages[0]
    sw,sh=float(src.mediabox.width),float(src.mediabox.height)
    avail_w=W-2*M;avail_h=H-180;s=min(avail_w/sw,avail_h/sh)
    x=(W-sw*s)/2;y=47+(avail_h-sh*s)/2
    target.merge_transformed_page(src,Transformation().scale(s).translate(x,y),over=True)

figure_page('03  Moderate：去干扰前后','Moderate_preprocessing','同一色标展示输入、候选输出和减去分量。减去的输入范数比例为 7.72%，并非能量比例；低频相干干扰有所减弱。')
figure_page('04  Strong：去干扰前后','Strong_preprocessing','减去的输入范数比例为 12.07%；干扰建模残差范数为原相干干扰的 0.680 倍。相干干扰并未被完整消除。')
figure_page('05  Moderate：补偿结果对比','Moderate_compensation_comparison','比较加入预处理前后的原分步法与优化方法。共用幅度色标，不逐道归一化；整体相对误差分别降至 0.618010 和 0.563300。')
figure_page('06  Strong：补偿结果对比','Strong_compensation_comparison','优化方法整体相对误差 0.701367 → 0.684622，低能量区误差 1.07630 → 1.02139。改善有限，深部及弱反射残差仍然可见。')

buf,c=page('07  真实资料与结果解释');y=H-96
y=table(c,[['数据','实际使用','采样及几何'],['USGS Alaska NPRA 31-81','1981 年陆地炸药源；下载 30 炮 / 3030 道；选 111、120、130 炮各前 96 道','2 ms / 6 s；按发布者脚本重建道距 33.528 m，非原始导航坐标'],['Yilmaz wz.25 / Alberta','97 道；中央第 49 道为原示例插入的零道，保留并标记','2 ms / 4 s；采用头文件声明的 50 m 道距，历史几何不确定性保留']],y,[200,290,275])
y=table(c,[['炮集','主候选减去范数比例','软 f-k 基线减去范数比例'],['wz25','0','0.117725'],['Alaska 111','0','0.378872'],['Alaska 120','0','0.107144'],['Alaska 130','0','0.170972']],y,[200,290,275])
y=para(c,'主候选四炮均未触发，不能视为实测去噪成功。软 f-k 基线的速度过渡为 1200-2500 m/s，频率过渡为观测 f90 的 0.4-0.8 倍，最大减除 0.60。减去分量可能混入有效直达波、折射波或反射波。',y,small)
y=para(c,'没有独立 Q 标定、已确认层厚 / 速度和真实无衰减参考，且炮集未经动校正，因此本轮未开展实测 Q 补偿保幅验证。后续图中 t² 为统一显示增益，不参与数值处理；共用 99.5% 分位色标。',y,small)
save(buf,c)
for title,name in [('08  Alberta：真实炮集对比','wz25'),('09  Alaska：111 炮对比','alaska_111'),('10  Alaska：120 炮对比','alaska_120'),('11  Alaska：130 炮对比','alaska_130')]:
    figure_page(title,name+'_field_comparison','从左至右：原始炮集、主候选输出、软 f-k 基线输出、软 f-k 减去分量。主候选输出等于输入；真实资料缺少无噪声参考，不能根据减去量计算实际 SNR 改善。')

buf,c=page('12  数据来源与复现记录');y=H-98
for label,url in [
('USGS NAMSS 原调查目录','https://walrus.wr.usgs.gov/namss/survey/n-pr-81-ak/'),
('Alaska 发布者处理脚本与几何','https://ahay.org/RSF/book/data/alaska/line31-81-m8r.html'),
('实际 SEG-Y 文件','https://ahay.org/data/alaska/L23535.SGY'),
('Yilmaz 实测去面波复现代码','https://ahay.org/RSF/book/tccs/orthogroll/field.html'),
('实际 wz.25 示例文件','https://ahay.org/data/wz/wz.25.H')]:
    y=para(c,f'{label}<br/><link href="{url}" color="#087f8c">{url}</link>',y,small)
y=para(c,'原始下载、文件头、SHA256、解析脚本、参数、MATLAB 输出、误差表和绘图文件均保存在项目及完整 ZIP 包中。原文件不覆盖；Alaska 使用时注明 USGS 来源。wz25 下载页未发现独立数据许可证，不将软件许可证等同于数据许可证。',y,small)
y=para(c,'合成干扰去除前后的范数守恒、原始下载哈希、两次 MATLAB 退出状态及全部图像已检查。所有 117 个旧 MATLAB 文件保持原哈希；本轮仅增加试验入口。主候选未进入默认生产流程。',y,small)
y=para(c,'复现入口：RESULTS_ZH.md；analysis/compensation_metrics.csv；scripts/multiband_filter.py；scripts/fk_baseline.py；scripts/read_field.py。完整实验包：denoise_and_real_seismic_20261023.zip。',y,small)
save(buf,c)
dest=OUT/'相干干扰压制与真实地震资料实验报告.pdf'
writer.add_metadata({'/Title':'相干干扰压制与真实地震资料实验报告','/Author':'研究实验记录','/Subject':'合成受控对照与公开真实炮集去噪试验'})
with dest.open('wb') as f:writer.write(f)
doc=pdfium.PdfDocument(str(dest)); thumbs=[]
for i in range(len(doc)):
    im=doc[i].render(scale=1.25).to_pil().convert('RGB');im.save(TMP/f'page_{i+1:02d}.png')
    im.thumbnail((560,410));tile=Image.new('RGB',(580,440),'#e5e9ed');tile.paste(im,((580-im.width)//2,10));ImageDraw.Draw(tile).text((15,420),str(i+1),fill='black');thumbs.append(tile)
contact=Image.new('RGB',(580*3,440*((len(thumbs)+2)//3)),'white')
for i,t in enumerate(thumbs):contact.paste(t,((i%3)*580,(i//3)*440))
contact.save(TMP/'contact_sheet.png')
check=PdfReader(dest);assert len(check.pages)==13
assert all(len(p.extract_text() or '')>80 for p in check.pages)
print(json.dumps({'pdf':str(dest),'pages':len(check.pages),'bytes':dest.stat().st_size,'contact_sheet':str(TMP/'contact_sheet.png')},ensure_ascii=False))
