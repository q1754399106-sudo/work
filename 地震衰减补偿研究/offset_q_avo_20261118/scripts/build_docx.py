"""Build the Word manuscript for Oil Geophysical Prospecting from PAPER_DRAFT_ZH.md.

usage: build_docx.py   (needs pandoc on PATH)
-> submission/manuscript_ogp.md (intermediate) and submission/叠前Q补偿AVO误差预算_投稿稿.docx
Keeps only the manuscript part (drops everything after "# 以下为内部工作记录" and the editorial note under the title),
inserts the 8 figures (figures/*.png) before their captions, turns \\tag{n} into a right-hand "(n)", and converts
plain-text subscripts / superscripts (a_x, σ_lnQ, Q_模型, a_x^m, 10^−3.5) outside $...$ into Word sub/superscripts.
Display equations become editable Word (OMML) equations. Styles come from submission/reference_ogp.docx
(pandoc default template with body 宋体 + Times New Roman 10.5 pt, headings 黑体 in black).
"""
import re, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "submission"; OUT.mkdir(exist_ok=True)
NAMES = {1: "fig1_workflow", 2: "fig2_gamma_check", 3: "fig3_avo_vs_sigma", 4: "fig4_uq_cost_vs_lambda",
         5: "fig5_window_compression", 6: "fig6_chat_accuracy", 7: "fig7_end_to_end", 8: "fig8_example_gather"}
SUB = re.compile(r"(?<![A-Za-z0-9~])([A-Za-zα-ωΓΔΦΨ])_([A-Za-z0-9]+|模型|真)(\^([A-Za-z0-9*]+))?")
SUP = re.compile(r"(\d+)\^([−\-]?[0-9.]+)")


def text_scripts(t):
    t = SUB.sub(lambda m: f"{m.group(1)}~{m.group(2)}~" + (f"^{m.group(4)}^" if m.group(4) else ""), t)
    return SUP.sub(lambda m: f"{m.group(1)}^{m.group(2)}^", t)


def convert(s):
    s = s.split("# 以下为内部工作记录")[0]
    s = re.sub(r"\n（《石油地球物理勘探》投稿稿.*?）\n", "\n", s, count=1)
    s = re.sub(r"\\tag\{(\d+)\}", r"\\qquad (\1)", s)
    parts = re.split(r"(\$\$.*?\$\$|\$[^$\n]+\$)", s, flags=re.S)        # leave maths untouched
    s = "".join(p if p.startswith("$") else text_scripts(p) for p in parts)
    s = re.sub(r"^图 (\d)　.*$", lambda m: f"![](../figures/{NAMES[int(m.group(1))]}.png){{width={'8.4cm' if m.group(1) == '2' else '17.4cm'}}}\n\n" + m.group(0), s, flags=re.M)
    return s.rstrip().rstrip("-").rstrip() + "\n"


if __name__ == "__main__":
    md = OUT / "manuscript_ogp.md"; md.write_text(convert((ROOT / "PAPER_DRAFT_ZH.md").read_text(encoding="utf-8")), encoding="utf-8")
    docx = OUT / "叠前Q补偿AVO误差预算_投稿稿.docx"
    subprocess.run(["pandoc", md.name, "-f", "markdown+tex_math_dollars+subscript+superscript", "-o", docx.name, "--resource-path=.", "--reference-doc=reference_ogp.docx"],
                   cwd=OUT, check=True)
    print("written", docx)
