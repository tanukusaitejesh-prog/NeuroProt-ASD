"""Render a manuscript Markdown file to a submission-style PDF, with its figures embedded.

No pandoc or weasyprint on this machine, so this builds the document directly with reportlab platypus. It
handles the Markdown subset the manuscripts actually use: headings, paragraphs with bold/italic/code, pipe
tables, bullet lists, horizontal rules, and figure legends (the matching PNG is inserted above each legend).

Usage:  python scripts/62_make_pdf.py MANUSCRIPT_RECURRENCE.md [out.pdf]
"""
import html
import re
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (HRFlowable, Image, KeepTogether, PageBreak, Paragraph,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else config.ROOT / "MANUSCRIPT_RECURRENCE.md"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else SRC.with_suffix(".pdf")
FIGSETS = {
    "MANUSCRIPT_RECURRENCE": {1: "recur_fig1_rate.png", 2: "recur_fig2_homopolymer.png",
                              3: "recur_fig3_residual.png"},
    "MANUSCRIPT_MuST-VEP": {1: "fig1_workflow.png", 2: "fig2_dissociation.png"},
}
FIGS = FIGSETS.get(SRC.stem, {})

INK = colors.HexColor("#111111")
MUTED = colors.HexColor("#555555")
RULE = colors.HexColor("#cccccc")
BAND = colors.HexColor("#f0f0ee")

# The built-in Type 1 fonts are WinAnsi-only: Greek renders blank and so do several symbols. matplotlib ships
# DejaVu, which is full Unicode, so register that instead and write rho / Delta / arrows literally.
import matplotlib
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

_TTF = Path(matplotlib.__file__).parent / "mpl-data" / "fonts" / "ttf"
for _name, _file in [("DJSerif", "DejaVuSerif.ttf"), ("DJSerif-Bold", "DejaVuSerif-Bold.ttf"),
                     ("DJSerif-Italic", "DejaVuSerif-Italic.ttf"),
                     ("DJSerif-BoldItalic", "DejaVuSerif-BoldItalic.ttf"),
                     ("DJSans", "DejaVuSans.ttf"), ("DJSans-Bold", "DejaVuSans-Bold.ttf"),
                     ("DJSans-Oblique", "DejaVuSans-Oblique.ttf"),
                     ("DJSans-BoldOblique", "DejaVuSans-BoldOblique.ttf"),
                     ("DJMono", "DejaVuSansMono.ttf")]:
    pdfmetrics.registerFont(TTFont(_name, str(_TTF / _file)))
pdfmetrics.registerFontFamily("DJSerif", normal="DJSerif", bold="DJSerif-Bold",
                              italic="DJSerif-Italic", boldItalic="DJSerif-BoldItalic")
pdfmetrics.registerFontFamily("DJSans", normal="DJSans", bold="DJSans-Bold",
                              italic="DJSans-Oblique", boldItalic="DJSans-BoldOblique")

ss = getSampleStyleSheet()
BODY = ParagraphStyle("body", parent=ss["Normal"], fontName="DJSerif", fontSize=9.0, leading=13.2,
                      alignment=TA_JUSTIFY, textColor=INK, spaceAfter=5)
H1 = ParagraphStyle("h1", parent=BODY, fontName="DJSans-Bold", fontSize=14, leading=18.5, alignment=0,
                    spaceBefore=4, spaceAfter=9)
H2 = ParagraphStyle("h2", parent=BODY, fontName="DJSans-Bold", fontSize=11, leading=14.5, alignment=0,
                    spaceBefore=13, spaceAfter=5)
H3 = ParagraphStyle("h3", parent=BODY, fontName="DJSans-BoldOblique", fontSize=9.4, leading=13, alignment=0,
                    spaceBefore=10, spaceAfter=4)
META = ParagraphStyle("meta", parent=BODY, fontName="DJSerif-Italic", fontSize=8.0, leading=11.2,
                      textColor=MUTED, alignment=0, spaceAfter=3)
CELL = ParagraphStyle("cell", parent=BODY, fontSize=7.2, leading=9.6, alignment=0, spaceAfter=0)
CELLH = ParagraphStyle("cellh", parent=CELL, fontName="DJSans-Bold")
CAP = ParagraphStyle("cap", parent=BODY, fontSize=8.0, leading=11.2, textColor=MUTED, alignment=0)
BUL = ParagraphStyle("bul", parent=BODY, leftIndent=11, bulletIndent=2, spaceAfter=3)

SUBS = [("‘", "'"), ("’", "'"), ("“", '"'), ("”", '"')]
SUP = {"⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
       "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}
SUB = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5"}


def inline(t):
    """Markdown inline -> reportlab mini-HTML."""
    t = html.escape(t, quote=False)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<!\*)\*(?!\s)(.+?)(?<!\s)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"`(.+?)`", r'<font face="Courier" size="8.6">\1</font>', t)
    t = re.sub(r"\[(.+?)\]\((.+?)\)", r"\1", t)
    # runs of superscript / subscript digits
    t = re.sub("[" + "".join(SUP) + "]+", lambda m: "<super>" + "".join(SUP[c] for c in m.group()) + "</super>", t)
    t = re.sub("[" + "".join(SUB) + "]+", lambda m: "<sub>" + "".join(SUB[c] for c in m.group()) + "</sub>", t)
    for a, b in SUBS:
        t = t.replace(a, b)
    return t


def table_flow(rows, width):
    head, body = rows[0], rows[1:]
    ncol = len(head)
    data = [[Paragraph(inline(c), CELLH) for c in head]]
    data += [[Paragraph(inline(c), CELL) for c in r] for r in body]
    first = max(28 * mm, width * 0.30) if ncol <= 4 else width * 0.26
    rest = (width - first) / max(ncol - 1, 1)
    t = Table(data, colWidths=[first] + [rest] * (ncol - 1), repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), BAND),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, INK),
        ("LINEBELOW", (0, -1), (-1, -1), 0.7, INK),
        ("INNERGRID", (0, 1), (-1, -1), 0.25, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.6), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.6),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4)]))
    return t


def figure_flow(n, width):
    p = config.FIGURES / FIGS[n]
    if not p.exists():
        return None
    from PIL import Image as PILImage
    with PILImage.open(p) as im:
        w, h = im.size
    iw = min(width, 170 * mm)
    return Image(str(p), width=iw, height=iw * h / w)


def build(md_path, pdf_path):
    text = md_path.read_text(encoding="utf8")
    doc = SimpleDocTemplate(str(pdf_path), pagesize=A4, title=md_path.stem,
                            leftMargin=20 * mm, rightMargin=20 * mm, topMargin=18 * mm, bottomMargin=18 * mm)
    W = doc.width
    flow, lines, i = [], text.split("\n"), 0
    para, bullets = [], []

    def flush_para():
        nonlocal para
        if para:
            flow.append(Paragraph(inline(" ".join(para)), BODY))
            para = []

    def flush_bullets():
        nonlocal bullets
        for b in bullets:
            flow.append(Paragraph(inline(b), BUL, bulletText="•"))
        if bullets:
            flow.append(Spacer(1, 3))
        bullets = []

    while i < len(lines):
        ln = lines[i].rstrip()
        s = ln.strip()
        if not s:
            flush_para(); flush_bullets(); i += 1; continue
        if s.startswith("|") and i + 1 < len(lines) and set(lines[i + 1].strip()) <= set("|-: "):
            flush_para(); flush_bullets()
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                r = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not set("".join(r)) <= set("-: "):
                    rows.append(r)
                i += 1
            flow.append(Spacer(1, 3)); flow.append(table_flow(rows, W)); flow.append(Spacer(1, 7))
            continue
        if re.fullmatch(r"-{3,}|\*{3,}", s):
            flush_para(); flush_bullets()
            flow.append(Spacer(1, 4)); flow.append(HRFlowable(width="100%", color=RULE, thickness=0.6))
            flow.append(Spacer(1, 6)); i += 1; continue
        m = re.match(r"^(#{1,4})\s+(.*)", s)
        if m:
            flush_para(); flush_bullets()
            lvl, txt = len(m.group(1)), m.group(2)
            flow.append(Paragraph(inline(txt), {1: H1, 2: H2, 3: H3, 4: H3}[lvl]))
            i += 1; continue
        m = re.match(r"^[-*]\s+(.*)", s)
        if m:
            flush_para(); bullets.append(m.group(1)); i += 1; continue
        if s.startswith("*") and s.endswith("*") and not s.startswith("**") and len(s) > 2:
            flush_para(); flush_bullets()
            flow.append(Paragraph(inline(s.strip("*")), META)); i += 1; continue
        fm = re.match(r"^\*\*Figure (\d)\.\*\*", s)
        if fm:
            flush_para(); flush_bullets()
            n = int(fm.group(1))
            legend_lines = [s]
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].strip().startswith("**Figure"):
                legend_lines.append(lines[j].strip()); j += 1
            img = figure_flow(n, W)
            block = ([img, Spacer(1, 3)] if img else []) + [Paragraph(inline(" ".join(legend_lines)), CAP)]
            flow.append(Spacer(1, 6)); flow.append(KeepTogether(block)); flow.append(Spacer(1, 9))
            i = j; continue
        para.append(s); i += 1
    flush_para(); flush_bullets()

    def footer(canv, d):
        canv.saveState(); canv.setFont("Helvetica", 7.5); canv.setFillColor(MUTED)
        canv.drawCentredString(A4[0] / 2, 11 * mm, str(canv.getPageNumber()))
        canv.restoreState()

    doc.build(flow, onFirstPage=footer, onLaterPages=footer)
    return pdf_path


p = build(SRC, OUT)
print(f"wrote {p}  ({p.stat().st_size / 1024:.0f} KB)")
