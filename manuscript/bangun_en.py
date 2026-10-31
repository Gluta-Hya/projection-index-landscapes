"""Bangun naskah .docx: pandoc (persamaan Word asli) lalu perapian tata letak dengan python-docx."""
import copy
import subprocess
import sys
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Cm, RGBColor, Emu
from docx.text.paragraph import Paragraph

SUMBER, MENTAH, KELUARAN = "naskah_en.md", "mentah_en.docx", sys.argv[1] if len(sys.argv) > 1 else "naskah_en.docx"
subprocess.run(["pandoc", SUMBER, "-o", MENTAH, "--resource-path=."], check=True)
doc = Document(MENTAH)

HURUF = "Times New Roman"
HITAM, ABU = RGBColor(0, 0, 0), RGBColor(0x59, 0x59, 0x59)


def bersihkan_rpr(rpr, nama=HURUF):
    """Huruf eksplisit di keempat slot, tanpa huruf dan warna tema (tema mengalahkan nilai eksplisit)."""
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts"); rpr.insert(0, rfonts)
    for slot in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rfonts.set(qn(slot), nama)
    for tema in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
        if rfonts.get(qn(tema)) is not None:
            del rfonts.attrib[qn(tema)]
    for warna in rpr.findall(qn("w:color")):
        for atr in ("w:themeColor", "w:themeShade", "w:themeTint"):
            if warna.get(qn(atr)) is not None:
                del warna.attrib[qn(atr)]


def set_huruf(font, nama=HURUF):
    font.name = nama
    bersihkan_rpr(font._element.get_or_add_rPr(), nama)


def cari_gaya(nama):
    for s in doc.styles:
        if s.name and s.name.lower() == nama.lower():
            return s
    return None


def gaya(nama, ukuran=11, tebal=None, miring=None, rata=None, sebelum=None, sesudah=None,
         spasi=None, kiri=None, kanan=None, pertama=None, warna=HITAM, tetap_bersama=None):
    s = cari_gaya(nama)
    if s is None:
        print("PERINGATAN: gaya tidak ditemukan:", nama)
        return None
    s.font.name = HURUF
    s.font.size = Pt(ukuran)
    s.font.color.rgb = warna
    bersihkan_rpr(s.element.get_or_add_rPr())
    if tebal is not None: s.font.bold = tebal
    if miring is not None: s.font.italic = miring
    pf = s.paragraph_format
    if rata is not None: pf.alignment = rata
    if sebelum is not None: pf.space_before = Pt(sebelum)
    if sesudah is not None: pf.space_after = Pt(sesudah)
    if spasi is not None:
        pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE; pf.line_spacing = spasi
    if kiri is not None: pf.left_indent = Cm(kiri)
    if kanan is not None: pf.right_indent = Cm(kanan)
    if pertama is not None: pf.first_line_indent = Cm(pertama)
    if tetap_bersama is not None: pf.keep_with_next = tetap_bersama
    return s


J, C, L = WD_ALIGN_PARAGRAPH.JUSTIFY, WD_ALIGN_PARAGRAPH.CENTER, WD_ALIGN_PARAGRAPH.LEFT
for nama in ("Normal", "Body Text", "First Paragraph", "Compact"):
    gaya(nama, 12, rata=J, sebelum=0, sesudah=6, spasi=1.5)
gaya("Compact", 12, rata=L, sebelum=0, sesudah=3, spasi=1.5)
gaya("Heading 1", 12, tebal=True, miring=False, rata=L, sebelum=14, sesudah=6, spasi=1.0, tetap_bersama=True)
gaya("Heading 2", 12, tebal=True, miring=False, rata=L, sebelum=10, sesudah=4, spasi=1.0, tetap_bersama=True)
gaya("PaperTitle", 16, tebal=True, rata=C, sebelum=0, sesudah=6, spasi=1.1)
gaya("AuthorName", 12, rata=C, sebelum=0, sesudah=2)
gaya("AuthorAffiliation", 10, rata=C, sebelum=0, sesudah=0, spasi=1.1)
gaya("AbstractHeading", 11, tebal=True, rata=C, sebelum=6, sesudah=4, tetap_bersama=True)
gaya("AbstractBody", 11, rata=J, sebelum=0, sesudah=4, spasi=1.3, kiri=1.0, kanan=1.0)
gaya("KeywordLine", 10.5, rata=L, sebelum=0, sesudah=6, spasi=1.2, kiri=1.0, kanan=1.0)
gaya("TheoremHeading", 12, tebal=True, rata=L, sebelum=8, sesudah=2, tetap_bersama=True)
gaya("TheoremBody", 12, miring=True, rata=J, sebelum=0, sesudah=4, spasi=1.5, tetap_bersama=True)
gaya("ReferenceEntry", 11, rata=L, sebelum=0, sesudah=4, spasi=1.15, kiri=1.0, pertama=-1.0)
gaya("Image Caption", 10.5, miring=False, rata=J, sebelum=4, sesudah=12, spasi=1.05)
gaya("Table Caption", 10.5, miring=False, rata=J, sebelum=10, sesudah=4, spasi=1.05, tetap_bersama=True)
gaya("Captioned Figure", 11, rata=C, sebelum=6, sesudah=0, tetap_bersama=True)
vc = cari_gaya("Verbatim Char")
if vc is not None:
    vc.font.name, vc.font.size = "Courier New", Pt(10)
    bersihkan_rpr(vc.element.get_or_add_rPr(), "Courier New")

# bahasa dokumen: Indonesia (pemeriksa ejaan Word)
for s in doc.styles:
    try:
        rpr = s.element.get_or_add_rPr()
    except Exception:
        continue
    lang = rpr.find(qn("w:lang"))
    if lang is None:
        lang = OxmlElement("w:lang"); rpr.append(lang)
    lang.set(qn("w:val"), "en-US")

# ---------- halaman ----------
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
for sisi in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
    setattr(sec, sisi, Cm(2.5))
sec.header_distance, sec.footer_distance = Cm(1.2), Cm(1.2)
LEBAR_TEKS = Cm(16.0)

# nomor baris (setiap 5 baris, berlanjut) untuk memudahkan komentar pembaca
sectPr = sec._sectPr
ln = OxmlElement("w:lnNumType")
ln.set(qn("w:countBy"), "5"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "283")
pgnum = sectPr.find(qn("w:pgNumType"))
(sectPr.insert(sectPr.index(pgnum), ln) if pgnum is not None else sectPr.append(ln))


def tambah_field(paragraf, instr):
    for jenis, teks in (("begin", None), (None, instr), ("end", None)):
        r = OxmlElement("w:r")
        if jenis:
            fc = OxmlElement("w:fldChar"); fc.set(qn("w:fldCharType"), jenis); r.append(fc)
        else:
            it = OxmlElement("w:instrText"); it.set(qn("xml:space"), "preserve"); it.text = f" {teks} "; r.append(it)
        paragraf._p.append(r)


kaki = sec.footer.paragraphs[0]
kaki.alignment = C
tambah_field(kaki, "PAGE")
for r in kaki._p.iter(qn("w:r")):
    rpr = OxmlElement("w:rPr"); sz = OxmlElement("w:sz"); sz.set(qn("w:val"), "18"); rpr.append(sz); r.insert(0, rpr)
kepala = sec.header.paragraphs[0]
kepala.alignment = WD_ALIGN_PARAGRAPH.RIGHT
rk = kepala.add_run("Projection-index landscapes of delay-embedded attractors")
rk.font.size, rk.font.color.rgb, rk.font.italic = Pt(8), ABU, True
set_huruf(rk.font)
for par in (kaki, kepala):  # kepala dan kaki halaman tidak ikut bernomor baris
    par._p.get_or_add_pPr().append(OxmlElement("w:suppressLineNumbers"))
pgMar = sectPr.find(qn("w:pgMar"))
if pgMar is not None and pgMar.get(qn("w:gutter")) is None:
    pgMar.set(qn("w:gutter"), "0")


# ---------- persamaan bernomor: tabel tiga kolom tanpa garis ----------
def tanpa_garis(tbl):
    tblPr = tbl.find(qn("w:tblPr"))
    if tblPr is None:
        tblPr = OxmlElement("w:tblPr"); tbl.insert(0, tblPr)
    st = tblPr.find(qn("w:tblStyle"))
    if st is not None:
        tblPr.remove(st)
    b = tblPr.find(qn("w:tblBorders"))
    if b is not None:
        tblPr.remove(b)
    b = OxmlElement("w:tblBorders")
    for sisi in ("top", "left", "bottom", "right", "insideH", "insideV"):
        e = OxmlElement(f"w:{sisi}"); e.set(qn("w:val"), "nil"); b.append(e)
    tblPr.append(b)
    lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tblPr.append(lay)
    mar = OxmlElement("w:tblCellMar")
    for sisi in ("left", "right"):
        e = OxmlElement(f"w:{sisi}"); e.set(qn("w:w"), "0"); e.set(qn("w:type"), "dxa"); mar.append(e)
    tblPr.append(mar)


def jadikan_persamaan(p, label):
    lebar = [Cm(1.4), Cm(13.2), Cm(1.4)]
    t = doc.add_table(rows=1, cols=3)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    tanpa_garis(t._tbl)
    for i, w in enumerate(lebar):
        t.columns[i].width = w
        sel = t.rows[0].cells[i]
        sel.width = w
        sel.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        pf = sel.paragraphs[0].paragraph_format
        pf.space_before, pf.space_after = Pt(3), Pt(3)
    tengah = t.rows[0].cells[1].paragraphs[0]
    tengah.alignment = C
    for m in list(p._p):
        if m.tag in (qn("m:oMathPara"), qn("m:oMath")):
            tengah._p.append(m)
    kanan = t.rows[0].cells[2].paragraphs[0]
    kanan.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    r = kanan.add_run(f"({label})"); r.font.size = Pt(12); set_huruf(r.font)
    p._p.addprevious(t._tbl)
    p._p.getparent().remove(p._p)


n_utama = n_lamp = 0
for p in list(doc.paragraphs):
    if p.style.name == "Equation":
        n_utama += 1; jadikan_persamaan(p, str(n_utama))
    elif p.style.name == "EquationA":
        n_lamp += 1; jadikan_persamaan(p, f"A.{n_lamp}")
print(f"persamaan bernomor: {n_utama} utama, {n_lamp} lampiran")


# ---------- tabel data: garis gaya jurnal ----------
def garis_jurnal(tabel):
    tblPr = tabel._tbl.tblPr
    st = tblPr.find(qn("w:tblStyle"))
    if st is not None:
        tblPr.remove(st)
    for lama in tblPr.findall(qn("w:tblBorders")):
        tblPr.remove(lama)
    b = OxmlElement("w:tblBorders")
    for sisi, ukuran in (("top", 10), ("bottom", 10)):
        e = OxmlElement(f"w:{sisi}"); e.set(qn("w:val"), "single"); e.set(qn("w:sz"), str(ukuran)); e.set(qn("w:color"), "000000"); b.append(e)
    for sisi in ("left", "right", "insideV", "insideH"):
        e = OxmlElement(f"w:{sisi}"); e.set(qn("w:val"), "nil"); b.append(e)
    tblPr.append(b)
    tabel.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, baris in enumerate(tabel.rows):
        for sel in baris.cells:
            tcPr = sel._tc.get_or_add_tcPr()
            if i == 0:
                tb = OxmlElement("w:tcBorders"); e = OxmlElement("w:bottom")
                e.set(qn("w:val"), "single"); e.set(qn("w:sz"), "6"); e.set(qn("w:color"), "000000"); tb.append(e); tcPr.append(tb)
            for par in sel.paragraphs:
                par.paragraph_format.space_before, par.paragraph_format.space_after = Pt(2), Pt(2)
                par.paragraph_format.line_spacing = 1.0
                for run in par.runs:
                    run.font.size = Pt(9.5)
                    set_huruf(run.font)
                    if i == 0:
                        run.font.bold = True
    # baris judul diulang bila tabel terpotong halaman
    trPr = tabel.rows[0]._tr.get_or_add_trPr()
    h = OxmlElement("w:tblHeader"); h.set(qn("w:val"), "true"); trPr.append(h)


data_tabel = [t for t in doc.tables if len(t.columns) != 3 or len(t.rows) != 1]
for t in data_tabel:
    garis_jurnal(t)
print(f"tabel data: {len(data_tabel)}")

# matematika di sel tabel ikut berukuran 9 pt
for t in data_tabel:
    for baris in t.rows:
        for sel in baris.cells:
            for r in sel._tc.iter(qn("m:r")):
                rpr = r.find(qn("w:rPr"))
                if rpr is None:
                    rpr = OxmlElement("w:rPr"); r.insert(0, rpr)
                sz = OxmlElement("w:sz"); sz.set(qn("w:val"), "19"); rpr.append(sz)


# ---------- tabel data: lebar kolom, tidak terpotong halaman, jarak sesudah ----------
LEBAR_KOLOM = [  # cm, urutan kemunculan: Tabel 1, Tabel 2, Tabel 3, Tabel B1
    [1.4, 2.6, 2.4, 2.4, 2.4, 2.4, 2.4],
    [3.4, 2.8, 2.2, 2.4, 2.4, 2.8],
    [4.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
    [1.2, 3.4, 4.0, 4.4, 1.6, 1.4],
]
assert [len(t.columns) for t in data_tabel] == [len(x) for x in LEBAR_KOLOM], [len(t.columns) for t in data_tabel]


def atur_lebar(tabel, lebar):
    tbl = tabel._tbl
    tblPr = tbl.tblPr
    tblW = tblPr.find(qn("w:tblW"))
    if tblW is None:
        tblW = OxmlElement("w:tblW"); tblPr.append(tblW)
    tblW.set(qn("w:w"), str(round(sum(lebar) * 567))); tblW.set(qn("w:type"), "dxa")
    for lama in tblPr.findall(qn("w:tblLayout")):
        tblPr.remove(lama)
    lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tblPr.append(lay)
    for gc, w in zip(tbl.tblGrid.findall(qn("w:gridCol")), lebar):
        gc.set(qn("w:w"), str(round(w * 567)))
    for baris in tabel.rows:
        for sel, w in zip(baris.cells, lebar):
            sel.width = Cm(w)


def tidak_terpotong(tabel):
    n = len(tabel.rows)
    for i, baris in enumerate(tabel.rows):
        trPr = baris._tr.get_or_add_trPr()
        if trPr.find(qn("w:cantSplit")) is None:
            trPr.append(OxmlElement("w:cantSplit"))
        if i < n - 1:
            for sel in baris.cells:
                for par in sel.paragraphs:
                    par.paragraph_format.keep_with_next = True


for t, lebar in zip(data_tabel, LEBAR_KOLOM):
    atur_lebar(t, lebar)
    tidak_terpotong(t)
    berikut = t._tbl.getnext()
    if berikut is not None and berikut.tag == qn("w:p"):
        Paragraph(berikut, t._parent).paragraph_format.space_before = Pt(8)


# ---------- urutan elemen sesuai skema OOXML (Word ketat soal ini) ----------
W_ = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
M_ = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"
URUTAN = {
    W_ + "pPr": "pStyle keepNext keepLines pageBreakBefore framePr widowControl numPr suppressLineNumbers pBdr shd tabs "
                "suppressAutoHyphens kinsoku wordWrap overflowPunct topLinePunct autoSpaceDE autoSpaceDN bidi "
                "adjustRightInd snapToGrid spacing ind contextualSpacing mirrorIndents suppressOverlap jc textDirection "
                "textAlignment textboxTightWrap outlineLvl divId cnfStyle rPr sectPr pPrChange",
    W_ + "rPr": "ins del moveFrom moveTo rStyle rFonts b bCs i iCs caps smallCaps strike dstrike outline shadow emboss imprint noProof "
                "snapToGrid vanish webHidden color spacing w kern position sz szCs highlight u effect bdr shd fitText "
                "vertAlign rtl cs em lang eastAsianLayout specVanish oMath",
    W_ + "tblPr": "tblStyle tblpPr tblOverlap bidiVisual tblStyleRowBandSize tblStyleColBandSize tblW jc tblCellSpacing "
                  "tblInd tblBorders shd tblLayout tblCellMar tblLook tblCaption tblDescription tblPrChange",
    W_ + "tcPr": "cnfStyle tcW gridSpan hMerge vMerge tcBorders shd noWrap tcMar textDirection tcFitText vAlign hideMark "
                 "headers cellIns cellDel cellMerge tcPrChange",
    W_ + "tblBorders": "top start left bottom end right insideH insideV",
    W_ + "tcBorders": "top start left bottom end right insideH insideV tl2br tr2bl",
    M_ + "rPr": "lit nor scr sty brk aln",
    M_ + "dPr": "begChr sepChr endChr grow shp ctrlPr",
}
URUTAN = {induk: [(M_ if induk.startswith(M_) else W_) + n for n in s.split()] for induk, s in URUTAN.items()}
URUTAN[M_ + "r"] = [M_ + "rPr", W_ + "rPr"]


def urutkan(e, urutan):
    idx = {t: i for i, t in enumerate(urutan)}
    anak = list(e)
    kunci = sorted(range(len(anak)), key=lambda k: (idx.get(anak[k].tag, len(urutan)), k))
    if kunci != list(range(len(anak))):
        for c in anak:
            e.remove(c)
        for k in kunci:
            e.append(anak[k])


def buang_ganda(e, tag_tag):
    terakhir = {}
    for c in list(e):
        if c.tag in tag_tag:
            if c.tag in terakhir:
                e.remove(terakhir[c.tag])
            terakhir[c.tag] = c


def rapikan_skema(akar):
    for e in list(akar.iter()):  # daftar dibekukan dulu: anak dipindah selama perulangan
        if e.tag == W_ + "tblPr":
            buang_ganda(e, set(URUTAN[e.tag]))
        elif e.tag == W_ + "trPr":
            buang_ganda(e, {W_ + "tblHeader", W_ + "cantSplit"})
        elif e.tag == M_ + "rPr" and e.find(M_ + "nor") is not None:
            for x in e.findall(M_ + "sty") + e.findall(M_ + "scr"):
                e.remove(x)  # nor dan (scr, sty) saling meniadakan di skema
        if e.tag in URUTAN:
            urutkan(e, URUTAN[e.tag])


rapikan_skema(doc.element.body)
rapikan_skema(sec.header._element)
rapikan_skema(sec.footer._element)

# bagian lain dari dokumen acuan bawaan pandoc: gaya, pengaturan, penomoran
URUTAN[W_ + "style"] = [W_ + n for n in (
    "name aliases basedOn next link autoRedefine hidden uiPriority semiHidden unhideWhenUsed qFormat locked "
    "personal personalCompose personalReply rsid pPr rPr tblPr trPr tcPr tblStylePr").split()]
gaya_akar = doc.styles.element
for r in gaya_akar.iter(W_ + "rPr"):  # dokumen acuan pandoc memuat karakter '>' nyasar di satu rPr
    r.text = None
    for c in r:
        c.tail = None
rapikan_skema(gaya_akar)
URUTAN_SETTINGS = (
    "writeProtection view zoom removePersonalInformation removeDateAndTime doNotDisplayPageBoundaries "
    "displayBackgroundShape printPostScriptOverText printFractionalCharacterWidth printFormsData embedTrueTypeFonts "
    "embedSystemFonts saveSubsetFonts saveFormsData mirrorMargins alignBordersAndEdges bordersDoNotSurroundHeader "
    "bordersDoNotSurroundFooter gutterAtTop hideSpellingErrors hideGrammaticalErrors activeWritingStyle proofState "
    "formsDesign attachedTemplate linkStyles stylePaneFormatFilter stylePaneSortMethod documentType mailMerge "
    "revisionView trackRevisions doNotTrackMoves doNotTrackFormatting documentProtection autoFormatOverride "
    "styleLockTheme styleLockQFSet defaultTabStop autoHyphenation consecutiveHyphenLimit hyphenationZone "
    "doNotHyphenateCaps showEnvelope summaryLength clickAndTypeStyle defaultTableStyle evenAndOddHeaders "
    "bookFoldRevPrinting bookFoldPrinting bookFoldPrintingSheets drawingGridHorizontalSpacing "
    "drawingGridVerticalSpacing displayHorizontalDrawingGridEvery displayVerticalDrawingGridEvery "
    "doNotUseMarginsForDrawingGridOrigin drawingGridHorizontalOrigin drawingGridVerticalOrigin doNotShadeFormData "
    "noPunctuationKerning characterSpacingControl printTwoOnOne strictFirstAndLastChars noLineBreaksAfter "
    "noLineBreaksBefore savePreviewPicture doNotValidateAgainstSchema saveInvalidXml ignoreMixedContent "
    "alwaysShowPlaceholderText doNotDemarcateInvalidXml saveXmlDataOnly useXSLTWhenSaving saveThroughXslt "
    "showXMLTags alwaysMergeEmptyNamespace updateFields hdrShapeDefaults footnotePr endnotePr compat docVars rsids "
    "mathPr attachedSchema themeFontLang clrSchemeMapping doNotIncludeSubdocsInStats doNotAutoCompressPictures "
    "forceUpgrade captions readModeInkLockDown smartTagType schemaLibrary shapeDefaults doNotEmbedSmartTags "
    "decimalSymbol listSeparator").split()
urutkan(doc.settings.element, [(M_ if n == "mathPr" else W_) + n for n in URUTAN_SETTINGS])
try:
    for e in doc.part.numbering_part.element.iter(W_ + "nsid"):  # nsid wajib 8 digit heksadesimal
        e.set(qn("w:val"), e.get(qn("w:val")).upper().zfill(8))
    for e in doc.part.numbering_part.element.iter(W_ + "lvlText"):  # butir bertingkat: bulatan, bukan tanda pisah
        if e.get(qn("w:val")) == "\u2013":
            e.set(qn("w:val"), "◦")
except (KeyError, NotImplementedError):
    pass

# ---------- metadata ----------
cp = doc.core_properties
cp.title = "Projection-index landscapes of delay-embedded chaotic attractors: period scaling and its mechanism"
cp.author = "Athaillah Farrel Inandra"
cp.subject = "Manuscript submitted to Dynamical Systems: An International Journal"
cp.keywords = "delay embedding; projection pursuit; bimodality coefficient; Lorenz attractor; surrogate data"
cp.language = "en-US"

doc.save(KELUARAN)
print("tersimpan:", KELUARAN)
