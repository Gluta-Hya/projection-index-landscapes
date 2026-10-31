import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import re
s = open(_os.path.join(MANUSCRIPT, "naskah.md"), encoding="utf-8").read()
badan, sisa = s.split("# Daftar pustaka", 1)
pustaka = sisa.split("# Lampiran A", 1)[0]
lampiran = sisa.split("# Lampiran A", 1)[1]
teks = badan + lampiran
# entri pustaka: nama belakang penulis pertama + tahun (+ penulis kedua untuk dua penulis)
entri = []
for baris in pustaka.splitlines():
    m = re.match(r"^([A-ZÀ-ÿ][^,]+), .*?\((\d{4})\)", baris)
    if m:
        penulis = re.findall(r"(?:^|&\s|\.,\s)([A-ZÀ-ÿ][A-Za-zÀ-ÿ'\-]+), [A-Z]\.", baris.split("(" + m.group(2) + ")")[0])
        entri.append((m.group(1), m.group(2), penulis))
print("entri pustaka:", len(entri))
# sitasi di teks
pola = re.compile(r"([A-ZÀ-ÿ][A-Za-zÀ-ÿ'\-]+)(?: (?:&|dan) ([A-ZÀ-ÿ][A-Za-zÀ-ÿ'\-]+)| dkk\.)?,? \(?(\d{4})\)?")
dikutip = set()
for m in pola.finditer(teks):
    a, b, th = m.group(1), m.group(2), m.group(3)
    if not (1900 <= int(th) <= 2026):
        continue
    cocok = [e for e in entri if e[0] == a and e[1] == th]
    if cocok:
        dikutip.add((a, th))
    elif a not in ("Fase", "Tabel", "Gambar", "Bagian", "Lampiran", "Persamaan", "Draf", "Oktober", "ORCID", "Warwick"):
        print("SITASI TANPA ENTRI?", m.group(0))
for e in entri:
    if (e[0], e[1]) not in dikutip:
        print("ENTRI TIDAK DIKUTIP:", e[0], e[1])
print("dikutip:", len(dikutip), "dari", len(entri))
# periksa kecocokan penulis kedua untuk sitasi dua penulis
for m in re.finditer(r"([A-ZÀ-ÿ][A-Za-zÀ-ÿ'\-]+) & ([A-ZÀ-ÿ][A-Za-zÀ-ÿ'\-]+), (\d{4})", teks):
    e = [x for x in entri if x[0] == m.group(1) and x[1] == m.group(3)]
    ok = e and m.group(2) in e[0][2]
    print(("OK  " if ok else "BEDA") , m.group(0), "|", e[0][2] if e else None)
for m in re.finditer(r"([A-ZÀ-ÿ][A-Za-zÀ-ÿ'\-]+) dkk\.,? \(?(\d{4})", teks):
    e = [x for x in entri if x[0] == m.group(1) and x[1] == m.group(2)]
    print(("OK  " if e and len(e[0][2]) >= 3 else "BEDA"), m.group(0), "| jumlah penulis:", len(e[0][2]) if e else None)
