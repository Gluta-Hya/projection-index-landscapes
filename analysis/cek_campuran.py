import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import numpy as np, csv
import fase3_model_mainan as M, fase2_sistem as F
# urutan pemanggilan sama persis dengan skrip Fase 3 -> deret identik
x_sin = np.cos(2*np.pi*np.arange(60_000)*M.DT)
x_jarang = M.dua_sayap(200_000, D=40)
x_lor = M.dua_sayap(60_000, D=2.2)
f = lambda r: 2*(2+r)**2/(3*r*r+24*r+8)
al = np.linspace(0, np.pi, 400000, endpoint=False); q = np.sin(2*al)
def bmix(ac, p):
    b0 = f(ac**2); return np.mean(b0*(1+(2*p-1)*q)**2/(1+2*(2*p-1)*q+q**2))
rows = list(csv.DictReader(open(DATA + 'fase3_hasil.csv', encoding='utf-8')))
sim = {r['model']: float(r['I3_rata']) for r in rows if r['m']=='2' and abs(float(r['W_per_P'])-1.0)<1e-9}
for nama, x, D in (("dua_sayap_jarang", x_jarang, 40), ("dua_sayap_lorenz", x_lor, 2.2)):
    s = np.sign(x); beda = np.mean(s[100:] != s[:-100])   # W = P = 100 sampel
    print(f"{nama}: 1-p nominal = {1/D:.4f}, terealisasi = {beda:.4f};  rumus(nominal) = {bmix(0.6, 1-1/D):.4f}, "
          f"rumus(terealisasi) = {bmix(0.6, 1-beda):.4f}, simulasi = {sim[nama]:.4f}")
    # cek efek pengambilan sampel/arah: hitung ulang lanskap di W=P dengan semua titik dan 4000 arah seragam
    V = F.embed(x, 2, 100); ang = np.linspace(0, np.pi, 4000, endpoint=False); U = np.c_[np.cos(ang), np.sin(ang)]
    I3 = []
    for blok in np.array_split(U, 40):
        Y = V @ blok.T; Y = (Y-Y.mean(0))/Y.std(0); n = len(Y); Y2 = Y*Y; g = (Y2*Y).mean(0); k = (Y2*Y2).mean(0)-3
        I3.append((g**2+1)/(k+3*(n-1)**2/((n-2)*(n-3))))
    I3 = np.concatenate(I3)
    print(f"   hitung ulang (semua titik, 4000 arah seragam): {I3.mean():.4f}")
