import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import csv, numpy as np, json
U = DATA
rd = lambda f: list(csv.DictReader(open(U+f, encoding='utf-8')))
out = {}
# ---------- Fase 1 (Larry) ----------
r1 = rd('fase1_ringkasan.csv')
print("FASE 1: W* (s) per m, T")
for T in ('75.0','150.0'):
    for m in '234567':
        v = np.array([float(r['W_bintang']) for r in r1 if r['T']==T and r['m']==m])
        print(f"  T={T} m={m}: {v.mean():.3f} ± {v.std():.3f}  (n={len(v)})")
d = []
for m in '234567':
    a = np.array([float(r['W_bintang']) for r in r1 if r['T']=='75.0' and r['m']==m])
    b = np.array([float(r['W_bintang']) for r in r1 if r['T']=='150.0' and r['m']==m])
    d.append(abs(a.mean()-b.mean()))
print("  maks |W*(75)-W*(150)| =", max(d))
h1 = rd('fase1_hasil.csv')
sur = np.array([float(r['I3_rata']) for r in h1 if r['sistem']=='surrogate'])
print("  surrogate I3_rata fase1: mean %.4f min %.4f max %.4f" % (sur.mean(), sur.min(), sur.max()))
# ---------- Fase 2 (Larry) ----------
r2 = rd('fase2_ringkasan.csv'); h2 = rd('fase2_hasil.csv')
print("\nFASE 2")
for s in ('lorenz28','lorenz35','lorenz45','rossler','henon'):
    P = np.array(sorted({(r['realisasi'], float(r['P'])) for r in r2 if r['sistem']==s}))[:,1].astype(float)
    print(f" {s}: P = {P.mean():.4f} ± {P.std():.4f}")
    for m in '23456':
        v = np.array([float(r['Wbintang_per_P']) for r in r2 if r['sistem']==s and r['m']==m])
        pk = np.array([float(r['puncak']) for r in r2 if r['sistem']==s and r['m']==m])
        print(f"   m={m}: c_m = {v.mean():.3f} ± {v.std():.3f}, puncak = {pk.mean():.3f}")
    lem = np.array([float(r['lembah_per_P']) for r in r2 if r['sistem']==s and r['m']=='2' and r['lembah_per_P']!='nan'])
    if len(lem): print(f"   lembah m=2: {lem.mean():.3f} ± {lem.std():.3f}; nilai {np.mean([float(r['lembah_nilai']) for r in r2 if r['sistem']==s and r['m']=='2']):.3f}")
    su = np.array([float(r['I3_rata']) for r in h2 if r['sistem']==s and r['sumber']=='surrogate'])
    print(f"   surrogate: {su.mean():.3f} (rentang {su.min():.3f} sampai {su.max():.3f})")
# lintas rho: rentang c_m per m
print(" rentang c_m lintas rho per m:")
for m in '23456':
    mus = [np.mean([float(r['Wbintang_per_P']) for r in r2 if r['sistem']==s and r['m']==m]) for s in ('lorenz28','lorenz35','lorenz45')]
    print(f"   m={m}: {[round(x,3) for x in mus]} rentang {max(mus)-min(mus):.3f}")
# W* dalam detik untuk m=3 rho=28 dan MI rule
P28 = np.mean([float(r['P']) for r in r2 if r['sistem']=='lorenz28'])
c3 = np.mean([float(r['Wbintang_per_P']) for r in r2 if r['sistem']=='lorenz28' and r['m']=='3'])
print(f" rho=28: W*(m=3) = {c3*P28:.3f} s; aturan MI tau=0.164 s -> W(m=3) = {2*0.164:.3f} s")
# ---------- Fase 3 (Larry) ----------
h3 = rd('fase3_hasil.csv')
s = np.array([float(r['I3_rata']) for r in h3 if r['model']=='sinus'])
print(f"\nFASE 3: sinus I3 = {s.mean():.4f} (rentang {s.min():.4f} sampai {s.max():.4f})")
for mod in ('dua_sayap_jarang','dua_sayap_lorenz'):
    rr = [r for r in h3 if r['model']==mod and r['m']=='2']
    W = np.array([float(r['W_per_P']) for r in rr]); y = np.array([float(r['I3_rata']) for r in rr])
    print(f"  {mod} m=2: maks {y.max():.4f} di W/P={W[y.argmax()]:.2f}; nilai di W=P {y[np.argmin(abs(W-1))]:.4f}")
for m in '23456':
    rr = [r for r in h3 if r['model']=='dua_sayap_lorenz' and r['m']==m]
    W = np.array([float(r['W_per_P']) for r in rr]); y = np.array([float(r['I3_rata']) for r in rr])
    print(f"  model D=2.2 m={m}: W*/P (grid) = {W[y.argmax()]:.2f}")
# ---------- Fase 3b (run saya, identik dengan Larry) ----------
r3b = list(csv.DictReader(open(DATA + 'fase3b_ringkasan.csv', encoding='utf-8')))
print("\nFASE 3b")
labs = sorted({(r['jenis'], r['label']) for r in r3b}, key=lambda t: (t[0]!='konstan', t[1], t[0]))
for j, l in labs:
    st = [np.mean([float(r['c_m']) for r in r3b if r['jenis']==j and r['label']==l and r['m']==m]) for m in '23456']
    D = np.mean([float(r['D_terukur']) for r in r3b if r['jenis']==j and r['label']==l])
    tepi = any(r['di_tepi_grid']=='True' for r in r3b if r['jenis']==j and r['label']==l)
    print(f"  {j:8s} {l:10s} D={D:.2f}: " + " ".join(f"{x:.3f}" for x in st) + f"  kemiringan {st[4]-st[1]:+.3f}" + ("  TEPI" if tepi else ""))
# ---------- Fase 3c (Larry) ----------
r3c = rd('fase3c_ringkasan.csv')
print("\nFASE 3c (Larry)")
for j in ('A cos','B harmonik 1-2','C bentuk rata2','D + a/c acak','E putaran nyata'):
    st = [(np.mean(v), np.std(v)) for v in ([float(r['c_m']) for r in r3c if r['jenis']==j and r['m']==m] for m in '23456')]
    print(f"  {j:16s} " + " ".join(f"{a:.3f}±{b:.3f}" for a,b in st) + f"  kemiringan {st[4][0]-st[1][0]:+.3f}")
lor = [np.mean([float(r['Wbintang_per_P']) for r in r2 if r['sistem']=='lorenz28' and r['m']==m]) for m in '23456']
A = [np.mean([float(r['c_m']) for r in r3c if r['jenis']=='A cos' and r['m']==m]) for m in '23456']
for j in ('B harmonik 1-2','C bentuk rata2','E putaran nyata'):
    X = [np.mean([float(r['c_m']) for r in r3c if r['jenis']==j and r['m']==m]) for m in '23456']
    print(f"  fraksi selisih tertutup {j}: " + " ".join(f"m{m}:{(X[i]-A[i])/(lor[i]-A[i]):.2f}" for i,m in enumerate('23456') if m!='2'))
