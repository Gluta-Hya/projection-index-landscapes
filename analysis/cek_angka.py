import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import csv, numpy as np
U = DATA
r2 = list(csv.DictReader(open(U + "fase2_ringkasan.csv", encoding="utf-8")))
M = (2, 3, 4, 5, 6)
def cm(sistem, m):
    v = [float(r["Wbintang_per_P"]) for r in r2 if r["sistem"] == sistem and int(r["m"]) == m]
    return np.mean(v), np.std(v), len(v)
print("Tabel 1 (tak dibulatkan), selisih antar rho:")
for m in M:
    vals = [cm(s, m)[0] for s in ("lorenz28", "lorenz35", "lorenz45")]
    sds = [cm(s, m)[1] for s in ("lorenz28", "lorenz35", "lorenz45")]
    print(f"  m={m}: " + "  ".join(f"{v:.5f}" for v in vals) + f"   spread={max(vals)-min(vals):.5f}   sd={['%.4f'%s for s in sds]}")
P = {s: np.mean([float(r["P"]) for r in r2 if r["sistem"] == s]) for s in sorted({r["sistem"] for r in r2})}
Psd = {s: np.std([float(r["P"]) for r in r2 if r["sistem"] == s]) for s in sorted({r["sistem"] for r in r2})}
print("P:", {k: (round(v, 4), round(Psd[k], 4)) for k, v in P.items()})
for s in ("rossler",):
    for m in M:
        rows = [r for r in r2 if r["sistem"].startswith("ross") and int(r["m"]) == m]
        print("  rossler m", m, [(round(float(r["Wbintang_per_P"]), 3), r["di_tepi_grid"]) for r in rows])
print("sistem:", sorted({r["sistem"] for r in r2}))
# lembah m=2
for s in ("lorenz28", "lorenz35", "lorenz45"):
    v = [float(r["lembah_per_P"]) for r in r2 if r["sistem"] == s and int(r["m"]) == 2]
    w = [float(r["lembah_nilai"]) for r in r2 if r["sistem"] == s and int(r["m"]) == 2]
    print(f"  lembah {s}: {np.mean(v):.3f} ± {np.std(v):.3f}, nilai {np.mean(w):.3f}")

r3c = list(csv.DictReader(open(U + "fase3c_ringkasan.csv", encoding="utf-8")))
def c3c(j, m):
    v = [float(r["c_m"]) for r in r3c if r["jenis"] == j and int(r["m"]) == m]
    return np.mean(v), np.std(v)
J = ("A cos", "B harmonik 1-2", "C bentuk rata2", "D + a/c acak", "E putaran nyata")
L = {m: cm("lorenz28", m)[0] for m in M}
print("\nTabel 3 (tak dibulatkan):")
for j in J:
    v = [c3c(j, m)[0] for m in M]
    sd = [c3c(j, m)[1] for m in M]
    print(f"  {j:18s} " + " ".join(f"{x:.4f}" for x in v) + f"  c6-c3={v[4]-v[1]:+.4f}  maxsd={max(sd):.4f}")
print(f"  Lorenz             " + " ".join(f"{L[m]:.4f}" for m in M) + f"  c6-c3={L[6]-L[3]:+.4f}")
for j in ("B harmonik 1-2", "C bentuk rata2", "E putaran nyata"):
    fr = [(c3c(j, m)[0] - c3c("A cos", m)[0]) / (L[m] - c3c("A cos", m)[0]) for m in M]
    print(f"  tutup {j:18s}: " + " ".join(f"m{m}={f*100:.1f}%" for m, f in zip(M, fr)))
print("  sisa E-Lorenz:", " ".join(f"{L[m]-c3c('E putaran nyata', m)[0]:.4f}" for m in M))

r3b = list(csv.DictReader(open(DATA + "fase3b_ringkasan.csv", encoding="utf-8")))
print("\nFase 3b:")
for lab in sorted({(r["jenis"], r["label"]) for r in r3b}):
    rows = [r for r in r3b if (r["jenis"], r["label"]) == lab]
    out = []
    for m in M:
        v = [float(r["c_m"]) for r in rows if int(r["m"]) == m]
        t = [r["di_tepi_grid"] for r in rows if int(r["m"]) == m]
        out.append(f"m{m}={np.mean(v):.4f}{'*' if 'True' in t else ''}")
    D = np.mean([float(r["D_terukur"]) for r in rows])
    v3 = np.mean([float(r["c_m"]) for r in rows if int(r["m"]) == 3]); v6 = np.mean([float(r["c_m"]) for r in rows if int(r["m"]) == 6])
    sd = max(np.std([float(r["c_m"]) for r in rows if int(r["m"]) == m]) for m in M)
    print(f"  {lab[0]:8s} {lab[1]:10s} D={D:.2f} " + " ".join(out) + f"  c6-c3={v6-v3:+.4f} maxsd={sd:.4f}")

print("\nPresisi 5 desimal (pemeriksaan pembulatan tabel):")
for j in J:
    print(f"  {j:18s} " + " ".join(f"{c3c(j, m)[0]:.5f}" for m in M))
for lab in [("konstan","a/c=0.60"),("konstan","a/c=0.75"),("konstan","a/c=0.90"),("spiral","λ=1.25"),("acak","λ=1.25"),("spiral","λ=1.6"),("acak","λ=1.6"),("spiral","λ=2.0"),("acak","λ=2.0")]:
    rows = [r for r in r3b if (r["jenis"], r["label"]) == lab]
    v3 = np.mean([float(r["c_m"]) for r in rows if int(r["m"]) == 3]); v6 = np.mean([float(r["c_m"]) for r in rows if int(r["m"]) == 6])
    D = np.mean([float(r["D_terukur"]) for r in rows])
    print(f"  {lab[0]:8s} {lab[1]:10s} D={D:.4f} c3={v3:.5f} c6={v6:.5f} d={v6-v3:+.5f}")
for lam in ("λ=1.25","λ=1.6","λ=2.0"):
    g = lambda j, m: np.mean([float(r["c_m"]) for r in r3b if r["jenis"]==j and r["label"]==lam and int(r["m"])==m])
    print(f"  {lam}: spiral-acak c3={g('spiral',3)-g('acak',3):+.5f}  c6={g('spiral',6)-g('acak',6):+.5f}")
for s in ("lorenz28", "lorenz35", "lorenz45"):
    print(f"  {s} " + " ".join(f"{cm(s, m)[0]:.5f}±{cm(s, m)[1]:.5f}" for m in M))
