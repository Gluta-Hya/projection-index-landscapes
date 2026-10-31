import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import csv, numpy as np
from gaya import *
U = DATA
f = lambda r: 2*(2+r)**2/(3*r*r+24*r+8)
al = np.linspace(0, np.pi, 100000, endpoint=False); s2 = np.sin(2*al); ok = (1+s2) > 1e-12
def L(ac, W):
    rho = ac**2*(1+s2[ok]*np.cos(2*np.pi*W))/(1+s2[ok]); return f(rho).mean()

fig, axs = plt.subplots(1, 3, figsize=(LEBAR, 2.45), gridspec_kw=dict(width_ratios=[1, 1, 1.25]))
# (a) b(rho)
ax = axs[0]
r = np.logspace(-2, 2, 400)
ax.plot(r, f(r), color=BIRU, lw=1.5)
for v, t, xpos, ha in ((1, "dua titik 1", 1.1e-2, "left"), (2/3, "sinus 2/3", 95, "right"),
                       (5/9, "seragam 5/9", 1.1e-2, "left"), (1/3, "Gaussian 1/3", 1.1e-2, "left")):
    ax.axhline(v, color=SUMBU, lw=0.6)
    if t.startswith("seragam"):
        ax.text(xpos, v - 0.01, t, fontsize=6.2, color=TINTA2, va="top", ha=ha)
    else:
        ax.text(xpos, v + 0.01, t, fontsize=6.2, color=TINTA2, va="bottom", ha=ha)
ax.plot([8/3], [7/15], marker="o", ms=5, color=BIRU, mec="white", mew=0.8, ls="none")
ax.annotate("min $7/15$\ndi $\\rho=8/3$", xy=(8/3, 7/15), xytext=(9, 0.39), fontsize=6.2, color=TINTA2, ha="left",
            arrowprops=dict(arrowstyle="-", color=REDUP, lw=0.6))
ax.set_xscale("log"); ax.set_ylim(0.3, 1.06); ax.set_xlabel(r"$\rho = R^2/C^2$"); ax.set_ylabel(r"$b(\rho)$")
ax.set_title(r"(a) Fungsi $b(\rho)$", loc="left", fontweight="bold")
# (b) ambang a/c
ax = axs[1]
ac = np.linspace(0.05, 1.0, 120)
LP2 = np.array([L(a, 0.5) for a in ac]); LP = f(ac**2)
ax.plot(ac, LP2, color=BIRU, lw=1.5, label=r"$W = P/2$")
ax.plot(ac, LP, color=ORANYE, lw=1.5, ls=(0, (4, 2)), label=r"$W \to 0$ dan $W = P$")
ax.axvline(0.4573, color=SUMBU, lw=0.7)
ax.text(0.47, 0.985, r"$a/c = 0{,}457$", fontsize=6.3, color=TINTA2, va="top")
ax.set_xlabel(r"$a/c$"); ax.set_ylabel(r"rata-rata $b$ atas arah ($m=2$)")
ax.set_title("(b) Ambang puncak di $P/2$", loc="left", fontweight="bold"); ax.set_ylim(0.45, 1.02)
ax.legend(loc="lower left", handlelength=2.2, fontsize=6.3)
# (c) validasi
ax = axs[2]
h3 = list(csv.DictReader(open(U+'fase3_hasil.csv', encoding='utf-8')))
Wg = np.linspace(0, 2, 401)
ax.plot(Wg, [L(0.6, w) for w in Wg], color=TINTA, lw=1.0, label="rumus tanpa pindah")
for mod, lab, c, pred in (("dua_sayap_jarang", r"sim. $D=40$", BIRU, 0.5131), ("dua_sayap_lorenz", r"sim. $D=2{,}2$", ORANYE, 0.4629)):
    rr = [b for b in h3 if b['model']==mod and b['m']=='2']
    ax.plot([float(b['W_per_P']) for b in rr], [float(b['I3_rata']) for b in rr], color=c, lw=1.5, label=lab)
    ax.plot([1.0], [pred], marker="D", ms=5, mfc="white", mec=c, mew=1.2, ls="none")
ax.plot([], [], marker="D", ms=5, mfc="white", mec=TINTA2, mew=1.2, ls="none", label="Proposisi 3")
ax.set_xlabel(r"$W/P$"); ax.set_ylabel(r"rata-rata $I_3$ ($m=2$, $a/c=0{,}6$)")
ax.set_title("(c) Rumus dan simulasi", loc="left", fontweight="bold"); ax.set_ylim(0.43, 0.84); ax.set_xlim(0, 2)
ax.legend(loc="upper center", ncol=2, handlelength=1.4, fontsize=6.0, borderaxespad=0.15, columnspacing=0.6, handletextpad=0.4)
fig.tight_layout(w_pad=0.9)
koma_semua(fig); fig.savefig("gambar/gambar3_analitik.png")
print("ok")
