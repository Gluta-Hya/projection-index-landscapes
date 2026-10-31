import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import csv, numpy as np
from gaya import *
U = DATA
h = list(csv.DictReader(open(U+'fase2_hasil.csv', encoding='utf-8')))
r = list(csv.DictReader(open(U+'fase2_ringkasan.csv', encoding='utf-8')))
SIS = [("lorenz28", r"$\rho=28$", BIRU, "o"), ("lorenz35", r"$\rho=35$", ORANYE, "s"), ("lorenz45", r"$\rho=45$", AQUA, "^")]

def kurva(sistem, m, sumber):
    d = [b for b in h if b['sistem']==sistem and b['m']==str(m) and b['sumber']==sumber]
    taus = sorted({int(b['tau']) for b in d})
    x = np.array([np.mean([float(b['W_per_P']) for b in d if int(b['tau'])==t]) for t in taus])
    y = np.array([[float(b['I3_rata']) for b in d if int(b['tau'])==t] for t in taus])
    return x, y.mean(1), y.std(1)

fig, axs = plt.subplots(1, 3, figsize=(LEBAR, 2.35), gridspec_kw=dict(width_ratios=[1, 1, 0.9]))
for ax, m, huruf in ((axs[0], 2, "(a)"), (axs[1], 4, "(b)")):
    for s, lab, c, mk in SIS:
        x, mu, sd = kurva(s, m, 'asli')
        ax.fill_between(x, mu-sd, mu+sd, color=c, alpha=0.15, lw=0)
        ax.plot(x, mu, color=c, lw=1.3, label=lab)
    xs, mus, _ = kurva('lorenz28', m, 'surrogate')
    ax.plot(xs, mus, color=REDUP, lw=1.0, ls=(0, (3, 2)), label="surrogate")
    ax.axhline(1/3, color=SUMBU, lw=0.7)
    for v, t in ((0.5, "$P/2$"), (1.0, "$P$")):
        ax.axvline(v, color=SUMBU, lw=0.7)
        ax.text(v + 0.03, 0.592, t, ha="left", va="top", fontsize=7, color=TINTA2)
    ax.set_xlim(0, 2.0); ax.set_ylim(0.30, 0.60)
    ax.set_xlabel(r"$W/P$"); ax.set_title(f"{huruf} $m={m}$", loc="left", fontweight="bold")
axs[0].set_ylabel(r"mean $I_3$ over directions")
axs[1].legend(loc="lower center", ncol=2, handlelength=1.6, columnspacing=1.0)
ax = axs[2]
M = [2, 3, 4, 5, 6]
for k, (s, lab, c, mk) in enumerate(SIS):
    mu = [np.mean([float(b['Wbintang_per_P']) for b in r if b['sistem']==s and b['m']==str(m)]) for m in M]
    sd = [np.std([float(b['Wbintang_per_P']) for b in r if b['sistem']==s and b['m']==str(m)]) for m in M]
    ax.errorbar(np.array(M)+(k-1)*0.09, mu, yerr=sd, color=c, marker=mk, ms=4.5, mec="white", mew=0.6,
                lw=1.2, capsize=2, label=lab)
ax.set_xticks(M); ax.set_xlabel(r"$m$"); ax.set_ylabel(r"$c_m = W^*/P$")
ax.set_title("(c) Peak location", loc="left", fontweight="bold"); ax.set_ylim(0.44, 0.78)
ax.legend(loc="lower right", handlelength=1.6)
fig.tight_layout(w_pad=1.0)
fig.savefig("gambar_en/fig2_scaling.png")
print("ok")
