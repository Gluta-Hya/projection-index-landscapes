import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import csv, json, numpy as np
from scipy.signal import find_peaks
from gaya import *
import fase2_sistem as F, fase3c_gelombang as G
U = DATA
NPH = G.NPH
fase = np.arange(NPH)/NPH

def ukur(ic):
    dt = 0.005
    f = lambda t, s: [10*(s[1]-s[0]), s[0]*(28-s[2])-s[1], s[0]*s[1]-8/3*s[2]]
    simpan = F.N_SAMPEL; F.N_SAMPEL = 200_000
    S = F._flow(f, ic, dt, 30.0); F.N_SAMPEL = simpan
    x = S[0]; pm, _ = find_peaks(-np.abs(x)); bentuk, rasio = [], []
    for a, b in zip(pm[:-1], pm[1:]):
        seg = x[a:b+1]
        if not (np.all(seg > 0) or np.all(seg < 0)): continue
        g = np.interp(fase, np.linspace(0, 1, len(seg)), np.abs(seg)); c = g.mean(); amp = (g.max()-g.min())/2
        bentuk.append(np.roll((g-c)/amp, NPH//2)); rasio.append(amp/c)
    return np.array(bentuk), np.array(rasio)

rng = np.random.default_rng(7)
stat = []
for k in range(5):
    ic = [1.0, 1.0, 20.0] if k == 0 else list(rng.normal([1, 1, 20], 5))
    b, r = ukur(ic)
    Gm = G.normalisasi(b.mean(0)); H = np.abs(np.fft.rfft(Gm))[1:4]
    stat.append(dict(putaran=len(b), ac_median=float(np.median(r)), ac_q1=float(np.percentile(r, 25)), ac_q3=float(np.percentile(r, 75)),
                     h2=float(H[1]/H[0]), h3=float(H[2]/H[0]), korelasi=float(np.corrcoef(r, b.min(1))[0, 1]),
                     fraksi_D_buang=float(np.mean(1 + r*Gm.min() <= 0.02))))
    if k == 0: b0, r0 = b, r
ring = {key: (float(np.mean([s[key] for s in stat])), float(np.std([s[key] for s in stat]))) for key in stat[0]}
print(json.dumps(ring, indent=1)); json.dump(dict(per_lintasan=stat, ringkasan=ring), open("gambar_en/stat_waveform.json", "w"), indent=1)

r3c = list(csv.DictReader(open(U+'fase3c_ringkasan.csv', encoding='utf-8')))
r2 = list(csv.DictReader(open(U+'fase2_ringkasan.csv', encoding='utf-8')))
M = [2, 3, 4, 5, 6]
fig, axs = plt.subplots(1, 2, figsize=(LEBAR, 2.5), gridspec_kw=dict(width_ratios=[1.05, 1]))
ax = axs[0]
for g in b0[np.random.default_rng(1).choice(len(b0), 40, replace=False)]:
    ax.plot(fase, g, color="#d9d8d2", lw=0.5)
Gm0 = G.normalisasi(b0.mean(0))
ax.plot(fase, G.normalisasi(np.cos(2*np.pi*fase)), color=BIRU, lw=1.4, ls=(0, (4, 2)), label="cosine (model A)")
ax.plot(fase, G.potong_harmonik(Gm0, 2), color=ORANYE, lw=1.4, label="harmonics 1 and 2 (model B)")
ax.plot(fase, Gm0, color=TINTA, lw=1.6, label="Lorenz mean (model C)")
ax.set_xlabel(r"rotation phase ($|x|$ minimum at 0.5)"); ax.set_ylabel(r"$g = (|x| - c)/a$")
ax.set_title("(a) In-wing waveform", loc="left", fontweight="bold")
ax.set_ylim(-1.25, 1.95); ax.legend(loc="upper center", fontsize=6.3, handlelength=2.0, ncol=1)
ax = axs[1]
for jenis, lab, c, mk in (("A cos", "A, cosine", BIRU, "o"), ("B harmonik 1-2", "B, harmonics 1 and 2", ORANYE, "s"), ("E putaran nyata", "E, real rotations", AQUA, "^")):
    mu = [np.mean([float(d['c_m']) for d in r3c if d['jenis']==jenis and d['m']==str(m)]) for m in M]
    sd = [np.std([float(d['c_m']) for d in r3c if d['jenis']==jenis and d['m']==str(m)]) for m in M]
    ax.errorbar(M, mu, yerr=sd, color=c, marker=mk, ms=4.5, mec="white", mew=0.6, lw=1.3, capsize=2, label=lab)
mu = [np.mean([float(d['Wbintang_per_P']) for d in r2 if d['sistem']=='lorenz28' and d['m']==str(m)]) for m in M]
sd = [np.std([float(d['Wbintang_per_P']) for d in r2 if d['sistem']=='lorenz28' and d['m']==str(m)]) for m in M]
ax.errorbar(M, mu, yerr=sd, color=TINTA, marker="D", ms=4.2, mec="white", mew=0.6, lw=1.3, capsize=2, label=r"Lorenz $\rho=28$")
ax.set_xticks(M); ax.set_xlabel(r"$m$"); ax.set_ylabel(r"$c_m = W^*/P$")
ax.set_title("(b) Model ladder", loc="left", fontweight="bold"); ax.legend(loc="upper left", fontsize=6.3)
ax.set_ylim(0.44, 0.80)
fig.tight_layout(w_pad=1.2)
fig.savefig("gambar_en/fig4_ladder.png")
print("ok")
