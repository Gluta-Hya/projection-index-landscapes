import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import numpy as np, json
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr
from gaya import *

dt = 0.01
f = lambda t, s: [10*(s[1]-s[0]), s[0]*(28-s[2])-s[1], s[0]*s[1]-8/3*s[2]]
t = np.arange(0, 120, dt)
S = solve_ivp(f, (0, 120), [1, 1, 1], t_eval=t, rtol=1e-9, atol=1e-9).y[:, 2000:]
x = S[0]
tau, m = 12, 3
N = len(x) - (m-1)*tau
V = np.stack([x[(m-1-k)*tau:(m-1-k)*tau+N] for k in range(m)], 1)
V = (V - V.mean(0))/V.std(0)
sub = V[::5]; n = len(sub)
th = np.linspace(0, np.pi, 91); ph = np.linspace(0, 2*np.pi, 181)
TT, PP = np.meshgrid(th, ph, indexing="ij")
U = np.stack([np.sin(TT)*np.cos(PP), np.sin(TT)*np.sin(PP), np.cos(TT)], -1).reshape(-1, 3)
Y = sub @ U.T; Y = (Y - Y.mean(0))/Y.std(0); Y2 = Y*Y
g = (Y2*Y).mean(0); k = (Y2*Y2).mean(0) - 3
I3 = ((g**2 + 1)/(k + 3*(n-1)**2/((n-2)*(n-3)))).reshape(TT.shape)
EG = np.mean(np.log(np.cosh(np.random.default_rng(0).standard_normal(400000))))
I1 = ((np.log(np.cosh(Y)).mean(0) - EG)**2).reshape(TT.shape)
w = np.sin(TT)                                   # bobot luas
frac = float((w*(I3 > 5/9)).sum()/w.sum())
rho_s = float(spearmanr(I1.ravel(), I3.ravel())[0])
i3 = np.unravel_index(I3.argmax(), I3.shape); i1 = np.unravel_index(I1.argmax(), I1.shape); j3 = np.unravel_index(I3.argmin(), I3.shape)
u3 = U.reshape(*TT.shape, 3)[i3]; u1 = U.reshape(*TT.shape, 3)[i1]; v3 = U.reshape(*TT.shape, 3)[j3]
sudut = float(np.degrees(np.arccos(min(1, abs(u3 @ u1)))))
stat = dict(I3_maks=float(I3.max()), I3_min=float(I3.min()), I3_rata_luas=float((w*I3).sum()/w.sum()),
            fraksi_bimodal=frac, spearman_I1_I3=rho_s, sudut_argmax_I1_I3_derajat=sudut,
            arah_maks_theta_phi=[float(np.degrees(th[i3[0]])), float(np.degrees(ph[i3[1]]))])
print(json.dumps(stat, indent=1))
json.dump(stat, open("gambar/stat_gambar1.json", "w"), indent=1)

fig = plt.figure(figsize=(LEBAR, 4.3))
ax = fig.add_subplot(2, 3, 1, projection="3d")
ax.plot(*S[:, :6000], lw=0.25, color=BIRU)
ax.set_title("(a) Attractor Lorenz $(x, y, z)$", loc="left", pad=0, fontweight="bold")
ax.set_xlabel("$x$", labelpad=-14); ax.set_ylabel("$y$", labelpad=-14); ax.set_zlabel("$z$", labelpad=-14)
ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
for a in (ax.xaxis, ax.yaxis, ax.zaxis): a.pane.set_facecolor("white"); a.pane.set_edgecolor(GRID)

ax = fig.add_subplot(2, 3, 2)
tt = np.arange(3000)*dt
ax.plot(tt, x[:3000], lw=0.7, color=BIRU)
ax.set_xlabel("$t$"); ax.set_ylabel("$x(t)$"); ax.set_title("(b) Deret yang diamati, $x(t)$", loc="left", fontweight="bold")

ax = fig.add_subplot(2, 3, 3, projection="3d")
ax.plot(V[:6000, 0], V[:6000, 1], V[:6000, 2], lw=0.25, color=BIRU)
ax.set_title("(c) Rekonstruksi dari $x(t)$", loc="left", pad=0, fontweight="bold")
ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
for a in (ax.xaxis, ax.yaxis, ax.zaxis): a.pane.set_facecolor("white"); a.pane.set_edgecolor(GRID)

ax = fig.add_subplot(2, 3, 4)
im = ax.imshow(I3, extent=[0, 360, 180, 0], aspect="auto", cmap=CMAP_BIRU, interpolation="nearest")
ax.contour(np.degrees(ph), np.degrees(th), I3, [5/9], colors=TINTA, linewidths=0.6)
ax.plot(np.degrees(ph[i3[1]]), np.degrees(th[i3[0]]), marker="*", ms=9, color=ORANYE, mec="white", mew=0.8, ls="none")
ax.plot(np.degrees(ph[j3[1]]), np.degrees(th[j3[0]]), marker="o", ms=5, color="white", mec=TINTA, mew=0.8, ls="none")
ax.set_xlabel(r"$\varphi$ (derajat)"); ax.set_ylabel(r"$\theta$ (derajat)"); ax.grid(False)
ax.set_title(r"(d) Lanskap $I_3(u)$ di $S^2$", loc="left", fontweight="bold")
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6.5)

gx = np.linspace(-4, 4, 300); gauss = np.exp(-gx**2/2)/np.sqrt(2*np.pi)
for idx, (u, warna, judul, huruf) in enumerate(((u3, BIRU, r"(e) Arah maksimum $I_3$", "(e)"), (v3, REDUP, r"(f) Arah minimum $I_3$", "(f)"))):
    ax = fig.add_subplot(2, 3, 5+idx)
    y = V @ u; y = (y - y.mean())/y.std()
    ax.hist(y, bins=90, density=True, color=warna, alpha=0.55, edgecolor="none")
    ax.plot(gx, gauss, color=TINTA, lw=0.9, ls="--", label="Gaussian")
    ax.set_xlabel(r"$u^\top v$ (distandarkan)"); ax.set_title(judul, loc="left", fontweight="bold"); ax.set_xlim(-4, 4)
    if idx == 0: ax.set_ylabel("kepadatan"); ax.legend(loc="upper right")
fig.tight_layout(w_pad=1.2, h_pad=1.0)
koma_semua(fig); fig.savefig("gambar/gambar1_lanskap.png")
print("tersimpan")
