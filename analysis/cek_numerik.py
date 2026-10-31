import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import numpy as np, sys
from scipy.optimize import brentq
f = lambda r: 2*(2+r)**2/(3*r*r+24*r+8)
al = np.linspace(0, np.pi, 200000, endpoint=False); s2 = np.sin(2*al); ok = (1+s2) > 1e-12
def L(ac, W):  # m=2 tanpa perpindahan
    rho = ac**2*(1+s2[ok]*np.cos(2*np.pi*W))/(1+s2[ok]); return f(rho).mean()
Ws = np.linspace(0, 0.5, 501)
print("a/c | argmax W/P di [0,0.5] | L(0) | L(P/2) | min interior di W/P")
for ac in (0.3, 0.4, 0.45, 0.457, 0.46, 0.5, 0.6, 0.8, 1.0, 1.63, 2.0, 3.0, 5.0):
    v = np.array([L(ac, w) for w in Ws]); i = v.argmax(); j = v.argmin()
    # local maxima interior
    locmax = [Ws[k] for k in range(1, len(Ws)-1) if v[k] > v[k-1] and v[k] > v[k+1]]
    print(f"{ac:5.3f} | {Ws[i]:.3f} | {v[0]:.4f} | {v[-1]:.4f} | {Ws[j]:.3f} | lokal maks interior: {locmax}")
g = lambda ac: L(ac, 0.5) - L(ac, 0.0)
print("ambang a/c (L(P/2)=L(0)):", brentq(g, 0.4, 0.5, xtol=1e-6))
# apakah ada ambang atas?
for ac in (1.63, 3, 10, 30):
    print(ac, 'L(P/2)-L(0) =', round(g(ac), 4))
