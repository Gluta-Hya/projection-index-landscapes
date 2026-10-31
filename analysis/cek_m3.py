import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import numpy as np
rng = np.random.default_rng(5)
def b(r): return 2*(2+r)**2/(3*r**2+24*r+8)
for m in (3, 4):
    U = rng.standard_normal((40000, m)); U /= np.linalg.norm(U, axis=1, keepdims=True)
    Wg = np.round(np.arange(0.02, 2.0001, 0.01), 3)
    for ac in (0.3, 0.45, 0.6, 0.75, 0.9, 1.0):
        L = []
        for W in Wg:
            psi = 2*np.pi*(W/(m-1))          # tau/P = W/((m-1)P)
            k = np.arange(m)
            Z = U @ np.exp(-1j*k*psi)
            C = U.sum(1)
            rho = ac**2*np.abs(Z)**2/C**2
            L.append(np.mean(b(rho)))
        L = np.array(L); i = np.argmax(L)
        # local maxima in interior
        loc = [Wg[j] for j in range(1, len(L)-1) if L[j] > L[j-1] and L[j] >= L[j+1]]
        print(f"m={m} a/c={ac:.2f}: argmax W/P={Wg[i]:.2f} (L={L[i]:.4f}); L(0.02)={L[0]:.4f}; maksimum lokal di {loc}")
