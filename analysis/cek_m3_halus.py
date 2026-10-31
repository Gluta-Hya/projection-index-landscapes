import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import numpy as np
rng = np.random.default_rng(11)
def b(r): return 2*(2+r)**2/(3*r**2+24*r+8)
ok = True
for m in (3, 4):
    U = rng.standard_normal((40000, m)); U /= np.linalg.norm(U, axis=1, keepdims=True)
    C2 = U.sum(1)**2
    # simetri L(tau) = L(P - tau): cukup tau/P di (0, 1/2]; W/P = (m-1) tau/P
    tg = np.round(np.arange(0.002, 0.5001, 0.002), 4)
    E = np.exp(-1j*np.outer(2*np.pi*tg, np.arange(m)))      # (ntau, m)
    Z2 = np.abs(U @ E.T)**2                                   # (K, ntau)
    for ac in np.round(np.arange(0.60, 1.0001, 0.05), 2):
        L = np.mean(b(ac**2*Z2/C2[:, None]), axis=0)
        i = int(np.argmax(L))
        # bandingkan juga dengan limit tau -> 0 (rho = (a/c)^2 untuk semua arah)
        L0 = b(ac**2)
        global_ok = (tg[i] == 0.5) and (L[i] > L0)
        ok &= global_ok
        print(f"m={m} a/c={ac:.2f}: argmax tau/P={tg[i]:.3f} (W/P={(m-1)*tg[i]:.2f}), L={L[i]:.4f}, L(tau->0)={L0:.4f}  {'OK' if global_ok else 'TIDAK'}")
print("semua OK:", ok)
