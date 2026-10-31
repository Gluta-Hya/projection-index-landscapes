"""
Eksperimen invariansi lanskap I(u) terhadap delay tau (dan dimensi m).

Pertanyaan: fitur lanskap mana yang milik attractor itu sendiri,
dan mana yang hanya artefak dari pilihan koordinat (tau, m)?

Karena koordinat bola berubah setiap kali tau berubah, lanskap tidak
dibandingkan titik-demi-titik. Yang dibandingkan adalah FITUR yang
tidak bergantung pada koordinat:
  - maks, median, rasio maks/median
  - fraksi luas bola di atas ambang (untuk I3: 0.555)
  - jumlah "pulau terang" (komponen terhubung dari 10% arah teratas)

Kontrol: surrogate acak-fase. Surrogate punya spektrum (autokorelasi)
yang sama dengan x(t) tapi struktur nonlinearnya dihancurkan. Fitur yang
sama antara data asli dan surrogate BUKAN milik attractor chaos.

Jalankan:  python invarians_lanskap.py
Kebutuhan: pip install numpy scipy matplotlib scikit-learn
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import skew, kurtosis
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import kneighbors_graph
import matplotlib.pyplot as plt
import csv

RNG = np.random.default_rng(42)
DT = 0.01


# ---------- 1. Data ----------
def lorenz_x(T=150.0, buang=20.0):
    f = lambda t, s: [10 * (s[1] - s[0]), s[0] * (28 - s[2]) - s[1], s[0] * s[1] - 8 / 3 * s[2]]
    t = np.arange(0, T, DT)
    sol = solve_ivp(f, (0, T), [1, 1, 1], t_eval=t, rtol=1e-9, atol=1e-9)
    return sol.y[0, int(buang / DT):]


def surrogate_acak_fase(x):
    X = np.fft.rfft(x - x.mean())
    fase = RNG.uniform(0, 2 * np.pi, len(X))
    fase[0] = 0
    return np.fft.irfft(np.abs(X) * np.exp(1j * fase), n=len(x)) + x.mean()


# ---------- 2. Pengangkatan phi (Takens) ----------
def embed(x, m, tau):
    N = len(x) - (m - 1) * tau
    V = np.stack([x[(m - 1 - k) * tau:(m - 1 - k) * tau + N] for k in range(m)], axis=1)
    return (V - V.mean(0)) / V.std(0)


# ---------- 3. Arah-arah di S^{m-1} (cukup setengah bola, karena I(u)=I(-u)) ----------
def arah(m, K):
    if m == 3:  # bola Fibonacci: sebaran hampir seragam
        i = np.arange(K) + 0.5
        z = 1 - 2 * i / (2 * K)          # hanya z>0 (setengah bola)
        r = np.sqrt(1 - z ** 2)
        a = np.pi * (1 + 5 ** 0.5) * i
        return np.stack([r * np.cos(a), r * np.sin(a), z], 1)
    U = RNG.standard_normal((K, m))       # m>3: acak seragam
    U /= np.linalg.norm(U, axis=1, keepdims=True)
    U[U[:, -1] < 0] *= -1
    return U


# ---------- 4. Indeks pola ----------
G = lambda y: np.log(np.cosh(y))
EGN = G(RNG.standard_normal(400_000)).mean()
H_GAUSS = 0.5 * np.log(2 * np.pi * np.e)
EDGES = np.linspace(-5, 5, 201)
W = EDGES[1] - EDGES[0]


def indeks(y):
    y = (y - y.mean()) / y.std()
    n = len(y)
    I1 = (G(y).mean() - EGN) ** 2
    p, _ = np.histogram(y, EDGES, density=True)
    q = p[p > 0]
    I2 = H_GAUSS + (q * np.log(q)).sum() * W
    I3 = (skew(y) ** 2 + 1) / (kurtosis(y) + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    return I1, I2, I3


def lanskap(V, U, sub=4):
    Y = V[::sub] @ U.T
    return np.array([indeks(Y[:, j]) for j in range(U.shape[0])])  # (K, 3)


# ---------- 5. Fitur yang bebas koordinat ----------
def jumlah_pulau(U, nilai, q=0.90, k=8):
    """Komponen terhubung dari arah dengan nilai di atas kuantil q.
    Antipodal digabung: u dan -u dianggap tetangga."""
    atas = nilai >= np.quantile(nilai, q)
    P = U[atas]
    if len(P) < 2:
        return len(P)
    P2 = np.vstack([P, -P])
    A = kneighbors_graph(P2, n_neighbors=min(k, len(P2) - 1), include_self=False)
    n_comp, lab = connected_components(A + A.T, directed=False)
    # pasangkan u dengan -u
    pasangan = {tuple(sorted((lab[i], lab[i + len(P)]))) for i in range(len(P))}
    induk = {c: c for c in range(n_comp)}
    def cari(c):
        while induk[c] != c:
            c = induk[c]
        return c
    for a, b in pasangan:
        induk[cari(a)] = cari(b)
    return len({cari(c) for c in range(n_comp)})


def fitur(U, L):
    out = {}
    for j, nama in enumerate(["I1", "I2", "I3"]):
        v = L[:, j]
        out[f"{nama}_maks"] = v.max()
        out[f"{nama}_median"] = np.median(v)
        out[f"{nama}_rasio"] = v.max() / (np.median(v) + 1e-12)
        out[f"{nama}_pulau"] = jumlah_pulau(U, v)
    out["I3_fraksi_bimodal"] = (L[:, 2] > 0.555).mean()
    return out


# ---------- 6. Eksperimen ----------
def jalankan(m_list=(3, 4, 5), tau_list=range(4, 33, 2), K=1500):
    x = lorenz_x()
    sumber = {"lorenz": x, "surrogate": surrogate_acak_fase(x)}
    baris = []
    for m in m_list:
        U = arah(m, K)
        for tau in tau_list:
            for nama, s in sumber.items():
                f = fitur(U, lanskap(embed(s, m, tau), U))
                f.update(sistem=nama, m=m, tau=tau, tau_detik=tau * DT)
                baris.append(f)
                print(f"m={m} tau={tau:2d} {nama:9s} I3_bimodal={f['I3_fraksi_bimodal']:.2f} "
                      f"pulau(I2)={f['I2_pulau']}")
    with open("hasil_invarians.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(baris[0].keys()))
        w.writeheader()
        w.writerows(baris)
    return baris


def gambar(baris, m_list):
    kolom = ["I1_rasio", "I2_maks", "I3_fraksi_bimodal", "I2_pulau"]
    fig, axs = plt.subplots(len(m_list), len(kolom), figsize=(16, 3.6 * len(m_list)), squeeze=False)
    for r, m in enumerate(m_list):
        for c, k in enumerate(kolom):
            ax = axs[r, c]
            for nama, gaya in [("lorenz", "-o"), ("surrogate", "--s")]:
                d = [b for b in baris if b["m"] == m and b["sistem"] == nama]
                ax.plot([b["tau_detik"] for b in d], [b[k] for b in d], gaya, ms=3, label=nama)
            ax.set_title(f"m={m}: {k}")
            ax.set_xlabel(r"$\tau$ (detik)")
            if r == 0 and c == 0:
                ax.legend()
    plt.tight_layout()
    plt.savefig("hasil_invarians.png", dpi=120)
    print("Tersimpan: hasil_invarians.csv dan hasil_invarians.png")


if __name__ == "__main__":
    M = (3, 4, 5)
    gambar(jalankan(m_list=M), M)
